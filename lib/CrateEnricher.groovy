import groovy.json.JsonOutput
import groovy.json.JsonSlurper

import java.security.MessageDigest

/**
 * Completes the Workflow Run RO-Crate that nf-prov writes, from the run
 * manifest and the published files.
 *
 * Why this exists, and why it is Groovy in lib/ rather than a script the user
 * runs:
 *
 * nf-prov writes ro-crate-metadata.json roughly 1.6 seconds AFTER the entry
 * workflow's onComplete: handler returns, measured on a clean run, so no
 * Nextflow hook ever sees the file. The only thing that runs later than every
 * hook is a JVM shutdown hook, which main.nf registers from onComplete:. That
 * hook calls enrich() here.
 *
 * It is Groovy rather than Python because the host running Nextflow is only
 * guaranteed to have Java. Every tool the pipeline uses lives in a container;
 * adding a python3 dependency on the launching machine to finish a metadata
 * file would be a poor trade.
 *
 * What it adds, none of which nf-prov can express:
 *
 *   - sha256, contentSize and dateModified on every published file. The plugin
 *     records none, leaving a crate that names files a reader cannot verify.
 *   - hasPart on directory-valued outputs (deseq2_output and friends), which
 *     the plugin leaves as a bare name, so the result tables are described.
 *   - every contributor as a Person behind Role entities carrying CRediT URIs.
 *     The plugin's config holds one agent and no roles.
 *   - the project as a ResearchProject with its RAiD, the funder, the grant and
 *     the rights holder.
 *
 * Two term choices the RO-Crate 1.1 context forces. It has no checksum term and
 * schema.org has no checksum property, so @context is extended with one
 * pointing at SPDX; the alternative, a PropertyValue entity per file, costs
 * several hundred entities to say the same thing. `rightsHolder` is not in the
 * context either, so `copyrightHolder` is used.
 *
 * Idempotent: a second pass over an already-enriched crate changes nothing.
 */
class CrateEnricher {

    static final String SPDX_CHECKSUM = 'http://spdx.org/rdf/terms#checksumValue'
    static final String CREDIT = 'https://credit.niso.org/contributor-roles/'
    static final String CRATE = 'ro-crate-metadata.json'

    /**
     * @return a one-line summary for the console, or null when there was no
     *         crate to enrich (the plugin disabled, or it failed to write one).
     */
    static String enrich(String outdir, String manifestPath) {
        def crateFile = new File(outdir, CRATE)
        def manifestFile = new File(manifestPath)
        if (!crateFile.exists() || !manifestFile.exists()) return null

        def doc = new JsonSlurper().parse(crateFile)
        def manifest = new JsonSlurper().parse(manifestFile)
        def graph = doc['@graph']
        def byId = [:]
        graph.each { e -> byId[e['@id']] = e }

        extendContext(doc, 'sha256', SPDX_CHECKSUM)

        def root = byId['./']
        def attribution = manifest['attribution'] ?: [:]
        if (root != null) {
            // Licence and citation come from the manifest so the crate has one
            // source rather than repeating them in the plugin's own config.
            if (attribution['license'])
                root['license'] = 'https://spdx.org/licenses/' + attribution['license']
            if (attribution['doi'])
                root['identifier'] = 'https://doi.org/' + attribution['doi']
        }

        int hashed = addIntegrity(graph, outdir)
        def expanded = expandDirectories(graph, byId, outdir)
        int people = addPeople(graph, byId, root, attribution)
        boolean project = addProject(graph, byId, root, manifest['project'])

        crateFile.setText(JsonOutput.prettyPrint(JsonOutput.toJson(doc)) + '\n', 'UTF-8')

        return "RO-Crate: ${graph.size()} entities, ${hashed} files hashed, " +
               "${expanded[0]} result folders expanded (${expanded[1]} files), " +
               "${people} contributor(s)${project ? ', project' : ''}"
    }

    // ---- crate helpers -----------------------------------------------------

    private static void extendContext(doc, String term, String uri) {
        def ctx = doc['@context']
        if (ctx instanceof String) ctx = [ctx]
        def extra = ctx.find { it instanceof Map }
        if (extra == null) {
            extra = [:]
            ctx << extra
        }
        extra[term] = uri
        doc['@context'] = ctx
    }

    private static Map entity(graph, Map byId, String id, String type) {
        def e = byId[id]
        if (e == null) {
            e = ['@id': id]
            if (type) e['@type'] = type
            graph << e
            byId[id] = e
        }
        else if (type && !e.containsKey('@type')) {
            e['@type'] = type
        }
        return e
    }

    /** Add a reference without duplicating it, keeping a single value scalar. */
    private static void link(Map entity, String prop, String targetId) {
        def ref = ['@id': targetId]
        def cur = entity[prop]
        if (cur == null) {
            entity[prop] = ref
        }
        else if (cur instanceof Map) {
            if (cur['@id'] != targetId) entity[prop] = [cur, ref]
        }
        else if (cur instanceof List) {
            if (!cur.any { it instanceof Map && it['@id'] == targetId }) cur << ref
        }
    }

    private static void unlink(Map entity, String prop, String targetId) {
        def cur = entity[prop]
        if (cur instanceof Map && cur['@id'] == targetId) {
            entity.remove(prop)
        }
        else if (cur instanceof List) {
            def kept = cur.findAll { !(it instanceof Map && it['@id'] == targetId) }
            if (kept.size() == 1) entity[prop] = kept[0]
            else if (kept) entity[prop] = kept
            else entity.remove(prop)
        }
    }

    private static boolean isLocal(String id) {
        return !(id.startsWith('http://') || id.startsWith('https://')
                 || id.startsWith('file:') || id.startsWith('#'))
    }

    private static List typesOf(e) {
        def t = e['@type']
        return (t instanceof List) ? t : [t]
    }

    // ---- integrity ---------------------------------------------------------

    static String sha256(File f) {
        def md = MessageDigest.getInstance('SHA-256')
        f.withInputStream { ins ->
            byte[] buf = new byte[1 << 20]
            int n
            while ((n = ins.read(buf)) != -1) md.update(buf, 0, n)
        }
        return md.digest().encodeHex().toString()
    }

    private static void integrity(Map e, File f) {
        e['sha256'] = sha256(f)
        e['contentSize'] = String.valueOf(f.length())
        e['dateModified'] = new Date(f.lastModified())
            .format("yyyy-MM-dd'T'HH:mm:ss'Z'", TimeZone.getTimeZone('UTC'))
    }

    private static int addIntegrity(graph, String outdir) {
        int n = 0
        graph.toList().each { e ->
            if (!typesOf(e).contains('File')) return
            String id = e['@id']
            if (!isLocal(id)) return
            def f = new File(outdir, id)
            if (f.isFile()) {
                integrity(e, f)
                n++
            }
        }
        return n
    }

    /**
     * Directory-valued outputs arrive as a Dataset carrying only a name, so the
     * result tables inside them are undescribed. Give each one hasPart and a
     * File entity per file.
     */
    private static List expandDirectories(graph, Map byId, String outdir) {
        int dirs = 0, files = 0
        graph.toList().each { e ->
            if (e['@type'] != 'Dataset' || e['@id'] == './') return
            String id = e['@id'].replaceAll('/$', '')
            if (!isLocal(id)) return
            def dir = new File(outdir, id)
            if (!dir.isDirectory()) return

            def parts = []
            dir.eachFileRecurse { f ->
                if (!f.isFile()) return
                String rel = new File(outdir).toPath().relativize(f.toPath())
                    .toString().replace('\\', '/')
                def child = entity(graph, byId, rel, 'File')
                if (!child.containsKey('name')) child['name'] = f.name
                integrity(child, f)
                parts << ['@id': rel]
                files++
            }
            if (parts) {
                e['hasPart'] = parts.sort { it['@id'] }
                dirs++
            }
        }
        return [dirs, files]
    }

    // ---- credit ------------------------------------------------------------

    /**
     * Contributors as Person entities behind Role entities carrying CRediT
     * URIs. The Role indirection is the schema.org pattern: the property points
     * at a Role, which repeats the property to reach the Person. It attaches a
     * role without inventing a vocabulary, and the CRediT URI is what makes it
     * machine-readable rather than a label.
     */
    private static int addPeople(graph, Map byId, Map root, Map attribution) {
        if (root == null) return 0
        int n = 0
        def contributors = attribution['contributors'] ?: []
        contributors.eachWithIndex { c, i ->
            String name = c['name']
            if (!name) return
            // Someone recorded only because they happened to be the login gets
            // an entity so the run is attributable, but no roles.
            String pid = c['orcid'] ?: "#contributor-${i + 1}"
            def person = entity(graph, byId, pid, 'Person')
            person['name'] = name
            n++

            String aff = c['affiliation']
            if (aff?.startsWith('https://ror.org/')) {
                def org = entity(graph, byId, aff, 'Organization')
                if (!org.containsKey('name')) org['name'] = aff
                link(person, 'affiliation', aff)
            }

            if (c['facility']) {
                String fid = (aff ?: '#facility') + '#' +
                             c['facility'].toLowerCase().replace(' ', '-')
                def fac = entity(graph, byId, fid, 'Organization')
                fac['name'] = c['facility']
                if (c['facility_id']) fac['identifier'] = c['facility_id']
                if (aff?.startsWith('https://ror.org/'))
                    link(fac, 'parentOrganization', aff)
                link(person, 'memberOf', fid)
            }

            def tokens = c['role_tokens'] ?: []
            if (!tokens) {
                link(root, 'author', pid)
                return
            }
            // nf-prov's own agent is already a direct author. Leaving it there
            // once the person has Roles makes root.author a mix of Person and
            // Role, which a consumer expecting authors has to pick apart. The
            // Role reaches the Person anyway.
            unlink(root, 'author', pid)
            tokens.each { tok ->
                String rid = "#role-${i + 1}-${tok}"
                def role = entity(graph, byId, rid, 'Role')
                role['roleName'] = ['@id': CREDIT + tok + '/']
                role['name'] = tok.replace('-', ' ')
                link(role, 'author', pid)
                link(root, 'author', rid)
            }
        }
        return n
    }

    private static boolean addProject(graph, Map byId, Map root, Map project) {
        if (root == null || !project) return false
        String pid = project['raid'] ?: '#project'
        def proj = entity(graph, byId, pid, 'ResearchProject')
        if (project['name']) proj['name'] = project['name']
        if (project['raid']) proj['identifier'] = project['raid']
        // The project is the organisation on whose behalf the work was done,
        // which is what sourceOrganization means. ResearchProject is an
        // Organization subtype, so the reference is type-correct.
        link(root, 'sourceOrganization', pid)

        def funder = project['funder'] ?: [:]
        if (funder['ror'] || funder['name']) {
            String fid = funder['ror'] ?: '#funder'
            def f = entity(graph, byId, fid, 'Organization')
            if (funder['name']) f['name'] = funder['name']
            // RO-Crate's guidance puts the funder on the root as well as on the
            // project, so a reader need not traverse to find it.
            link(root, 'funder', fid)
            link(proj, 'funder', fid)
            if (project['grant']) {
                String gid = '#grant-' + project['grant'].replace('/', '-')
                def g = entity(graph, byId, gid, 'Grant')
                g['identifier'] = project['grant']
                link(g, 'funder', fid)
                link(proj, 'funding', gid)
            }
        }

        def holder = project['rights_holder'] ?: [:]
        String hid = holder['orcid'] ?: holder['ror']
        if (hid) {
            def h = entity(graph, byId, hid, holder['orcid'] ? 'Person' : 'Organization')
            if (holder['name']) h['name'] = holder['name']
            link(root, 'copyrightHolder', hid)
        }
        else if (holder['name']) {
            def h = entity(graph, byId, '#rights-holder', 'Organization')
            h['name'] = holder['name']
            link(root, 'copyrightHolder', '#rights-holder')
        }
        return true
    }
}
