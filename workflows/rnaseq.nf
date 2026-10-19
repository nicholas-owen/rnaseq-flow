/*
 * Import modules
 */
include { FASTQC        } from '../modules/local/fastqc'
include { FASTP         } from '../modules/local/fastp'
include { MULTIQC       } from '../modules/local/multiqc'
include { QUARTO_REPORT } from '../modules/local/quarto'
include { STAR_ALIGN    } from '../modules/local/star'
include { HISAT2_ALIGN  } from '../modules/local/hisat2'
include { SALMON_QUANT  } from '../modules/local/salmon'
include { KALLISTO_QUANT} from '../modules/local/kallisto'
include { FEATURECOUNTS } from '../modules/local/featurecounts'
include { DESEQ2        } from '../modules/local/deseq2'
include { EDGER         } from '../modules/local/edger'
include { GSEA          } from '../modules/local/gsea'
include { GPROFILER     } from '../modules/local/gprofiler'
include { GTF2BED       } from '../modules/local/gtf2bed'
include { RSEQC         } from '../modules/local/rseqc'
include { SAMTOOLS_INDEX} from '../modules/local/samtools_index'
include { BAM_TO_BIGWIG } from '../modules/local/bam_to_bigwig'
include { RMATS         } from '../modules/local/rmats'
include { STAR_FUSION   } from '../modules/local/star_fusion'
include { ISOFORM_SWITCH} from '../modules/local/isoform_switch'
include { TXIMPORT      } from '../modules/local/tximport'
include { DTU           } from '../modules/local/dtu'
include { FEATURECOUNTS_EXON } from '../modules/local/featurecounts_exon'
include { DIFFSPLICE    } from '../modules/local/diffsplice'
include { GTF2GENEINFO  } from '../modules/local/gtf2geneinfo'

/*
 * Resolve the DESeq2/edgeR model formula:
 *   - an explicit --design wins (a leading '~' is added if missing);
 *   - else '~ batch + condition' when the samplesheet has a 'batch' column;
 *   - else the plain '~ condition'.
 */
def resolveDesign(samplesheet) {
    if (params.design) {
        def d = params.design.toString().trim()
        return d.startsWith('~') ? d : "~ ${d}"
    }
    def cols = file(samplesheet).splitCsv(header: true).first().keySet()
    return cols.contains('batch') ? '~ batch + condition' : '~ condition'
}

/*
 * An ORCID in its canonical URL form, from either the bare identifier or a URL.
 *
 * The schema accepts both because people copy whichever their profile shows.
 * Recording one shape means a reader does not have to normalise before
 * comparing, and the URL form is what CITATION.cff, .zenodo.json and the
 * Workflow Run RO-Crate agent field all use.
 */
def normaliseOrcid(orcid) {
    if (!orcid) return null
    def id = orcid.toString().trim()
    return id.startsWith('http') ? id : "https://orcid.org/${id}"
}

/*
 * A ROR organisation identifier in URL form, from either the bare id or a URL.
 * A ROR id is 0 followed by eight lowercase alphanumerics, with the last two a
 * checksum, so anything else is not one and is returned unchanged for the
 * caller to reject.
 */
def normaliseRor(ror) {
    if (!ror) return null
    def id = ror.toString().trim().replaceFirst(/^(https?:\/\/)?ror\.org\//, '')
    return id ==~ /0[a-z0-9]{8}/ ? "https://ror.org/${id}" : null
}

/*
 * A RAiD in URL form. RAiD names (ISO 23527:2022) are DataCite DOIs with a
 * raid.org URL form, so 10.12345/abc123 and https://raid.org/10.12345/abc123
 * name the same project, and either also resolves through doi.org or
 * handle.net. The URL form is stored because it is the one a reader can follow.
 *
 * The suffix is specified as alphanumeric ASCII; a dot, hyphen and underscore
 * are tolerated here so that a validly minted identifier is never rejected,
 * while a value that is not an identifier at all still is.
 */
def normaliseRaid(raid) {
    if (!raid) return null
    def id = raid.toString().trim()
                 .replaceFirst(/^(https?:\/\/)?(raid\.org|doi\.org|hdl\.handle\.net)\//, '')
    return id ==~ /10\.\d{4,9}\/[A-Za-z0-9._-]+/ ? "https://raid.org/${id}" : null
}

/*
 * A map entry as a trimmed string, or null when absent or blank.
 *
 * A plain function rather than a closure assigned to a variable: Nextflow's
 * strict parser rejects calling such a closure as `name(args)`.
 */
def trimmedField(m, k) {
    def v = m[k]?.toString()?.trim()
    return v ?: null
}

/*
 * The project this run belongs to, from params.project.
 *
 * Project-level metadata is deliberately separate from the per-sample and
 * per-run kinds. It changes rarely, it describes the research activity rather
 * than this execution of the pipeline, and it is the level at which a RAiD and
 * a funder make sense.
 *
 * Every identifier is validated and normalised. An identifier that does not
 * resolve is worse than an absent one: it reads as provenance while pointing
 * nowhere, and a reader cannot tell it from a good one without trying it.
 */
def resolveProject() {
    if (!params.containsKey('project') || !params.project) return null
    def p = params.project
    def errors = []

    def raid = trimmedField(p, 'raid')
    def raid_url = normaliseRaid(raid)
    if (raid && !raid_url) {
        errors << "project.raid: '${raid}' is not a RAiD. Expected a DOI-form identifier " +
                  "such as 10.12345/abc123, with or without the https://raid.org/ prefix"
    }

    // An identifier resolves, a name reads. Recording both means the report can
    // say "Medical Research Council" and still link the ROR behind it, rather
    // than showing a reader nine characters of hex and expecting them to click.
    def funder = trimmedField(p, 'funder')
    def funder_url = normaliseRor(funder)
    if (funder && !funder_url) {
        errors << "project.funder: '${funder}' is not a ROR identifier. Expected something " +
                  "like 03x94j517 or https://ror.org/03x94j517, from ror.org"
    }
    def funder_name = trimmedField(p, 'funder_name')
    def funder_rec = (funder_url || funder_name)
        ? [ name: funder_name, ror: funder_url ]
        : null

    // The data's owner, which is a rights relation and so deliberately not a
    // CRediT role: CRediT describes contributions to a work, not control over
    // it. A person (ORCID) or an organisation (ROR) can hold those rights, so
    // both are accepted, with rights_holder_name alongside for display.
    //
    // Where rights_holder is neither an ORCID nor a ROR it is taken as a name,
    // so naming an owner who has no identifier still works with one field.
    def holder = trimmedField(p, 'rights_holder')
    def holder_name = trimmedField(p, 'rights_holder_name')
    def holder_rec = null
    if (holder || holder_name) {
        def as_orcid = holder && holder ==~ /(https:\/\/orcid\.org\/)?\d{4}-\d{4}-\d{4}-\d{3}[\dX]/
            ? normaliseOrcid(holder) : null
        def as_ror = (holder && !as_orcid) ? normaliseRor(holder) : null
        holder_rec = [ name : holder_name ?: ((as_orcid || as_ror) ? null : holder),
                       orcid: as_orcid,
                       ror  : as_ror ]
    }

    if (errors) {
        error "Project metadata problems:\n  " + errors.join("\n  ")
    }

    return [ name         : trimmedField(p, 'name'),
             raid         : raid_url,
             rights_holder: holder_rec,
             funder       : funder_rec,
             grant        : trimmedField(p, 'grant') ]
}

/*
 * The fourteen CRediT contributor roles (ANSI/NISO Z39.104-2022), as the tokens
 * that form their URIs.
 *
 * A role is declared as the short token, because that is what anyone will type,
 * and stored in the manifest as the full URI so a reader can resolve it. The
 * two writing roles' tokens are not derivable from their display names
 * ("Writing - original draft", "Writing - review & editing"), so both were
 * checked against credit.niso.org rather than guessed.
 */
def creditRoles() {
    [ 'conceptualization', 'data-curation', 'formal-analysis',
      'funding-acquisition', 'investigation', 'methodology',
      'project-administration', 'resources', 'software', 'supervision',
      'validation', 'visualization', 'writing-original-draft',
      'writing-review-editing' ]
}

/*
 * Who worked on this run, with their CRediT roles.
 *
 * Three ways in, highest priority first: `run_contributors`, a list of maps for
 * the case of several people with different roles; `--run_by` and
 * `--run_by_orcid`, shorthand for the one-person case because a list is awkward
 * to pass on a command line; and failing both, the launching account's user
 * name.
 *
 * None of those params is declared in the project's nextflow.config, and that
 * is deliberate. Nextflow merges config files with $HOME/.nextflow/config at
 * the LOWEST priority, below the project's own nextflow.config, so a default
 * declared there would silently override whatever someone had set once in their
 * home config, which is the natural place to record who you are. They are
 * declared in nextflow_schema.json instead, which is what the typo check reads,
 * so `--run_by` is still a known parameter without a project-level default
 * outranking the user's.
 *
 * An unknown role token is an error, not a warning. A misspelled token produces
 * metadata that claims CRediT conformance for a role that does not exist, and
 * nothing downstream would catch it. Every problem is reported at once, as the
 * samplesheet validator does.
 */
def resolveContributors() {
    def valid  = creditRoles()
    def errors = []

    def entries = null
    if (params.containsKey('run_contributors') && params.run_contributors) {
        entries = params.run_contributors instanceof List
            ? params.run_contributors
            : [params.run_contributors]
    }
    if (!entries) {
        def named = params.containsKey('run_by') ? params.run_by : null
        entries = named
            ? [[ 'name' : named,
                 'orcid': params.containsKey('run_by_orcid') ? params.run_by_orcid : null ]]
            : []
    }

    // Whether anybody identified themselves at all. The fallback records the
    // login that happened to be in use, which is not a claim anyone made.
    def declared = !entries.isEmpty()
    if (!declared) entries = [[ 'name': workflow.userName ]]

    def contributors = []
    entries.eachWithIndex { e, i ->
        def name = e['name']?.toString()?.trim()
        if (!name) errors << "contributor ${i + 1}: 'name' is required"

        def roles = e['roles'] ?: []
        if (!(roles instanceof List)) roles = [roles]
        roles = roles.collect { it.toString().trim().toLowerCase() }.findAll { it }
        roles.each { r ->
            if (!valid.contains(r)) {
                errors << "contributor ${i + 1} (${name ?: 'unnamed'}): '${r}' is not a CRediT role"
            }
        }

        // Running the pipeline is formal analysis under CRediT, which covers
        // applying computational techniques to analyse data. So default to it
        // for someone who identified themselves, never for the fallback, and
        // record that it was inferred rather than claimed.
        def roles_declared = !roles.isEmpty()
        if (!roles_declared && declared) roles = ['formal-analysis']

        contributors << [
            name          : name,
            orcid         : normaliseOrcid(e['orcid']),
            roles         : roles.collect { "https://credit.niso.org/contributor-roles/${it}/" },
            role_tokens   : roles,
            roles_declared: roles_declared,
            // A ROR is normalised to its URL so it resolves; anything else is
            // kept as typed, because an institution without a ROR is still
            // worth naming.
            affiliation   : normaliseRor(e['affiliation']) ?: e['affiliation']?.toString(),
            facility      : e['facility']?.toString(),
            facility_id   : e['facility_id']?.toString(),
            source        : declared ? 'declared' : 'system_user',
        ]
    }

    if (errors) {
        error "Run contributor problems:\n  " + errors.join("\n  ") +
              "\n\nValid CRediT roles:\n  " + valid.join(', ')
    }
    return contributors
}

/*
 * Resolve where the annotation came from, for the analysis report.
 *
 * `--download_refs` writes reference_metadata.json beside the files it fetches
 * (assets/download_refs.py). An analysis run never downloads anything -- main.nf
 * makes --download_refs, --build_indices and the analysis run mutually exclusive
 * branches -- so that record has to be *found on disk* next to the GTF rather
 * than handed over in a channel.
 *
 * Always returns a map with the same keys and an `origin` of:
 *   'pipeline_download' - a sidecar was found and it describes this annotation
 *   'user_supplied'     - no sidecar, or one that describes something else
 *
 * Fields that could not be established are null, and `note` says why. Callers
 * should render those as "not recorded" rather than omitting them or guessing:
 * an Ensembl GTF's header carries the assembly, the accession and the genebuild
 * date but never the release, and the trailing number in the filename is an
 * annotation version rather than the release (release-116 ships
 * Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz). A release inferred from a
 * filename is a false record that outlives the run and cannot afterwards be
 * told apart from a true one.
 */
def resolveReferenceProvenance(gtf) {
    def base = [
        origin            : 'user_supplied',
        file              : null,
        dir               : null,
        bytes             : null,
        species           : null,
        source            : null,
        release           : null,
        release_requested : null,
        assembly          : null,
        annotation_version: null,
        downloaded_at     : null,
        url               : null,
        sha256            : null,
        note              : null,
    ]

    if (!gtf) {
        return base + [note: 'no --gtf supplied']
    }

    // file() returns a *list* when the argument globs, and the documented usage
    // does exactly that (--gtf refs/yeast/v116/*.gtf.gz). Anything other than a
    // single match leaves us unable to say which file this is about.
    def resolved = file(gtf)
    if (resolved instanceof List) {
        if (resolved.size() != 1) {
            return base + [note: "--gtf matched ${resolved.size()} files; cannot identify one annotation"]
        }
        resolved = resolved[0]
    }
    if (!resolved.exists()) {
        return base + [note: 'annotation not found on disk']
    }

    def found = base + [file : resolved.name,
                        dir  : resolved.parent?.toString(),
                        bytes: resolved.size()]

    def sidecar = file("${resolved.parent}/reference_metadata.json")
    if (!sidecar.exists()) {
        return found + [note: 'no reference_metadata.json beside the annotation']
    }

    def meta = null
    try {
        meta = new groovy.json.JsonSlurper().parseText(sidecar.text)
    } catch (Exception e) {
        return found + [note: "reference_metadata.json could not be read: ${e.message}"]
    }

    // The sidecar must be shown to describe *this* file. A reference directory
    // someone has added a hand-fetched GTF to would otherwise have confident,
    // wrong provenance attached to it -- worse than reporting nothing.
    def ann = meta instanceof Map ? meta.files?.annotation : null
    if (!(ann instanceof Map) || !ann.file) {
        return found + [note: 'reference_metadata.json has no annotation record']
    }
    if (ann.file != resolved.name) {
        return found + [note: "reference_metadata.json describes '${ann.file}', not this annotation"]
    }
    // Size, not checksum: hashing a multi-gigabyte GTF on the head node would
    // add real time to the launch of every run. The name-plus-size pair is a
    // cheap discriminator, and the recorded sha256 is carried through for
    // anyone who wants to verify deliberately.
    if (ann.bytes != null && resolved.size() != ann.bytes) {
        return found + [note: "annotation is ${resolved.size()} bytes but " +
                              "reference_metadata.json records ${ann.bytes}; treating as user-supplied"]
    }

    return found + [
        origin            : 'pipeline_download',
        species           : meta.species,
        source            : meta.source,
        release           : meta.release,
        release_requested : meta.release_requested,
        assembly          : meta.assembly,
        annotation_version: meta.annotation_version,
        downloaded_at     : meta.downloaded_at,
        url               : ann.url,
        sha256            : ann.sha256,
    ]
}

/*
 * Describe a GMT gene-set file, as far as it can be described.
 *
 * GMT is a bare tab-separated format with no header, so a gene-set file carries
 * no version, no build date and no provenance of any kind -- the first line is
 * already a pathway. The collection name can only be taken from the filename,
 * which is inference, and is flagged as such. `--download_gmt` does not yet
 * record what it fetched (msigdbr knows its own version at download time but
 * writes no metadata), so `version` is always null for now.
 */
def describeGeneSets(gmt) {
    def base = [path: null, file: null, bytes: null, collection: null,
                collection_inferred: false, version: null, note: null]

    if (!gmt) {
        return base + [note: 'no --gmt supplied; GSEA did not run']
    }

    def resolved = file(gmt)
    if (resolved instanceof List) {
        if (resolved.size() != 1) {
            return base + [path: gmt.toString(),
                           note: "--gmt matched ${resolved.size()} files"]
        }
        resolved = resolved[0]
    }
    if (!resolved.exists()) {
        return base + [path: gmt.toString(), note: 'gene-set file not found on disk']
    }

    return base + [
        path               : resolved.toString(),
        file               : resolved.name,
        bytes              : resolved.size(),
        collection         : resolved.name.replaceFirst(/\.gmt$/, ''),
        collection_inferred: true,
        note               : 'GMT files carry no version metadata; the collection is read from the filename',
    ]
}

/*
 * Assemble the run manifest: the record of *how this run was configured*,
 * for the analysis report.
 *
 * The report is the artefact that travels -- it gets emailed to people with no
 * access to the run directory or the command line. Without this it cannot say
 * what was aligned, against what, or with which gene sets, because none of that
 * reaches the reporting step otherwise.
 *
 * This is deliberately a record of *configuration*, not of outcome: it is built
 * before any process runs, so it describes what was asked for. What actually
 * got produced is something the report already establishes for itself, from
 * which result directories exist. Keeping the two apart avoids the manifest
 * claiming a stage ran when it failed or was skipped.
 */
/*
 * Whether the run is paired- or single-end, read from the samplesheet.
 *
 * Taken from the sheet rather than from the read channel because the manifest
 * is assembled synchronously, before any channel has emitted. A cohort may
 * legitimately mix the two, and that is reported as 'mixed' rather than
 * flattened to whichever layout the first row happens to have -- a mixed cohort
 * is worth a reader noticing.
 *
 * Returns null when there is no sheet to read: a glob input carries no R2
 * column, and a --counts run has no reads at all.
 */
def resolveLibraryLayout(samplesheet) {
    if (!samplesheet) {
        return null
    }
    def sheet = samplesheet.toString()
    if (!sheet.endsWith('.csv')) {
        return null
    }
    def rows = file(sheet).splitCsv(header: true)
    if (!rows) {
        return null
    }
    if (!rows.first().containsKey('R2')) {
        return 'single-end'
    }
    def paired = rows.count { row -> row.R2?.toString()?.trim() }
    if (paired == 0)          return 'single-end'
    if (paired == rows.size()) return 'paired-end'
    return "mixed (${paired} of ${rows.size()} paired)"
}

/*
 * Summarise the experimental design from the samplesheet: how many samples, in
 * which conditions, with what replication.
 *
 * Conditions are kept in first-appearance order rather than sorted. The `REF`
 * baseline convention means the control is written first, and every fold change
 * in the report is oriented against it, so alphabetical order would present the
 * design in a sequence that contradicts how the results were computed.
 *
 * Returns null when there is nothing to read -- a glob input, or a sheet
 * missing the columns this describes.
 */
def summariseSamples(samplesheet) {
    if (!samplesheet) {
        return null
    }
    def sheet = samplesheet.toString()
    if (!sheet.endsWith('.csv')) {
        return null
    }
    def rows = file(sheet).splitCsv(header: true)
    if (!rows) {
        return null
    }
    def cols = rows.first().keySet()
    if (!cols.contains('sample') || !cols.contains('condition')) {
        return null
    }

    def order  = []
    def byCond = [:]
    rows.each { row ->
        def cond = row.condition?.toString()?.trim()
        def name = row.sample?.toString()?.trim()
        if (cond) {
            if (!byCond.containsKey(cond)) {
                byCond[cond] = []
                order.add(cond)
            }
            byCond[cond].add(name)
        }
    }

    def batches = []
    if (cols.contains('batch')) {
        batches = rows.collect { row -> row.batch?.toString()?.trim() }
                      .findAll { value -> value }
                      .unique()
    }

    return [
        total     : rows.size(),
        conditions: order.collect { cond ->
            [name: cond, n: byCond[cond].size(), samples: byCond[cond]]
        },
        has_batch : cols.contains('batch'),
        batches   : batches,
    ]
}

def buildRunManifest(design, notices) {
    def from_counts = params.counts != null

    def quantification = null
    if (from_counts) {
        quantification = 'supplied count matrix'
    } else if (params.aligner == 'star' || params.aligner == 'hisat2') {
        quantification = 'featureCounts, gene level from BAM'
    } else if (params.aligner == 'salmon' || params.aligner == 'kallisto') {
        quantification = 'tximport, transcript quantification summarised to gene level'
    }

    def index = null
    if (params.aligner == 'star')          index = params.star_index
    else if (params.aligner == 'hisat2')   index = params.hisat2_index
    else if (params.aligner == 'salmon')   index = params.salmon_index
    else if (params.aligner == 'kallisto') index = params.kallisto_index

    // Only the optional analyses actually switched on, so the report can list
    // them without a row of "false" for everything nobody asked for.
    def optional = []
    if (params.isoform_switch) optional.add('isoform switching (IsoformSwitchAnalyzeR)')
    if (params.dtu)            optional.add('differential transcript usage (DEXSeq)')
    if (params.diffsplice)     optional.add('differential splicing (edgeR diffSplice)')
    if (params.ctat_lib)       optional.add('fusion calling (STAR-Fusion)')

    return [
        // 2 adds the `attribution` block. Additive, so a reader written against
        // 1 still works.
        schema_version: 2,
        pipeline: [
            name            : workflow.manifest.name,
            version         : workflow.manifest.version,
            home_page       : workflow.manifest.homePage,
            // null for a local run; set when launched with `nextflow run <repo> -r <rev>`.
            revision        : workflow.revision,
            commit_id       : workflow.commitId,
            run_name        : workflow.runName,
            session_id      : workflow.sessionId?.toString(),
            nextflow_version: nextflow.version?.toString(),
            started_at      : workflow.start?.toString(),
            profile         : workflow.profile,
            container_engine: workflow.containerEngine,
            // The verbatim invocation. Structured fields above are what the
            // report renders; this is the fallback that answers questions the
            // schema did not anticipate.
            command_line    : workflow.commandLine,
            // Nextflow exposes no per-parameter provenance, so there is no way
            // to record which file a contributor's details came from. The list
            // of config files that were in play is the next best thing: a
            // reader chasing an attribution at least knows where to look.
            config_files    : workflow.configFiles?.collect { it.toString() },
        ],
        // Who to credit, and under what terms. The results travel: a reader who
        // has only a results directory can otherwise tell what was run but not
        // whose it is or what they may do with it.
        //
        // Two different claims, kept apart. `authors` and `license` describe the
        // pipeline and come from the manifest scope, so they are the same for
        // every run. `contributors` describes the people who did this analysis,
        // each with their CRediT roles, and each carrying `source` so a name
        // somebody declared is distinguishable from the login that happened to
        // be in use. See resolveContributors().
        attribution: [
            license     : workflow.manifest.license,
            doi         : workflow.manifest.doi,
            authors     : (workflow.manifest.contributors ?: []).collect { c ->
                [ name        : c?.name?.toString(),
                  affiliation : c?.affiliation?.toString(),
                  orcid       : c?.orcid?.toString(),
                  // Nextflow returns these as enum names (AUTHOR, MAINTAINER);
                  // lower-cased to match how they are written in the manifest
                  // scope and everywhere else they appear.
                  contribution: c?.contribution?.collect { it.toString().toLowerCase() } ]
            },
            contributors: resolveContributors(),
        ],
        // The research activity this run belongs to, where one has been
        // declared. Null when it has not, and the report skips the group
        // rather than showing a row of blanks.
        project: resolveProject(),
        analysis: [
            entry         : from_counts ? 'count matrix' : 'sequencing reads',
            // The samplesheet or read glob, and the matrix for a --counts run.
            // Named so a reader can tell which cohort produced the report.
            input         : params.input?.toString(),
            counts        : params.counts?.toString(),
            library_layout: from_counts ? null : resolveLibraryLayout(params.input),
            aligner       : from_counts ? null : params.aligner,
            quantification: quantification,
            strandedness  : from_counts ? null : params.strandedness,
            design        : design,
            // The denominator of every contrast, so it decides which way every
            // fold change and PSI difference points -- worth recording in the
            // artefact that travels.
            reference_level: params.reference_level?.toString(),
            organism      : params.organism,
            stop_at       : params.stop_at,
            optional      : optional,
            outdir        : params.outdir?.toString(),
        ],
        // The experimental design, read from the samplesheet. The report has no
        // other route to it -- the samplesheet is never staged for rendering.
        samples: summariseSamples(params.input),
        // Analyses that were requested and will not run. Empty on a clean run,
        // and the report renders nothing for an empty list.
        notices: notices ?: [],
        references: [
            annotation      : resolveReferenceProvenance(params.gtf),
            genome_fasta    : params.genome_fasta?.toString(),
            transcript_fasta: params.transcript_fasta?.toString(),
            index           : index?.toString(),
        ],
        gene_sets: describeGeneSets(params.gmt),
    ]
}

/*
 * Main workflow
 */
workflow RNASEQ {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()

    // Notices about analyses that were asked for and will not happen.
    //
    // These are the dangerous class of message: the run succeeds, every process
    // ticks green, and part of what was requested is quietly absent. A failure
    // announces itself; this does not. Each is logged as a warning *and*
    // recorded here, so it reaches the report rather than living only in
    // terminal scrollback that has scrolled away by the time anyone reads the
    // results.
    //
    // Accumulated during workflow construction, which is synchronous, so the
    // list is complete by the time buildRunManifest() reads it.
    def notices = []

    // Result directories fed to the Quarto analysis report. Each defaults to an
    // empty list and is reassigned to the real channel if that stage runs, so
    // the report renders whichever sections have data.
    ch_quarto_deseq2     = channel.value([])
    ch_quarto_edger      = channel.value([])
    ch_quarto_gsea       = channel.value([])
    ch_quarto_gprofiler  = channel.value([])
    ch_quarto_rmats      = channel.value([])
    ch_quarto_dtu        = channel.value([])
    ch_quarto_diffsplice = channel.value([])
    ch_quarto_isoform    = channel.value([])

    // Count-matrix entry: when --counts is set the pipeline skips QC and
    // alignment entirely and enters at differential expression (Level 3),
    // consuming the pre-computed matrix instead of FASTQs.
    def from_counts = params.counts != null

    // Define execution levels
    def run_level = 100 // Default: run all
    if (params.stop_at == 'preQC')       run_level = 1
    else if (params.stop_at == 'postQC') run_level = 2
    else if (params.stop_at == 'DE')     run_level = 3
    else if (params.stop_at == 'GSEA')   run_level = 4
    
    //
    // MODULE: FastQC (Raw) (Level >= 1)
    //
    if (run_level >= 1 && !from_counts) {
        FASTQC ( ch_reads )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        
        //
        // MODULE: fastp (Trimming)
        //
        FASTP ( ch_reads )
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }

    //
    // MODULE: Alignment / Quantification (Level >= 2)
    //
    ch_align_results = channel.empty()
    ch_bams          = channel.empty()

    // Initialize downstream channels to empty in case skipped.
    // ch_multiqc_files is declared here (not later) so it can be safely
    // referenced from inside conditional blocks below.
    ch_featurecounts_results = channel.empty()  // count files -> DESeq2 / edgeR
    ch_featurecounts_mqc     = channel.empty()  // .summary files -> MultiQC
    ch_featurecounts_exon    = channel.empty()  // exon counts -> diffSplice
    ch_rseqc_results         = channel.empty()
    ch_starfusion_results    = channel.empty()
    ch_multiqc_files         = channel.empty()
    ch_tx_quant              = channel.empty()  // Salmon/Kallisto quant dirs -> tximport

    if (run_level >= 2 && !from_counts) {
        if (params.aligner == 'star') {
            if (!params.star_index) error "STAR index not provided via --star_index"
            STAR_ALIGN( FASTP.out.reads, file(params.star_index), file(params.gtf) ) 
            ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
            ch_align_results = ch_align_results.mix(STAR_ALIGN.out.log_final.collect { row -> row[1] })
            ch_bams = STAR_ALIGN.out.bam
            
            // STAR-Fusion (runs if Fusion is requested AND we are at Level 4+ OR just run it if CTAT provided? 
            // Usually fusion is advanced analysis. Let's put it in Level 4 as per plan).
            if (run_level >= 4 && params.ctat_lib) {
                STAR_FUSION (
                    STAR_ALIGN.out.chimeric_junction,
                    file(params.ctat_lib)
                )
                ch_versions = ch_versions.mix(STAR_FUSION.out.versions.first())
                ch_starfusion_results = ch_starfusion_results.mix(STAR_FUSION.out.fusions_abridged)
            }
        }
        else if (params.aligner == 'hisat2') {
            if (!params.hisat2_index) error "HISAT2 index not provided via --hisat2_index"
            HISAT2_ALIGN( FASTP.out.reads, file(params.hisat2_index) )
            ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
            ch_align_results = ch_align_results.mix(HISAT2_ALIGN.out.summary.collect { row -> row[1] })
            ch_bams = HISAT2_ALIGN.out.bam
        }
        else if (params.aligner == 'salmon') {
            if (!params.salmon_index) error "Salmon index not provided via --salmon_index"
            SALMON_QUANT( FASTP.out.reads, file(params.salmon_index) )
            ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
            ch_align_results = ch_align_results.mix(SALMON_QUANT.out.results.collect { row -> row[1] })
            ch_tx_quant = SALMON_QUANT.out.results.collect { row -> row[1] }
            
            // Isoform Switch (Level 4+)
            if (run_level >= 4 && params.isoform_switch && params.transcript_fasta) {
                ISOFORM_SWITCH (
                    file(params.input),
                    file(params.transcript_fasta),
                    file(params.gtf),
                    SALMON_QUANT.out.results.collect { row -> row[1] },
                    file("${projectDir}/assets/isoform_switch.R")
                )
                ch_versions = ch_versions.mix(ISOFORM_SWITCH.out.versions)
                ch_quarto_isoform = ISOFORM_SWITCH.out.results
            }
        }
        else if (params.aligner == 'kallisto') {
            if (!params.kallisto_index) error "Kallisto index not provided via --kallisto_index"
            if (params.strandedness == 'auto') {
                def msg = "kallisto cannot infer strandedness (it produces no BAM for RSeQC). " +
                          "Running library-type-agnostic; pass --strandedness forward or reverse " +
                          "explicitly if your library is stranded."
                log.warn msg
                notices.add(msg)
            }
            KALLISTO_QUANT( FASTP.out.reads, file(params.kallisto_index) )
            ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
            ch_align_results = ch_align_results.mix(KALLISTO_QUANT.out.results.collect { row -> row[1] })
            ch_tx_quant = KALLISTO_QUANT.out.results.collect { row -> row[1] }
        }
        
        //
        // Post Alignment QC / Processing (Star/Hisat)
        //
        ch_bam_bai = channel.empty()
        
        if (params.aligner == 'star' || params.aligner == 'hisat2') {
             SAMTOOLS_INDEX ( ch_bams )
             ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
             ch_bam_bai = SAMTOOLS_INDEX.out.bam_bai
             
             // RSeQC (PostQC)
             if (params.gtf) {
                 GTF2BED ( file(params.gtf), file("${projectDir}/assets/gtf2bed12.py") )
                 ch_versions = ch_versions.mix(GTF2BED.out.versions)
                 
                 RSEQC ( ch_bam_bai, GTF2BED.out.bed, file("${projectDir}/assets/parse_strandedness.py") )
                 ch_versions = ch_versions.mix(RSEQC.out.versions.first())
                 ch_rseqc_results = ch_rseqc_results.mix(RSEQC.out.infer_experiment.collect { row -> row[1] })
                 ch_rseqc_results = ch_rseqc_results.mix(RSEQC.out.read_distribution.collect { row -> row[1] })
                 ch_rseqc_results = ch_rseqc_results.mix(RSEQC.out.genebody_coverage.collect { row -> row[1] })
             }
             
             // BigWig
             BAM_TO_BIGWIG( ch_bam_bai )
             ch_versions = ch_versions.mix(BAM_TO_BIGWIG.out.versions.first())
             
             // rMATS (Level 4+) - collapse all per-sample BAMs into one task.
             if (run_level >= 4 && params.input.endsWith('.csv')) {
                ch_rmats_input = ch_bam_bai
                    // meta.id is carried through so the module can emit an
                    // authoritative sample->BAM mapping. It used to be dropped
                    // here, which forced run_rmats.py to reconstruct the pairing
                    // from filenames with a substring match -- and sample
                    // 'ctrl1' matched 'ctrl10.bam' (C3): wrong BAM, wrong
                    // splicing calls, no error.
                    .map { meta, bam, bai -> [ meta.id, bam, bai, meta.single_end ] }
                    // flat:false keeps each [id, bam, bai, single_end] as a
                    // sub-list instead of flattening everything into one list.
                    .collect(flat: false)
                    .map { rows ->
                         // rMATS runs one -t mode for the whole cohort, so a
                         // samplesheet mixing single- and paired-end cannot be
                         // analysed correctly -- previously row 0's layout was
                         // silently applied to everyone (H7).
                         def layouts = rows.collect { row -> row[3] }.unique()
                         if (layouts.size() > 1) {
                             error "rMATS needs a single library layout, but the samplesheet " +
                                   "mixes single-end and paired-end samples."
                         }
                         def ids  = rows.collect { row -> row[0] }
                         def bams = rows.collect { row -> row[1] }
                         def bais = rows.collect { row -> row[2] }
                         [ [ id:'all_samples', single_end: rows[0][3], ids: ids ], bams, bais ]
                    }
                RMATS ( file(params.input), ch_rmats_input, file(params.gtf), file("${projectDir}/assets/run_rmats.py") )
                ch_versions = ch_versions.mix(RMATS.out.versions)
                ch_quarto_rmats = RMATS.out.results
             }
         
             // FeatureCounts gene counts (input to DESeq2 / edgeR at Level >= 3).
             if (params.gtf) {
                 def gtf_file = file(params.gtf)

                 // Per-sample strandedness for featureCounts: the RSeQC
                 // inference when --strandedness is 'auto', otherwise the
                 // user-supplied --strandedness value applied to every sample.
                 ch_strand = (params.strandedness == 'auto')
                     ? RSEQC.out.strandedness.map { meta, txt -> [ meta, txt.text.trim() ] }
                     : ch_bams.map { meta, _bam -> [ meta, params.strandedness ] }

                 ch_fc_input = ch_bams.join(ch_strand)
                     .map { meta, bam, strand -> [ meta, bam, gtf_file, strand ] }

                 FEATURECOUNTS ( ch_fc_input )
                 ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions.first())
                 ch_featurecounts_results = FEATURECOUNTS.out.counts.collect { row -> row[1] }
                 ch_featurecounts_mqc     = FEATURECOUNTS.out.summary.collect { row -> row[1] }

                 // Exon-level counts for the optional edgeR diffSplice test.
                 if (run_level >= 3 && params.diffsplice && params.input.endsWith('.csv')) {
                     FEATURECOUNTS_EXON ( ch_fc_input )
                     ch_versions = ch_versions.mix(FEATURECOUNTS_EXON.out.versions.first())
                     ch_featurecounts_exon = FEATURECOUNTS_EXON.out.counts.collect { row -> row[1] }
                 }
             }
        }
    }
    
    //
    // MODULE: Differential Expression (Level >= 3)
    //
    // STAR/HISAT2 feed DESeq2 + edgeR via the featureCounts gene matrix.
    // Salmon/Kallisto feed them via TXIMPORT, which summarises transcript-level
    // quantification to gene level (this is what enables transcript-based
    // aligners to run gene-level differential expression).
    //
    if (run_level >= 3 && params.input.endsWith('.csv')) {
        ch_de_counts = channel.empty()
        ch_gene_info = channel.empty()
        def run_de   = false

        // Gene-symbol / biotype table parsed once from the GTF; used to
        // annotate the DESeq2 / edgeR result tables with readable symbols.
        if (params.gtf) {
            GTF2GENEINFO ( file(params.gtf), file("${projectDir}/assets/gtf2geneinfo.py") )
            ch_versions  = ch_versions.mix(GTF2GENEINFO.out.versions)
            ch_gene_info = GTF2GENEINFO.out.gene_info
        }

        if (from_counts) {
            // Pre-computed matrix: feed it straight to DESeq2/edgeR (both R
            // scripts read it in --matrix mode). No aligner ran.
            ch_de_counts = channel.fromPath(params.counts)
            run_de = true
        }
        else if (params.aligner == 'star' || params.aligner == 'hisat2') {
            ch_de_counts = ch_featurecounts_results
            run_de = true
            if (params.dtu) {
                def msg = "--dtu (differential transcript usage) needs a pseudo-aligner; " +
                          "it has no effect with --aligner ${params.aligner}."
                log.warn msg
                notices.add(msg)
            }

            // Optional: edgeR diffSplice exon-usage test (genome aligners).
            if (params.diffsplice) {
                if (params.gtf) {
                    DIFFSPLICE (
                        file(params.input),
                        ch_featurecounts_exon,
                        file(params.gtf),
                        'exon',
                        params.aligner,
                        file("${projectDir}/assets/diffsplice_edger.R")
                    )
                    ch_versions = ch_versions.mix(DIFFSPLICE.out.versions)
                    // Assigned in both DIFFSPLICE branches. They are mutually
                    // exclusive at run time (a run has one aligner), so only
                    // one ever executes, but leaving either out would drop the
                    // section from the report on that route.
                    ch_quarto_diffsplice = DIFFSPLICE.out.results
                } else {
                    def msg = "--diffsplice needs --gtf for exon-level counting; skipping."
                    log.warn msg
                    notices.add(msg)
                }
            }
        }
        else if (params.aligner == 'salmon' || params.aligner == 'kallisto') {
            if (params.gtf) {
                TXIMPORT (
                    file(params.input),
                    ch_tx_quant,
                    file(params.gtf),
                    params.aligner,
                    file("${projectDir}/assets/tximport.R")
                )
                ch_versions  = ch_versions.mix(TXIMPORT.out.versions)
                ch_de_counts = TXIMPORT.out.txi
                run_de = true

                // Optional: differential transcript usage (DEXSeq).
                if (params.dtu) {
                    DTU (
                        file(params.input),
                        ch_tx_quant,
                        file(params.gtf),
                        params.aligner,
                        file("${projectDir}/assets/dtu_dexseq.R")
                    )
                    ch_versions = ch_versions.mix(DTU.out.versions)
                    ch_quarto_dtu = DTU.out.results
                }

                // Optional: edgeR diffSplice transcript-usage test.
                if (params.diffsplice) {
                    DIFFSPLICE (
                        file(params.input),
                        ch_tx_quant,
                        file(params.gtf),
                        'transcript',
                        params.aligner,
                        file("${projectDir}/assets/diffsplice_edger.R")
                    )
                    ch_versions = ch_versions.mix(DIFFSPLICE.out.versions)
                    ch_quarto_diffsplice = DIFFSPLICE.out.results
                }
            } else {
                def msg = "Differential expression for --aligner ${params.aligner} needs --gtf " +
                          "(transcript-to-gene map); skipping DESeq2/edgeR."
                log.warn msg
                notices.add(msg)
            }
        }

        if (run_de) {
            def de_design = resolveDesign(params.input)
            log.info "Differential-expression design: ${de_design}"

            DESEQ2 ( file(params.input), de_design, ch_gene_info, ch_de_counts, file("${projectDir}/assets/deseq2.R") )
            ch_versions = ch_versions.mix(DESEQ2.out.versions)
            ch_quarto_deseq2 = DESEQ2.out.results

            EDGER ( file(params.input), de_design, ch_gene_info, ch_de_counts, file("${projectDir}/assets/edger.R") )
            ch_versions = ch_versions.mix(EDGER.out.versions)
            ch_quarto_edger = EDGER.out.results

            //
            // MODULE: Enrichment (Level >= 4)
            //
            if (run_level >= 4) {
                if (params.gmt) {
                     GSEA ( file(params.gmt), DESEQ2.out.results, file("${projectDir}/assets/gsea.R") )
                     ch_versions = ch_versions.mix(GSEA.out.versions)
                     ch_quarto_gsea = GSEA.out.results
                }
                GPROFILER ( params.organism, DESEQ2.out.results, file("${projectDir}/assets/gprofiler.R") )
                ch_versions = ch_versions.mix(GPROFILER.out.versions)
                ch_quarto_gprofiler = GPROFILER.out.results
            }
        }
    }
    else if (run_level >= 3 && !params.input.endsWith('.csv')) {
        def msg = "Differential expression was skipped: the input is a file glob, not a CSV " +
                  "samplesheet, so there are no conditions to contrast."
        log.warn msg
        notices.add(msg)
    }

    //
    // MODULE: MultiQC
    //
    // Skipped entirely for a --counts run: there are no QC / alignment logs to
    // aggregate, so `multiqc .` would find no modules and fail. The Quarto
    // report then receives an empty MultiQC directory and renders without its
    // QC section (the .qmd already guards every QC chunk on the data existing).
    //
    // ch_multiqc_files was initialised to channel.empty() at the top of the
    // workflow, so the conditional mixes below are always safe.
    ch_multiqc_data = channel.value([])
    if (!from_counts) {
        if (run_level >= 1) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { row -> row[1] })
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.html.collect { row -> row[1] })
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { row -> row[1] })
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { row -> row[1] })
        }

        ch_multiqc_files = ch_multiqc_files.mix(ch_align_results)
        ch_multiqc_files = ch_multiqc_files.mix(ch_featurecounts_mqc)
        ch_multiqc_files = ch_multiqc_files.mix(ch_rseqc_results)
        ch_multiqc_files = ch_multiqc_files.mix(ch_starfusion_results)

        // .collect() so MultiQC runs ONCE over all staged files, not once per file.
        // The custom config (auto-detected by `multiqc .`) styles the report; the
        // CSS + logo it references are staged alongside it so MultiQC finds them.
        MULTIQC (
            file("${projectDir}/assets/multiqc_config.yml"),
            [ file("${projectDir}/assets/multiqc_custom.css"),
              file("${projectDir}/assets/multiqc_logo.png") ],
            ch_versions.unique().collectFile(name: 'software_versions.yml'),
            [],
            ch_multiqc_files.collect(),
            // The samplesheet is passed so MULTIQC can map each tool's own
            // sample naming onto the sample ids; see the note in the module.
            // Only a CSV samplesheet carries that mapping -- a glob input has
            // no sample column, so nothing is staged and MultiQC is left alone.
            params.input && params.input.endsWith('.csv') ? file(params.input) : []
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions.first())
        ch_multiqc_data = MULTIQC.out.data
    }

    //
    // Run manifest: how this run was configured, for the report.
    //
    // The design is resolved here only when the input is a samplesheet -- a
    // glob input has no condition column to read one from.
    //
    // collectFile writes the JSON and hands back its path, so the record stays
    // inside Nextflow's staging rather than being written to disk on the side;
    // storeDir also publishes it to pipeline_info/, where it is independently
    // useful to anything scripting across runs. .first() makes it a value
    // channel, matching the other inputs to QUARTO_REPORT.
    def manifest_design = params.input && params.input.toString().endsWith('.csv')
        ? resolveDesign(params.input)
        : null
    ch_run_manifest = channel
        .of(groovy.json.JsonOutput.prettyPrint(
                groovy.json.JsonOutput.toJson(buildRunManifest(manifest_design, notices))))
        .collectFile(name: 'run_manifest.json', storeDir: "${params.tracedir}")
        .first()

    //
    // MODULE: Quarto analysis report
    //
    // Always runs; the DE/enrichment result directories are passed when those
    // stages ran, else an empty list, and the report renders whichever sections
    // have data. multiqc_data is empty for a --counts run (see above).
    QUARTO_REPORT (
        ch_multiqc_data,
        file("${projectDir}/assets/analysis_report.qmd"),
        ch_run_manifest,
        file("${projectDir}/assets/rnaseq-flow_logo.svg"),
        ch_quarto_deseq2,
        ch_quarto_edger,
        ch_quarto_gsea,
        ch_quarto_gprofiler,
        ch_quarto_rmats,
        ch_quarto_dtu,
        ch_quarto_diffsplice,
        ch_quarto_isoform
    )
    // Single unnamed emit (the strict parser wants the name omitted when there
    // is only one): the accumulated version channel, plus the reporting steps'.
    emit:
    ch_versions.mix(QUARTO_REPORT.out.versions)
}
