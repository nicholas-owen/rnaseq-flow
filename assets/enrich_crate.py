"""Add integrity and credit metadata to the Workflow Run RO-Crate that nf-prov
writes, using the run manifest as the source of both.

    python3 assets/enrich_crate.py <outdir>

nf-prov writes the expensive half of the crate: one CreateAction per task, the
file and dataset entities, and conformance to the three Workflow Run RO-Crate
profiles. What it cannot write is everything its config has no slot for. Its
wrroc format accepts exactly one `agent`, one `organization`, a `license` and a
`publisher`, so a run with several contributors, CRediT roles, a project, a
funder or a rights holder loses all of it. It also records no checksums, sizes
or modification times at all, which leaves a crate that says which files exist
and gives a reader no way to check the bytes in hand are those files.

This fills both gaps from pipeline_info/run_manifest.json, which already holds
every identifier in normalised, resolvable form.

What it adds
------------
1. sha256, contentSize and dateModified on every File entity that resolves to a
   real file, and on every file inside a directory-valued output.
2. hasPart on the directory outputs (deseq2_output and friends), which nf-prov
   leaves as bare Dataset entities with only a name, so the result tables are
   described rather than referenced as opaque folders.
3. Every contributor as a Person behind a Role carrying its CRediT URI, their
   affiliation as a ROR Organization, and their core facility as an
   Organization carrying its RRID.
4. The project as a ResearchProject keyed by its RAiD, the funder as a ROR
   Organization, and the rights holder.

Term choices, and why
---------------------
The RO-Crate 1.1 context is schema.org plus a little, 2627 terms, and two of
the things needed are not among them.

`sha256` is not defined, and schema.org has no checksum property at all, so the
crate's @context is extended with one term pointing at SPDX, a long-standing
vocabulary that does define it. The alternative, a PropertyValue entity per
file, would have added several hundred entities to say the same thing. If the
handoff profile in the design notes is ever written, this is a term it should
pin rather than leave to each implementation.

`rightsHolder` is not defined either; `copyrightHolder` is, and is used instead.

A project is given @type ResearchProject, which is in the context and is a
subtype of Organization, so it satisfies RO-Crate's own guidance that a project
be modelled as an Organization while still saying what it actually is.

Idempotent: running it twice produces the same crate, so it is safe to re-run
after a resume.
"""
import argparse
import datetime
import hashlib
import json
import os
import sys

SPDX_CHECKSUM = "http://spdx.org/rdf/terms#checksumValue"
CREDIT = "https://credit.niso.org/contributor-roles/"
CRATE = "ro-crate-metadata.json"
MANIFEST = os.path.join("pipeline_info", "run_manifest.json")


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def iso_mtime(path):
    """Modification time as a UTC ISO 8601 instant, seconds resolution."""
    ts = datetime.datetime.fromtimestamp(os.path.getmtime(path),
                                         datetime.timezone.utc)
    return ts.replace(microsecond=0).isoformat().replace("+00:00", "Z")


class Crate:
    """The crate's @graph, indexed by @id, with helpers that overwrite rather
    than append so a second run changes nothing."""

    def __init__(self, doc):
        self.doc = doc
        self.graph = doc["@graph"]
        self.by_id = {e["@id"]: e for e in self.graph}

    def entity(self, eid, etype=None):
        e = self.by_id.get(eid)
        if e is None:
            e = {"@id": eid}
            if etype:
                e["@type"] = etype
            self.graph.append(e)
            self.by_id[eid] = e
        elif etype and "@type" not in e:
            e["@type"] = etype
        return e

    def root(self):
        return self.by_id["./"]

    def unlink(self, entity, prop, target_id):
        cur = entity.get(prop)
        if isinstance(cur, dict) and cur.get("@id") == target_id:
            del entity[prop]
        elif isinstance(cur, list):
            kept = [x for x in cur
                    if not (isinstance(x, dict) and x.get("@id") == target_id)]
            if len(kept) == 1:
                entity[prop] = kept[0]
            elif kept:
                entity[prop] = kept
            else:
                del entity[prop]

    def link(self, entity, prop, target_id):
        """Add a reference, keeping the property a list only when it needs to
        be, and never duplicating."""
        ref = {"@id": target_id}
        cur = entity.get(prop)
        if cur is None:
            entity[prop] = ref
        elif isinstance(cur, dict):
            if cur.get("@id") != target_id:
                entity[prop] = [cur, ref]
        elif isinstance(cur, list):
            if not any(x.get("@id") == target_id for x in cur if isinstance(x, dict)):
                cur.append(ref)

    def extend_context(self, term, uri):
        ctx = self.doc.get("@context")
        if isinstance(ctx, str):
            ctx = [ctx]
        extra = next((c for c in ctx if isinstance(c, dict)), None)
        if extra is None:
            extra = {}
            ctx.append(extra)
        extra[term] = uri
        self.doc["@context"] = ctx


def integrity(entity, path):
    entity["sha256"] = sha256_of(path)
    entity["contentSize"] = str(os.path.getsize(path))
    entity["dateModified"] = iso_mtime(path)


def add_integrity(crate, outdir):
    """Integrity fields on every File entity that resolves to a real file."""
    done = 0
    for e in list(crate.graph):
        t = e.get("@type")
        t = t if isinstance(t, list) else [t]
        if "File" not in t:
            continue
        eid = e["@id"]
        if eid.startswith(("http://", "https://", "file:", "#")):
            continue
        p = os.path.join(outdir, eid)
        if os.path.isfile(p):
            integrity(e, p)
            done += 1
    return done


def expand_directories(crate, outdir):
    """hasPart for the directory-valued outputs, which nf-prov leaves bare.

    Each file inside gets its own File entity with integrity fields, so the
    result tables are described rather than hidden behind a folder name.
    """
    expanded, added = 0, 0
    for e in list(crate.graph):
        if e.get("@type") != "Dataset" or e["@id"] == "./":
            continue
        eid = e["@id"].rstrip("/")
        if eid.startswith(("http://", "https://", "file:", "#")):
            continue
        d = os.path.join(outdir, eid)
        if not os.path.isdir(d):
            continue
        parts = []
        for root, _, names in os.walk(d):
            for n in sorted(names):
                full = os.path.join(root, n)
                rel = os.path.relpath(full, outdir).replace(os.sep, "/")
                child = crate.entity(rel, "File")
                child.setdefault("name", n)
                integrity(child, full)
                parts.append({"@id": rel})
                added += 1
        if parts:
            e["hasPart"] = parts
            expanded += 1
    return expanded, added


def add_people(crate, attribution):
    """Contributors as Person entities behind Role entities carrying CRediT
    URIs, plus affiliation and facility organisations.

    The Role indirection is the schema.org pattern: the property points at a
    Role, which repeats the property to reach the Person. It is how a role is
    attached to a person without inventing a vocabulary, and the CRediT URI is
    what makes the role machine-readable rather than a label.
    """
    root = crate.root()
    people = 0
    for i, c in enumerate(attribution.get("contributors") or []):
        name = c.get("name")
        if not name:
            continue
        # Someone recorded only because they happened to be the login gets an
        # entity so the run is attributable, but no roles: see run_manifest.
        pid = c.get("orcid") or "#contributor-%d" % (i + 1)
        person = crate.entity(pid, "Person")
        person["name"] = name
        people += 1

        if c.get("affiliation", "").startswith("https://ror.org/"):
            org = crate.entity(c["affiliation"], "Organization")
            crate.link(person, "affiliation", c["affiliation"])
            org.setdefault("name", c["affiliation"])

        if c.get("facility"):
            fid = (c.get("affiliation") or "#facility") + "#" + \
                  c["facility"].lower().replace(" ", "-")
            fac = crate.entity(fid, "Organization")
            fac["name"] = c["facility"]
            if c.get("facility_id"):
                fac["identifier"] = c["facility_id"]
            if c.get("affiliation", "").startswith("https://ror.org/"):
                crate.link(fac, "parentOrganization", c["affiliation"])
            crate.link(person, "memberOf", fid)

        tokens = c.get("role_tokens") or []
        if not tokens:
            crate.link(root, "author", pid)
            continue
        # nf-prov's own agent is already on the root as a direct author. Once
        # that person has Role entities, leaving the direct reference too makes
        # root.author a mix of Person and Role, which a consumer expecting
        # authors has to pick apart. The Role indirection reaches the Person
        # anyway, so the direct reference goes.
        crate.unlink(root, "author", pid)
        for tok in tokens:
            rid = "#role-%d-%s" % (i + 1, tok)
            role = crate.entity(rid, "Role")
            role["roleName"] = {"@id": CREDIT + tok + "/"}
            role["name"] = tok.replace("-", " ")
            crate.link(role, "author", pid)
            crate.link(root, "author", rid)
    return people


def add_project(crate, project):
    """The research activity, its funder, its grant and its rights holder."""
    if not project:
        return False
    root = crate.root()
    pid = project.get("raid") or "#project"
    proj = crate.entity(pid, "ResearchProject")
    if project.get("name"):
        proj["name"] = project["name"]
    if project.get("raid"):
        proj["identifier"] = project["raid"]
    # The project is the organisation on whose behalf the work was done, which
    # is what sourceOrganization means. ResearchProject is an Organization
    # subtype, so the reference is type-correct.
    crate.link(root, "sourceOrganization", pid)

    funder = project.get("funder") or {}
    if funder.get("ror") or funder.get("name"):
        fid = funder.get("ror") or "#funder"
        f = crate.entity(fid, "Organization")
        if funder.get("name"):
            f["name"] = funder["name"]
        # RO-Crate's own guidance puts funder on the root as well as on the
        # project, so a reader does not have to traverse to find it.
        crate.link(root, "funder", fid)
        crate.link(proj, "funder", fid)
        if project.get("grant"):
            gid = "#grant-" + project["grant"].replace("/", "-")
            g = crate.entity(gid, "Grant")
            g["identifier"] = project["grant"]
            crate.link(g, "funder", fid)
            crate.link(proj, "funding", gid)

    holder = project.get("rights_holder") or {}
    hid = holder.get("orcid") or holder.get("ror")
    if hid:
        h = crate.entity(hid, "Person" if holder.get("orcid") else "Organization")
        if holder.get("name"):
            h["name"] = holder["name"]
        crate.link(root, "copyrightHolder", hid)
    elif holder.get("name"):
        h = crate.entity("#rights-holder", "Organization")
        h["name"] = holder["name"]
        crate.link(root, "copyrightHolder", "#rights-holder")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("outdir", help="the run's --outdir, holding " + CRATE)
    ap.add_argument("--crate", help="override the crate path")
    ap.add_argument("--manifest", help="override the run manifest path")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change and write nothing")
    a = ap.parse_args()

    crate_path = a.crate or os.path.join(a.outdir, CRATE)
    manifest_path = a.manifest or os.path.join(a.outdir, MANIFEST)
    for p in (crate_path, manifest_path):
        if not os.path.isfile(p):
            sys.exit("not found: %s" % p)

    doc = json.load(open(crate_path, encoding="utf-8"))
    manifest = json.load(open(manifest_path, encoding="utf-8"))
    crate = Crate(doc)
    crate.extend_context("sha256", SPDX_CHECKSUM)

    # The licence and the citation come from the manifest too, so the crate
    # has one source rather than repeating them in the plugin's own config
    # where they could drift.
    attribution = manifest.get("attribution") or {}
    root = crate.root()
    if attribution.get("license"):
        root["license"] = "https://spdx.org/licenses/" + attribution["license"]
    if attribution.get("doi"):
        root["identifier"] = "https://doi.org/" + attribution["doi"]

    hashed = add_integrity(crate, a.outdir)
    expanded, added = expand_directories(crate, a.outdir)
    people = add_people(crate, attribution)
    project = add_project(crate, manifest.get("project"))

    print("crate:     %s" % crate_path)
    print("hashed:    %d existing file entities" % hashed)
    print("expanded:  %d directory outputs, %d files described" % (expanded, added))
    print("people:    %d contributors" % people)
    print("project:   %s" % ("added" if project else "none declared"))
    print("entities:  %d" % len(crate.graph))

    if a.dry_run:
        print("\ndry run, nothing written")
        return
    with open(crate_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(doc, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    print("written")


if __name__ == "__main__":
    main()
