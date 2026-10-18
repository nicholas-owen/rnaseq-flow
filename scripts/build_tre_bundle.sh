#!/usr/bin/env bash
#
# Build the offline Singularity bundle for a TRE or air-gapped cluster.
#
#   scripts/build_tre_bundle.sh <version> [work-dir] [destination]
#   scripts/build_tre_bundle.sh 1.7.0
#
# Produces, in the destination directory (default $HOME):
#
#   singularity-<version>/                 the bundle directory
#     *.img                                every container the pipeline uses
#     plugins/nf-prov-<ver>/               the plugin, which Nextflow resolves
#                                          on every run and cannot download here
#     tre.config                           defines -profile tre
#     README.txt
#   rnaseq-flow-<version>-singularity.tar
#   rnaseq-flow-<version>-singularity.tar.sha256
#
# Two kinds of image go in, and they are found two different ways.
#
# The containers declared in modules/local/*.nf are pulled straight from their
# registry. The seven processes that declare a `conda` directive instead have
# no registry image at all: Wave builds one per run from the Conda spec, under
# an ephemeral URL. Those are taken from a completed run's task directories,
# which record the image in `.command.run`, and rewritten to stable local names
# that tre.config then points at.
#
# So this needs a run that exercised all seven first. The script says which are
# missing rather than producing a bundle with holes in it.

set -euo pipefail

VERSION="${1:-}"
WORK="${2:-work}"
[ -n "$VERSION" ] || { echo "usage: $0 <version> [work-dir] [destination]" >&2; exit 2; }

REPO="$(cd "$(dirname "$0")/.." && pwd)"
# Third argument, else $HOME. Not $REPO/.. : a repository cloned inside another
# directory of the same name puts the bundle somewhere surprising, and it
# should not land inside the checkout either.
DEST="${3:-$HOME}"
OUT="$DEST/singularity-$VERSION"
TAR="$DEST/rnaseq-flow-$VERSION-singularity.tar"

# The seven Conda-only processes, and the local image name each gets.
CONDA_PROCS="DESEQ2:deseq2 DIFFSPLICE:diffsplice DOWNLOAD_GMT:download-gmt \
DTU:dtu GPROFILER:gprofiler GSEA:gsea QUARTO_REPORT:quarto-report"

command -v singularity >/dev/null || { echo "singularity not found" >&2; exit 1; }

mkdir -p "$OUT/plugins"
cd "$REPO"

# ---- 1. the declared containers ------------------------------------------
echo "== registry containers =="
grep -hoE "^\s*container\s+'[^']+'" modules/local/*.nf \
  | sed -E "s/.*'([^']+)'.*/\1/" | sort -u | while read -r uri; do
    img="$OUT/$(echo "$uri" | tr '/:' '--').img"
    if [ -f "$img" ]; then
        echo "  have  $(basename "$img")"
    else
        echo "  pull  $uri"
        singularity build --force "$img" "docker://$uri" >/dev/null 2>&1 \
          || { echo "    FAILED to pull $uri" >&2; exit 1; }
    fi
done

# ---- 2. the Wave-built images --------------------------------------------
echo "== Wave-built images, from a completed run =="
missing=""
for pair in $CONDA_PROCS; do
    proc="${pair%%:*}"; name="${pair##*:}"
    # Newest task directory for this process wins, so a rerun after a spec
    # change supersedes the image built before it.
    # Match the process under any workflow prefix: DOWNLOAD_GMT runs as
    # DOWNLOAD:DOWNLOAD_GMT, not RNASEQ:DOWNLOAD_GMT, and a tagged task carries
    # a trailing " (sample)".
    run=$(grep -lsE "^### name: '[^']*:${proc}[ ']" "$WORK"/*/*/.command.run 2>/dev/null \
          | xargs -r ls -t 2>/dev/null | head -1 || true)
    if [ -z "$run" ]; then
        missing="$missing $proc"
        continue
    fi
    uri=$(sed -n "s/^### container: '\(.*\)'$/\1/p" "$run" | head -1)
    if [ -z "$uri" ]; then
        missing="$missing $proc"
        continue
    fi
    img="$OUT/rnaseq-flow-$name-$VERSION.img"

    # A Wave URL is per-run and does not stay pullable, so the one recorded in
    # .command.run is often already gone. The tag survives: its trailing hash
    # identifies the image contents, and every run that needed the same
    # environment pulled an image carrying it under a different URL. So if the
    # exact URL is not in the local Docker store, take any local image with the
    # same tag.
    tag="${uri##*:}"
    if docker image inspect "$uri" >/dev/null 2>&1; then
        local_img="$uri"
    else
        local_img=$(docker images --format '{{.Repository}}:{{.Tag}}' 2>/dev/null | grep -F ":$tag" | head -1 || true)
    fi

    echo "  $proc -> $(basename "$img")"
    if [ -n "$local_img" ]; then
        # Via docker save, not docker-daemon://. Singularity cannot read the
        # layers of some Wave images straight from the daemon and dies on a
        # blob it cannot find, while the same image runs fine under Docker. An
        # archive sidesteps that and worked for every image here.
        tmptar=$(mktemp "${TMPDIR:-/tmp}/rnaseq-flow-bundle.XXXXXX")
        if ! docker save -o "$tmptar" "$local_img" 2>/dev/null; then
            rm -f "$tmptar"; echo "    FAILED: could not export $local_img" >&2; exit 1
        fi
        if ! singularity build --force "$img" "docker-archive://$tmptar" >/dev/null 2>&1; then
            rm -f "$tmptar"
            echo "    FAILED to convert $local_img" >&2
            echo "    To see why, run it with the output shown:" >&2
            echo "      docker save -o /tmp/x.tar '"'"'$local_img'"'"'" >&2
            echo "      singularity build /tmp/x.img docker-archive:///tmp/x.tar" >&2
            exit 1
        fi
        rm -f "$tmptar"
    else
        # Nothing local carries this content, and a Wave URL is per-run, so
        # this rarely succeeds.
        singularity build --force "$img" "docker://$uri" >/dev/null 2>&1           || { echo "    FAILED: no local image carries tag $tag, and the" >&2
               echo "    Wave URL no longer resolves. Re-run the pipeline so" >&2
               echo "    $proc executes, then retry." >&2
               exit 1; }
    fi
done

if [ -n "$missing" ]; then
    echo >&2
    echo "No completed task found for:$missing" >&2
    echo "Run the pipeline so every Conda-built process executes, then retry." >&2
    echo "DOWNLOAD_GMT needs --download_refs --download_gmt; DTU and DIFFSPLICE" >&2
    echo "need a pseudo-aligner route with --dtu --diffsplice." >&2
    exit 1
fi

# ---- 3. the plugin --------------------------------------------------------
echo "== plugin =="
PLUGIN_DIR="${NXF_HOME:-$HOME/.nextflow}/plugins"
plugin=$(ls -d "$PLUGIN_DIR"/nf-prov-* 2>/dev/null | sort -V | tail -1 || true)
[ -n "$plugin" ] || { echo "nf-prov not in $PLUGIN_DIR; run the pipeline once" >&2; exit 1; }
cp -r "$plugin" "$OUT/plugins/"
echo "  $(basename "$plugin")"

# ---- 4. tre.config --------------------------------------------------------
{
  echo "profiles {"
  echo "  tre {"
  echo "    wave.enabled = false"
  echo "    singularity.enabled = true"
  echo "    singularity.autoMounts = true"
  echo "    process {"
  for pair in $CONDA_PROCS; do
      printf "      withName: '%s' { container = 'rnaseq-flow/%s:%s' }\n" \
             "${pair%%:*}" "${pair##*:}" "$VERSION"
  done
  echo "    }"
  echo "  }"
  echo "}"
} > "$OUT/tre.config"

# ---- 5. README ------------------------------------------------------------
cat > "$OUT/README.txt" <<EOF
rnaseq-flow $VERSION, offline Singularity image bundle

$(ls "$OUT"/*.img | wc -l) images, for running without internet access (e.g. in a TRE).

Usage:
  tar -xf rnaseq-flow-$VERSION-singularity.tar
  export NXF_SINGULARITY_CACHEDIR=\$PWD/singularity-$VERSION
  cp -r \$NXF_SINGULARITY_CACHEDIR/plugins/nf-prov-* ~/.nextflow/plugins/
  export NXF_OFFLINE=true
  nextflow run main.nf -c \$NXF_SINGULARITY_CACHEDIR/tre.config -profile tre,<your-profile> ...

The plugin must go in ~/.nextflow/plugins. Pointing NXF_PLUGINS_DIR at this
bundle does not work: Nextflow reports "No plugins", carries on, and the
RO-Crate is silently missing. NXF_OFFLINE=true is required, not optional.

tre.config disables Wave and pins the seven Conda-built processes to the image
names in this bundle. Those names are local to it and resolve nowhere online.

Must match pipeline release $VERSION: container references and Conda specs both
changed in 1.7.0, so earlier bundles do not work with it.
EOF

# ---- 6. pack --------------------------------------------------------------
echo "== packing =="
tar -cf "$TAR" -C "$(dirname "$OUT")" "$(basename "$OUT")"
( cd "$(dirname "$TAR")" && sha256sum "$(basename "$TAR")" > "$TAR.sha256" )

echo
echo "bundle:  $OUT  ($(du -sh "$OUT" | cut -f1), $(ls "$OUT"/*.img | wc -l) images)"
echo "tarball: $TAR  ($(du -h "$TAR" | cut -f1))"
echo "sha256:  $TAR.sha256"
