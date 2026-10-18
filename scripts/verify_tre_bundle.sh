#!/usr/bin/env bash
#
# Check an offline bundle is complete and self-consistent before it is uploaded.
#
#   scripts/verify_tre_bundle.sh <bundle-dir> [tarball]
#   scripts/verify_tre_bundle.sh ~/singularity-1.7.0 ~/rnaseq-flow-1.7.0-singularity.tar
#
# Inside a TRE there is no way to fetch a missing image, so every gap has to be
# found here. Checks that each container the modules declare has an image, that
# every name tre.config invents resolves to a file, that the plugin is present,
# and that the tarball matches its checksum.

set -uo pipefail

BUNDLE="${1:-}"
TARBALL="${2:-}"
[ -d "$BUNDLE" ] || { echo "usage: $0 <bundle-dir> [tarball]" >&2; exit 2; }
REPO="$(cd "$(dirname "$0")/.." && pwd)"

fail=0
note() { printf '  %-7s %s\n' "$1" "$2"; [ "$1" = "MISSING" ] && fail=1; return 0; }

echo "== containers declared in modules/local =="
while read -r uri; do
    img="$BUNDLE/$(echo "$uri" | tr '/:' '--').img"
    [ -f "$img" ] && note ok "$(basename "$img")" || note MISSING "$uri"
done < <(grep -hoE "^[[:space:]]*container[[:space:]]+'[^']+'" "$REPO"/modules/local/*.nf \
         | sed -E "s/.*'([^']+)'.*/\1/" | sort -u)

echo "== names tre.config points at =="
if [ -f "$BUNDLE/tre.config" ]; then
    while read -r name; do
        img="$BUNDLE/$(echo "$name" | tr '/:' '--').img"
        [ -f "$img" ] && note ok "$(basename "$img")" || note MISSING "$name"
    done < <(grep -oE "container = '[^']+'" "$BUNDLE/tre.config" | sed -E "s/.*'([^']+)'.*/\1/")
else
    note MISSING "tre.config"
fi

echo "== supporting files =="
[ -f "$BUNDLE/README.txt" ] && note ok README.txt || note MISSING README.txt
plug=$(ls -d "$BUNDLE"/plugins/nf-prov-* 2>/dev/null | head -1)
[ -n "$plug" ] && note ok "plugins/$(basename "$plug")" || note MISSING "plugins/nf-prov-*"

echo "== counts =="
echo "  images in bundle: $(ls "$BUNDLE"/*.img 2>/dev/null | wc -l)"
echo "  declared in repo: $(grep -hoE "^[[:space:]]*container[[:space:]]+'[^']+'" "$REPO"/modules/local/*.nf | sort -u | wc -l)"
echo "  named in tre.config: $(grep -c "container = " "$BUNDLE/tre.config" 2>/dev/null || echo 0)"

if [ -n "$TARBALL" ]; then
    echo "== tarball =="
    if [ -f "$TARBALL" ] && [ -f "$TARBALL.sha256" ]; then
        if ( cd "$(dirname "$TARBALL")" && sha256sum -c "$(basename "$TARBALL").sha256" >/dev/null 2>&1 ); then
            note ok "checksum matches"
        else
            note MISSING "checksum does NOT match"
        fi
    else
        note MISSING "tarball or .sha256 not found"
    fi
fi

echo
[ "$fail" -eq 0 ] && echo "bundle looks complete" || echo "BUNDLE INCOMPLETE, see MISSING above"
exit "$fail"
