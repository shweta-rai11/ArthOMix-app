#!/usr/bin/env bash
# Build the CODE-ONLY tree for the Hugging Face Docker Space. Local only: it never
# touches the network or the Hub.
#
#   deploy/huggingface/make_space_snapshot.sh [OUT_DIR]     (default ~/arthomix-hf-space-code)
#
# A Space repository is capped at 1 GB, so the ~3.8 GB data/ tree is NOT included:
# it is uploaded separately to a private Hugging Face bucket and mounted into the
# container at /data (see start.sh and README notes in DEPLOY steps).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP="$(cd "$HERE/../.." && pwd)"
OUT="${1:-$HOME/arthomix-hf-space-code}"

[ -f "$APP/global.R" ] || { echo "run from an ArthOMix checkout (global.R not found in $APP)" >&2; exit 1; }
if [ -e "$OUT" ] && [ -n "$(ls -A "$OUT" 2>/dev/null)" ]; then
  echo "refusing to overwrite non-empty $OUT - remove it or pass another path" >&2; exit 1
fi
mkdir -p "$OUT"

cd "$APP"
# Tracked files only (no caches/uploads), minus data and everything not needed at runtime.
git ls-files -z \
  | grep -zvE '^(data/|tests/|reproduce/|results/|deploy/|\.github/|test_results_thesis|docker-compose\.yml$|Dockerfile$|\.dockerignore$|\.gitattributes$)' \
  | rsync -a --from0 --files-from=- "$APP/" "$OUT/"

cp "$HERE/Dockerfile.space" "$OUT/Dockerfile"
cp "$HERE/README.space.md"  "$OUT/README.md"
cp "$HERE/start.sh"         "$OUT/start.sh"
chmod +x "$OUT/start.sh"
cat > "$OUT/.dockerignore" <<'EOF'
.git
.Rproj.user
.Rhistory
.RData
.env
data
EOF

# Hard check: the Space repo must stay far below its 1 GB cap.
SIZE_KB="$(du -sk "$OUT" | cut -f1)"
[ "$SIZE_KB" -lt 204800 ] || { echo "code snapshot is $((SIZE_KB/1024)) MB - too big for a Space repo, check what got copied" >&2; exit 1; }
BAD="$(find "$OUT" -type f -size +10M -print)"
[ -z "$BAD" ] || { echo "files >10MB in code snapshot (Hub rejects plain-git files this big):"; echo "$BAD"; exit 1; }

echo "Code snapshot ready: $OUT ($((SIZE_KB/1024)) MB, $(find "$OUT" -type f | wc -l | tr -d ' ') files)"
echo "Upload:  hf upload <user>/<space> $OUT . --repo-type space"
