#!/bin/sh
# Container entrypoint for the Hugging Face Space.
# The ~3.8 GB data tree is NOT in the image (a Space repo is capped at 1 GB); it is a
# private Hugging Face bucket mounted read-only at $ARTHOMIX_DATA_SRC (default /data).
# Link each top-level entry into /app/data, keeping .cache and uploads local and writable
# so user uploads and computed caches never reach the bucket.
set -eu

SRC="${ARTHOMIX_DATA_SRC:-/data}"

# The mount is lazy; give the platform a moment, then fail with a readable message.
i=0
while [ ! -d "$SRC/preloaded/transcriptomics" ] && [ "$i" -lt 60 ]; do
  i=$((i + 1)); sleep 2
done
if [ ! -d "$SRC/preloaded/transcriptomics" ]; then
  echo "ERROR: data volume not found at $SRC/preloaded/transcriptomics." >&2
  echo "Mount the bucket:  hf spaces volumes set <space> -v hf://buckets/<user>/arthomix-data:$SRC:ro" >&2
  exit 1
fi

mkdir -p /app/data /app/data/.cache /app/data/uploads
for entry in "$SRC"/* "$SRC"/.[!.]*; do
  [ -e "$entry" ] || continue
  name="$(basename "$entry")"
  case "$name" in .cache|uploads) continue ;; esac
  ln -sfn "$entry" "/app/data/$name"
done

# Explicit host/port override the 127.0.0.1:7788 dev options set in .Rprofile.
exec R -e "shiny::runApp(appDir = '/app', host = '0.0.0.0', port = 7860)"
