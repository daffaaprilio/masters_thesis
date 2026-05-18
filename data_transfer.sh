#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] <src> <dst>

Rsync files between servers via SSH. Server aliases are resolved via ~/.ssh/config.

Arguments:
  src       Source path (use 'host:/path' for remote)
  dst       Destination path (use 'host:/path' for remote)

Options:
  -n, --dry-run    Show what would be transferred without doing it
  -h, --help       Show this help message

Examples:
  # Push local dir to remote
  $(basename "$0") /local/data/ myserver:/remote/data/

  # Pull from remote to local
  $(basename "$0") myserver:/remote/data/ /local/data/

  # Dry run
  $(basename "$0") --dry-run /local/data/ myserver:/remote/data/
EOF
}

DRY_RUN=false
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

if [[ ${#POSITIONAL[@]} -ne 2 ]]; then
    echo "Error: expected 2 positional arguments, got ${#POSITIONAL[@]}" >&2
    usage >&2
    exit 1
fi

SRC="${POSITIONAL[0]}"
DST="${POSITIONAL[1]}"

RSYNC_OPTS=(-avz --progress -e "ssh")

if $DRY_RUN; then
    RSYNC_OPTS+=(--dry-run)
    echo "[dry-run] Would rsync: $SRC -> $DST"
else
    echo "Syncing: $SRC -> $DST"
fi

rsync "${RSYNC_OPTS[@]}" "$SRC" "$DST"
