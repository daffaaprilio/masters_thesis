#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] <server> <src> <dst>

Rsync files between servers via SSH. Server aliases are resolved via ~/.ssh/config.

Arguments:
  server    SSH host alias (from ~/.ssh/config)
  src       Source path (prefix with 'server:' to pull from remote, omit to push)
  dst       Destination path

Options:
  -n, --dry-run    Show what would be transferred without doing it
  -h, --help       Show this help message

Examples:
  # Push local dir to remote
  $(basename "$0") myserver /local/data/ myserver:/remote/data/

  # Pull from remote to local
  $(basename "$0") myserver myserver:/remote/data/ /local/data/

  # Dry run
  $(basename "$0") --dry-run myserver /local/data/ myserver:/remote/data/
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

if [[ ${#POSITIONAL[@]} -ne 3 ]]; then
    echo "Error: expected 3 positional arguments, got ${#POSITIONAL[@]}" >&2
    usage >&2
    exit 1
fi

SERVER="${POSITIONAL[0]}"
SRC="${POSITIONAL[1]}"
DST="${POSITIONAL[2]}"

RSYNC_OPTS=(-avz --progress -e "ssh")

if $DRY_RUN; then
    RSYNC_OPTS+=(--dry-run)
    echo "[dry-run] Would rsync: $SRC -> $DST via $SERVER"
else
    echo "Syncing: $SRC -> $DST via $SERVER"
fi

rsync "${RSYNC_OPTS[@]}" "$SRC" "$DST"
