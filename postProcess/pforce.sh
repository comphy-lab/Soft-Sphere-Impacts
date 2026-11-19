#!/usr/bin/env bash
set -euo pipefail

case "${1:-}" in
    bounce|bounce-ve) expected=2 ;;
    bounce-kv|bounce-ve-wi|bounce-bingham) expected=3 ;;
    *) expected=0 ;;
esac
if [[ $expected == 0 || $# -ne $expected ]]; then
    echo "Usage: $0 {bounce We | bounce-ve We | bounce-ve-wi We El | bounce-kv We Ec | bounce-bingham We J}" >&2
    exit 2
fi

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
# shellcheck source=src-local/build.sh
source "$source_dir/../src-local/build.sh"
build_basilisk "$source_dir/getBetaMax.c" getBetaMax

case "$1" in
    bounce|bounce-ve)
        python3 "$source_dir/getBetaMax.py" "$2"
        python3 "$source_dir/plot_force.py" "$2" "$1"
        python3 "$source_dir/f_vcm.py" "$2" "$1"
        ;;
    bounce-bingham)
        python3 "$source_dir/getBetaMax.py" "$2"
        python3 "$source_dir/plot_force.py" "$2" "$1" "$3"
        python3 "$source_dir/f_vcm.py" "$2" "$1" "$3"
        ;;
    bounce-kv)
        python3 "$source_dir/getBetaMax.py" "$2"
        python3 "$source_dir/plot_force_kv.py" "$2" "$3"
        python3 "$source_dir/f_vcm_kv.py" "$2" "$3"
        ;;
    bounce-ve-wi)
        python3 "$source_dir/getBetaMax_wi.py" "$3"
        python3 "$source_dir/plot_force_wi.py" "$2" "$3"
        python3 "$source_dir/f_vcm_wi.py" "$2" "$3"
        ;;
esac

case "$1" in
    bounce|bounce-ve|bounce-kv)
        build_basilisk "$source_dir/getRfoot.c" getRfoot
        python3 "$source_dir/getRfoot.py"
        ;;
esac

case "$1" in
    bounce) python3 "$source_dir/plot_Rfoot.py" "$2" ;;
    bounce-kv) python3 "$source_dir/plot_Rfoot_kv.py" "$2" "$3" ;;
esac
