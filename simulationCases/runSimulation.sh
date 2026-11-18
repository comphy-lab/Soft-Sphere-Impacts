#!/usr/bin/env bash
set -euo pipefail

case "${1:-}" in
    bounce|bounce-ve|bounce-kv)
        executable=$1
        source_name=$1.c
        flags=()
        ;;
    bounce-ve-wi)
        executable=$1
        source_name=bounce-ve.c
        flags=(-DBOUNCE_VE_EL_WI)
        ;;
    *)
        echo "Usage: $0 {bounce We | bounce-ve We Ec De | bounce-ve-wi We El Wi | bounce-kv We Ec}" >&2
        exit 2
        ;;
esac
shift

if [[ "${OPENMP:-0}" == 1 ]]; then
    flags+=(-fopenmp)
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
fi

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
# shellcheck source=src-local/build.sh
source "$source_dir/../src-local/build.sh"
build_basilisk "$source_dir/$source_name" "$executable" "${flags[@]+"${flags[@]}"}"
exec "./$executable" "$@"
