#!/usr/bin/env bash

# Compile in a temporary directory with a local source filename for qcc.
build_basilisk () (
    set -euo pipefail
    local source_file=$1 output_name=$2
    shift 2
    local repo_root run_dir expected_ref installed_ref build_dir qcc_path
    repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
    run_dir=$(pwd -P)
    case "$run_dir/" in
        "$repo_root/"*)
            echo "Run from an output directory outside the repository." >&2
            exit 1 ;;
    esac

    if [[ ! -r "$repo_root/.project_config" || ! -r "$repo_root/basilisk/.comphy-lock" ]]; then
        echo "Install the Basilisk release in basilisk.ref first." >&2
        exit 1
    fi
    expected_ref=$(cat "$repo_root/basilisk.ref")
    installed_ref=$(sed -n 's/^ref=//p' "$repo_root/basilisk/.comphy-lock")
    if [[ "$installed_ref" != "$expected_ref" ]]; then
        echo "Basilisk ref mismatch. Expected $expected_ref, found $installed_ref." >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    source "$repo_root/.project_config"
    qcc_path=$(command -v qcc || true)
    if [[ ! "$BASILISK" -ef "$repo_root/basilisk/src" ||
          ! "$qcc_path" -ef "$repo_root/basilisk/src/qcc" ]]; then
        echo "qcc must resolve to the project-local Basilisk installation." >&2
        exit 1
    fi

    if [[ -L "$output_name" || -d "$output_name" ]]; then
        printf 'Cannot replace %s with a compiled executable.\n' "$output_name" >&2
        exit 1
    fi
    build_dir=$(mktemp -d .build.XXXXXX)
    trap 'rm -rf -- "$build_dir"' EXIT
    cp -- "$source_file" "$repo_root"/src-local/*.h "$build_dir/"
    local compile_source=${source_file##*/}
    if ! (cd "$build_dir" && "$qcc_path" -O2 -Wall -disable-dimensions "$@" "$compile_source" -o "$output_name" -lm) \
        > "$build_dir/compile.log" 2>&1 ||
        grep -qi 'parse error' "$build_dir/compile.log"; then
        cat "$build_dir/compile.log" >&2
        exit 1
    fi
    cat "$build_dir/compile.log" >&2
    mv -f -- "$build_dir/$output_name" "./$output_name"
)
