#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 Ldomain GridsPerR We tsnap [output.mp4]" >&2
    exit 2
fi

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
# shellcheck source=src-local/build.sh
source "$source_dir/../src-local/build.sh"
build_basilisk "$source_dir/getFacet.c" getFacet
build_basilisk "$source_dir/getData.c" getData
python3 "$source_dir/VideoFullDomain.py" "$1" "$2" "$3" "$4"
ffmpeg -framerate 30 -pattern_type glob -i 'Video_v1/*.png' -vf scale=850:880 \
    -c:v libx264 -r 30 -pix_fmt yuv420p -y -- "${5:-vid_$3.mp4}"
