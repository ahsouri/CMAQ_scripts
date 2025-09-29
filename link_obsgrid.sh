#!/bin/bash
# link_obsgrid_files_to_cwd.sh
# Usage: ./link_obsgrid_files_to_cwd.sh /path/to/your/files

indir="${1:-.}"    # directory containing original files, default to current dir
outdir="$(pwd)"    # save links in the current working directory

cd "$indir" || { echo "Input directory not found: $indir"; exit 1; }

for f in SURFACE_OBS:??????????; do
    # Skip if no match
    [ -e "$f" ] || continue

    base="${f%%:*}"             # should be "SURFACE_OBS"
    rest="${f#*:}"              # YYYYMMDDHH part

    yyyy="${rest:0:4}"
    mm="${rest:4:2}"
    dd="${rest:6:2}"
    hh="${rest:8:2}"

    newname="${base}:${yyyy}-${mm}-${dd}_${hh}"

    # Create symlink in the current directory
    ln -sf "$indir/$f" "$outdir/$newname"
    echo "Linked $indir/$f -> $outdir/$newname"
done
