#!/bin/sh
# Stage the trusted example inputs at a fresh path and render without execution.
set -eu
if [ "$#" -ne 1 ]; then
    printf '%s\n' "Usage: sh render.sh /path/to/new-build-directory" >&2
    exit 2
fi

# Resolve the tools and source before creating the output directory.
command -v quarto >/dev/null 2>&1 || {
    printf '%s\n' "Quarto is required on PATH." >&2
    exit 2
}
source_root=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd -P)
destination=$1
case $destination in
    /*) ;;
    *) destination=$PWD/$destination ;;
esac
destination_parent=$(CDPATH='' cd -- "$(dirname -- "$destination")" && pwd -P)
case $destination_parent/ in
    "$source_root"/|"$source_root"/*) printf '%s\n' "Use a build directory outside the source template." >&2; exit 2 ;;
esac
destination=$destination_parent/$(basename -- "$destination")

# mkdir refuses to replace an existing build. The project config pins render targets.
mkdir -- "$destination"
cp "$source_root/_quarto.yml" "$source_root/report.qmd" "$source_root/slides.qmd" \
   "$source_root/slide-navigation.html" "$source_root/title-slide.html" "$source_root"/_*.md "$destination"
for folder in data figures models; do
    cp -R "$source_root/$folder" "$destination/"
done
cd -- "$destination"
quarto render --no-execute
printf '%s\n' "$destination/site/report.html" "$destination/site/slides.html"
