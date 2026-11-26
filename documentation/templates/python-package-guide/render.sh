#!/bin/sh
# Render a selected Markdown guide with the static Pandoc reader layout.
set -eu

if [ "$#" -ne 2 ]; then
  printf 'Usage: sh render.sh SOURCE.md OUTPUT.html\n' >&2
  exit 2
fi

layout=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd -P)
source_dir=$(CDPATH='' cd -- "$(dirname -- "$1")" && pwd -P)
output_dir=$(CDPATH='' cd -- "$(dirname -- "$2")" && pwd -P)
source_name=$(basename -- "$1")
output_name=$(basename -- "$2")

# A neighboring output keeps source, PDFs and other relative links portable.
case "$source_name" in
  *[!a-zA-Z0-9._-]*|.*|-*) printf 'Use a simple Markdown filename.\n' >&2; exit 2 ;;
  *.md) ;;
  *) printf 'The selected source must be Markdown.\n' >&2; exit 2 ;;
esac
case "$output_name" in
  *.html) ;;
  *) printf 'The output name must end in .html.\n' >&2; exit 2 ;;
esac
if [ "$source_dir" != "$output_dir" ] || [ ! -f "$source_dir/$source_name" ]; then
  printf 'Select an existing Markdown source and place the HTML beside it.\n' >&2
  exit 2
fi

cd -- "$source_dir"
html=$("${QUARTO:-quarto}" pandoc --sandbox --standalone \
  --from=markdown-raw_html-raw_attribute --to=html5 --section-divs \
  --toc --toc-depth=2 --no-highlight \
  --template="$layout/reader.html" --metadata="source-file:$source_name" \
  "./$source_name")
set -C
printf '%s\n' "$html" > "./$output_name"
