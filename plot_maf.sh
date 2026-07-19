#!/bin/sh

# Use:
# plot_maf.sh MAF BED REF REGION "SPECIES" "LABELS" OUTPUT
# - BED is made by ConservFinder.
# - REF is the first -s species.
# - REGION is scaffold:start-end.
# - SPECIES and LABELS stay in quotes.

maf=$1
bed=$2
reference=$3
region=$4
species=$5
labels=$6
output=$7

tracks="${output}.ini"
index="${output}.index"

printf '[conserved]\n' > "$tracks"
printf 'file = %s\n' "$bed" >> "$tracks"
printf 'file_type = bed\n' >> "$tracks"
printf 'title = ConservFinder\n' >> "$tracks"
printf 'height = 1\n' >> "$tracks"
printf 'color = #e8872d\n' >> "$tracks"
printf 'border_color = none\n' >> "$tracks"
printf 'display = collapsed\n' >> "$tracks"
printf 'labels = false\n\n' >> "$tracks"

printf '[maf]\n' >> "$tracks"
printf 'file = %s\n' "$maf" >> "$tracks"
printf 'file_index = %s\n' "$index" >> "$tracks"
printf 'reference = %s\n' "$reference" >> "$tracks"
printf 'species_order = %s\n' "$species" >> "$tracks"
printf 'species_labels = %s\n' "$labels" >> "$tracks"
printf 'species_order_only = true\n' >> "$tracks"
printf 'height = 4\n' >> "$tracks"
printf 'color_identical = black\n' >> "$tracks"
printf 'color_mismatch = grey\n' >> "$tracks"
printf 'color_gap = white\n' >> "$tracks"
printf 'rasterize = true\n' >> "$tracks"
printf 'title =\n\n' >> "$tracks"
printf '[x-axis]\n' >> "$tracks"

pyGenomeTracks --tracks "$tracks" --region "$region" --width 30 --fontSize 7 --dpi 150 -o "$output"
