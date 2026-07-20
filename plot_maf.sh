#!/bin/sh

# Use:
# plot_maf.sh MAF BED BEDGRAPH REF REGION "SPECIES" "LABELS" OUTPUT
# - BED and BEDGRAPH are made by ConservFinder.
# - REF is the first -s species.
# - REGION is scaffold:start-end.
# - SPECIES and LABELS stay in quotes.

maf=$1
bed=$2
score=$3
reference=$4
region=$5
species=$6
labels=$7
output=$8

tracks="${output}.ini"
index="${output}.index"

printf '[agreement]\n' > "$tracks"
printf 'file = %s\n' "$score" >> "$tracks"
printf 'file_type = bedgraph\n' >> "$tracks"
printf 'title = Agreement\n' >> "$tracks"
printf 'height = 2\n' >> "$tracks"
printf 'color = #55599b\n' >> "$tracks"
printf 'min_value = 0\n' >> "$tracks"
printf 'max_value = 1\n\n' >> "$tracks"

printf '[conserved]\n' >> "$tracks"
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
