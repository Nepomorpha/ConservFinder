# ConservFinder
The ConservFinder is a script designed to identify conserved sequences of a certain length that are near same (occurrence is above certain threshold) within alignment block of MAF file. The aim of this script is to assist in studies of genomic conservation across different species.

## Features
**Input**: MAF file, names of species of interest, threshold.
**Output**: (1) Number of alignment blocks and number of conserved sequences (2) .bed file with genomic ranges.

**Used packages**:
- collections.Counter: For counting nucleotide occurrences.
- Bio.AlignIO: For parsing MAF files.
- pybedtools: For writing .bed files.
- argparse: For interacting with script trough terminal.

## Functions

- **NtCounter** -> Calculates conserved nucleotides based on a threshold.
- **Indices_to_Ranges** -> Converts matching indices to conserved sequence ranges.
- **Ranges_to_Coordinates** -> Maps ranges of conserved sequences to their genomic coordinates. Write them in .bed file.
- **Main**

## Usage

Run the script with the following commands:
- **path** to your MAF file (-f),
- **threshold** for conserved nucleotides (-t). The default is 0.6.
- a list of **species/record.id** as in your .maf file (-s),
- a **name** of output file in .bed format that will be saved in current directory (-o).

Example:
```python
python conservfinder.py -f "path/to/maf_file.maf" -t 0.7 -s "Species1" "Species2" "Species2" -o "conserved_regions.bed"
```
