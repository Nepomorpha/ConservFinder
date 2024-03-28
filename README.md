# ConservFinder
The ConservFinder is a script designed to identify conserved sequences of certain leenght that are near same (occuance is above certain threshold) within alighment block of MAF file. The aim of this script is to assist in studies of genomic conservation across different species.

## Features
**Input**: MAF file, names of species of interest, threshold.
**Output**: (1) Conserved sequences with genomic ranges both relative to the beginning of alignment and specific chromosome positions. (2) (Soon be made) .bed file with ranges. 

**Used packages**:
- collections.Counter: For counting nucleotide occurrences.
- Bio.AlignIO: For parsing MAF files.

## Functions

- **NtCounter** -> Calculates conserved nucleotides based on a threshold.
- **Indices_to_Ranges** -> Converts matching indices to conserved sequence ranges.
- **Ranges_to_Coordinates** -> Maps ranges of conserved sequences to their genomic coordinates. Write them in .bed file. 
- **Main**

## Usage

Run the script with the following commands:
- **path** to your MAF file (-f),
- **threshold** for conserved nucleotides (-t). The default is 0.6.
- a list of **species/record.id** as in your .maf file (-s).

Example:
```python
python conservfinder.py -f "path/to/maf_file.maf" -t 0.7 -s "Species1" "Species2" "Species2"
```
