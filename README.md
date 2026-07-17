# ConservFinder
The ConservFinder is a script designed to identify conserved sequences of a certain length that are near same (occurrence is above certain threshold) within alignment block of MAF file. The aim of this script is to assist in studies of genomic conservation across different species.

## Features
**Input**: MAF file, names of species of interest, threshold.
**Output**: Number of alignment blocks, number of conserved sequences and an .rbh2 file.

**Used packages**:
- collections.Counter: For counting nucleotide occurrences.
- Bio.AlignIO: For parsing MAF files.
- pandas: For writing the .rbh2 table.
- argparse: For interacting with script trough terminal.

## Functions

- **NtCounter** -> Calculates conserved nucleotides based on a threshold.
- **Indices_to_Ranges** -> Converts matching indices to conserved sequence ranges.
- **maf_coords** -> Maps alignment ranges to genomic coordinates.
- **process_alignments** -> Reads the MAF and writes the .rbh2 file.
- **Main**

## Usage

Run the script with the following commands:
- **path** to your MAF file (-f),
- **threshold** for conserved nucleotides (-t). The default is 0.6.
- a list of **species/record.id** as in your .maf file (-s),
- a **name** of output file in .rbh2 format that will be saved in current directory (-o).
- the minimum conserved sequence length (-m). The default is 10.

Species can be given as `Species` or `Species.scaffold`.

Example:
```python
python conservfinder.py -f "path/to/maf_file.maf" -t 0.7 -s "Species1" "Species2" "Species3" -m 10 -o "conserved_regions.rbh2"
```

Install packages with:
```bash
pip install -r requirements.txt
```
