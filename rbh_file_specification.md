# Can you add this stuff to my repo with the specification for `rbh2`?
# This is the link to the repo: https://github.com/conchoecia/rbh-specs
# It would be better to just have one version of the file specification

# Overview
This file format is designed to store about conserved sequences across different species, including their locations on chromosomes (start and end points), species, scaffold and their orientation (strand +/-). This file format is second version of original .rbh format "LINK TO .rbh SPECIFICATION"
# File Format
- **Type:** Text

# File Structure

- **Header Line:**  first line of the file is a header. It consists of the following fields:
    - `Speciesname_gene`: Identifier for the alignment block or gene region (e.g., block_1).
    - `Speciesname_colour`: ???.
    - `Speciesname_scaf`: Scaffold identifier where the sequence is located.
    - `Speciesname_start`: Starting position of the conserved sequence on the chromosome.
    - `Speciesname_stop`: Ending position of the conserved sequence on the chromosome.
    - `Speciesname_strand`: Strand orientation, indicated by "+" for forward and "-" for reverse strands.
- **Data Lines:** Following the header, each line corresponds to a specific conserved sequence.

Sample File: 
```
Speciesname_gene Speciesname_colour Speciesname_scaf Speciesname_start Speciesname_stop Speciesname_strand
block_1 blue scaffold_23 12345 67890 +
block_2 red scaffold_11 98765 43210 -
```

For parsing:
- Each line in the file, after the header, represents a unique conserved sequence entry.
- Fields are separated by a single space or tab.

# Versioning
Actual version is 1.0. from April 2024. 
