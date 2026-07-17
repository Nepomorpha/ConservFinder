#!/usr/bin/env python
"""

Program     : ConservFinder for .maf files.
Language    : Python
Author      : Yehor Tertyshnyk, Darrin T. Schultz
Email       : egor.tertyshnyk@gmail.com
Github      : https://github.com/Nepomorpha/ConservFinder
License     : ?--?
Citation    : -TBS-

Description : This program takes .maf file and output conserved sequences as well as their genomic range in .rbh2 format.

Usage
instructions: see -TBS-

"""
from collections import Counter
from Bio import AlignIO
import argparse
import pandas as pd


def get_args():
    """
    The purpose of this function is accepting parameter such as path to the file, threshold or included species as command-line arguments.

    Template for input is:

    python conservfinder.py -f "path_to_maf_file.maf" -t 0.7 -s "species1" "species2" "species3" -o "name_of_example_file.rbh2"

    Example usage:

    python conservfinder.py -f "/Users/egortertyshnyk/Desktop/Simakov_Group/Conserved_regions/MOLLUSC_Chr10_small.maf" -t 0.6 -s "Octopusvulgaris6645.OX597823.1" "Octopusbimaculoides37653.NC_068990.1" "Octopussinensis2607531.NC_043005.1" -o "example.rbh2"

    """
    parser = argparse.ArgumentParser(description="Process MAF files to find conserved regions.")
    parser.add_argument('-f', '--file',      required=True,               help="Path to the MAF file")
    parser.add_argument('-t', '--threshold', type=float, default=0.6,     help="Threshold for freq. of conserved nts")
    parser.add_argument('-s', '--species',   nargs='+', required=True,    help="List of species to include")
    parser.add_argument('-o', '--output',    default='output.rbh2',        help="Output filename for the rbh2 file")
    parser.add_argument('-m', '--min-length', type=int, default=10,       help="Minimum length of a conserved region")
    parser.add_argument('-a', '--aligner',   type=str, default="cactus",  help="Tells the program which aligner was used to generate the MAF file. Default is 'cactus'. No other options are currently supported")
    args = parser.parse_args()
    # Check if the aligner is supported
    if args.aligner not in ["cactus"]:
        raise ValueError(f"Aligner '{args.aligner}' is not supported. Only 'cactus' is supported.")
    return args

def NtCounter(sequences, threshold):
    """
    The purpose of this function:
        The purpose of this function is to find the conserved nucleotides in a multiple sequence alignment.
        Conserved nucleotides -> nucleotides that occur in single column of alignment more often than in previously specified threshold.
    Example:
          ix:              1111111111222222
                 01234567890123456789012345
         %id:    11111    111      11111111
                 000006666000      00000000
                 00000666600000000000000000
        >70%:    *****    ***      ********
        Seq1:    ATGCTGAGTCGT------GTATCGAT
        Seq2:    ATGCTCCCCCGTAAAAAAGTATCGAT
        Seq3:    ATGCTGAGTCGT------GTATCGAT

    Example output:
       The conserved indices in the above example are:
          [0, 1, 2, 3, 4, 9, 10, 11, 18, 19, 20, 21, 22, 23, 24, 25]
        Therefore, the list of integers like above is returned.

    The arguments are:
      - sequences: an iterable of strings, each string is a sequence. The permissible characters are "ACGNT-".
        - example: ['ATGCATGA', 'ACCCATGA', 'ATGCACGA']
      - threshold: a float, the minimum fraction of a nucleotide to be considered conserved. The value should be between 0 and 1.

    The output is:
      - conserved_indices: a list of integers, the indices of the conserved nucleotides
    """
    conserved_indices = []
    for sn in range(len(sequences[0])):
        column = ''.join(seq[sn] for seq in sequences)
        c = Counter(column)
        total = sum(c.values())
        for nucleotide, count in c.items():
            if nucleotide in "ACGT" and count / total >= threshold:
                conserved_indices.append(sn)
                break
    return conserved_indices

def indices_to_ranges(matching_indices, min_match_len = 10): # Note: it was changed from 5 to 10
    """
    Converts a list of indices to a list of start and end ranges.

    Example:
        In the example from the function NtCounter, the list of matching indices is:
            [0, 1, 2, 3, 4, 9, 10, 11, 18, 19, 20, 21, 22, 23, 24, 25]
        The ranges of conserved nucleotides are:
            [(0, 4), (9, 11), (18, 25)]

    Arguments:
      - matching_indices: List of indices (int) where nucleotides are conserved.
      - min_match_len: Minimum length of a conserved region to be considered.
                       A value of 5 means that a conserved region must be at least 5 nucleotides long.

    Returns:
      - List of tuples representing start and end indices of conserved regions.
        Because these are indices for the sequence positions, we must change them later for rbh2 format.
    """
    range_indices = []
    if not matching_indices:  # Check if the list is empty
        return range_indices  # Return an empty list if there are no matching indices

    start = end = matching_indices[0]
    for index in matching_indices[1:] + [None]:
        if index is not None and index == end + 1:
            end = index
        else:
            if end - start >= (min_match_len - 1): # For example, start = 1, end = 5, is a match of length 5
                range_indices.append((start, end))
            if index is not None:
                start = end = index
    return range_indices

def maf_coords(record, start, end):
    """
    Converts ranges of conserved sequences into genomic coordinates.
    Each alignment comes from a unique scaffold and position from the source genomes.
    Therefore, we must convert the indices of the alignment to genomic coordinates.

    Arguments:
      - record: Sequence record corresponding to one species.
                It contains the fasta header, chromosome start and strand.
                It is parsed from the .maf file.
      - start: Start index of the conserved region in the alignment.
      - end: End index of the conserved region in the alignment.
    """
    sequence = str(record.seq)
    before = 0
    for nucleotide in sequence[:start]:
        if nucleotide != '-':
            before += 1

    length = 0
    for nucleotide in sequence[start:end + 1]:
        if nucleotide != '-':
            length += 1

    maf_start = record.annotations['start']

    if record.annotations['strand'] == 1:
        genomic_start = maf_start + before
        genomic_end = genomic_start + length
    else:
        genomic_end = record.annotations['srcSize'] - maf_start - before
        genomic_start = genomic_end - length

    return genomic_start, genomic_end


def process_alignments(maf_file, species_list, aligner, threshold, output_file, min_length=10):
    """
    Process the MAF file to find conserved regions and save them to an RBH2 file.

    Arguments:
    - maf_file:     Path to the MAF file.
    - species_list: List of species to include in the analysis.
    - aligner:      The aligner used to generate the MAF file. Currently, only "cactus" is supported.
    - threshold:    Threshold for the frequency of conserved nucleotides.
    - output_file:  Output filename for the RBH2 file.

    Packages used:
        - Bio.AlignIO
        - pandas
    """
    if aligner != "cactus":
        raise ValueError("Aligner not supported. Only 'cactus' is supported.")

    species_names = []
    for species in species_list:
        species_names.append(species.split('.')[0])
    species_list = species_names

    maf_count = 0
    skip_count = 0
    rbh2_entries = []

    for multiple_alignment in AlignIO.parse(maf_file, "maf"):
        maf_count += 1
        filtered_records = []
        species_in_alignment = []

        for record in multiple_alignment:
            species_identifier = record.id.split('.')[0]
            if species_identifier in species_list:
                filtered_records.append(record)
                species_in_alignment.append(species_identifier)

        if len(filtered_records) == len(species_list) and len(species_in_alignment) == len(set(species_in_alignment)):
            sequences = []
            for record in filtered_records:
                sequences.append(str(record.seq).upper())

            matching_indices = NtCounter(sequences, threshold)
            range_indices = indices_to_ranges(matching_indices, min_length)

            for start_index, end_index in range_indices:
                block = f"block_{len(rbh2_entries) + 1}"
                rbh2_entry = {
                    "rbh": block,
                    "gene_group": "",
                    "color": ""
                }

                for record in filtered_records:
                    species = record.id.split('.')[0]
                    genomic_start, genomic_end = maf_coords(record, start_index, end_index)
                    scaffold = record.id.split('.', 1)[1] if '.' in record.id else record.id
                    strand = '-' if record.annotations['strand'] == -1 else '+'
                    rbh2_entry.update({
                        f"{species}_scaf": scaffold,
                        f"{species}_gene": block,
                        f"{species}_strand": strand,
                        f"{species}_start": genomic_start,
                        f"{species}_stop": genomic_end
                    })

                rbh2_entries.append(rbh2_entry)
        else:
            skip_count += 1

    columns = ["rbh", "gene_group", "color"]
    for species in species_list:
        columns.extend([
            f"{species}_scaf",
            f"{species}_gene",
            f"{species}_strand",
            f"{species}_start",
            f"{species}_stop"
        ])

    df = pd.DataFrame(rbh2_entries, columns=columns)
    output_path = f"{output_file}.rbh2" if not output_file.endswith('.rbh2') else output_file
    df.to_csv(output_path, sep="\t", index=False)

    print(f"MAF blocks: {maf_count}")
    print(f"Skipped blocks: {skip_count}")
    print(f"Conserved regions: {len(rbh2_entries)}")
    print(f"Output: {output_path}")

def main():
    args = get_args()
    maf_file = args.file
    Bp_Threshold    = args.threshold
    include_species = args.species
    output_filename = args.output

    process_alignments(maf_file, include_species, args.aligner, Bp_Threshold, output_filename, args.min_length)

if __name__ == '__main__':
    main()
