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
import random

bed_intervals = [] # we need that for bed


def get_args():
    """
The purpose of this function is accepting parameter such as path to the file, threshold or included species as command-line arguments.

Template for input is:

python conservfinder.py -f "path_to_maf_file.maf" -t 0.7 -s "species1" "species2" "species3" -o "name_of_example_file.rbh2"

As for example:

python conservfinder.py -f "/Users/egortertyshnyk/Desktop/Simakov_Group/Conserved_regions/MOLLUSC_Chr10_small.maf" -t 0.6 -s "Octopusvulgaris6645.OX597823.1" "Octopusbimaculoides37653.NC_068990.1" "Octopussinensis2607531.NC_043005.1" -o "example.rbh2"

    """
    parser = argparse.ArgumentParser(description="Process MAF files to find conserved regions.")
    parser.add_argument('-f', '--file', required=True, help="Path to the MAF file")
    parser.add_argument('-t', '--threshold', type=float, default=0.6, help="Threshold for freq. of conserved nts")
    parser.add_argument('-s', '--species', nargs='+', help="List of species to include")
    parser.add_argument('-o', '--output', default='output.bed', help="Output filename for the BED file")
    return parser.parse_args()

def NtCounter(sequences, threshold):
    """
 The purpose of this function:
        The purporse of this function is to find the conserved nucleotides in a multiple sequence alignment.
        Conserved nucleotides -> nucleotides that occur in single column of alignment more often than in previously specified threshold. 

    The arguments are:
      - sequences: an iterable of strings, each string is a sequence
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
            if count / total >= threshold:
                conserved_indices.append(sn)
                break
    return conserved_indices

def indices_to_ranges(matching_indices):
    """
    Converts a list of indices to a list of start and end ranges.

    Arguments:
    - matching_indices: List of indices (int) where nucleotides are conserved.

    Returns:
    - List of tuples representing start and end indices of conserved regions.
    """
    range_indices = []
    if not matching_indices:  # Check if the list is empty
        return range_indices  # Return an empty list if there are no matching indices

    start = end = matching_indices[0]
    for index in matching_indices[1:] + [None]:
        if index is not None and index == end + 1:
            end = index
        else:
            if end - start >= 4:
                range_indices.append((start, end))
            if index is not None:
                start = end = index
    return range_indices

def ranges_to_coordinates(range_indices, sequences, records, Chrom_Position, ali_block_counter):
    """
    Converts ranges of conserved sequences into genomic coordinates and prints them.

    Arguments:
    - range_indices: List of tuples (start, end) of conserved regions.
    - sequences: List of sequences (str) analyzed.
    - records: List of record IDs (str) corresponding to sequences.
    - Chrom_Position: List of chromosome start positions (int) for each sequence.
    - ali_block_counter: counter (int) of .maf alighments.
    Note: for .bed format score of sequence is 0 by defauld and strand is unidentified (.).
    """
    for start, end in range_indices:
        conserved_sequence = sequences[0][start:end+1]
        for i, record_id in enumerate(records):
            genomic_start = Chrom_Position[i] + start + 1
            genomic_end = Chrom_Position[i] + end + 1
            bed_intervals.append([record_id, genomic_start, genomic_end, f"block_{ali_block_counter}", "0", "."])

def process_alignments(maf_file, species_list, threshold, output_bed):
    """
    Process the MAF file to find conserved regions and save them to an RBH2 file.

    Arguments:
    - maf_file: Path to the MAF file.
    - species_list: List of species to include in the analysis.
    - threshold: Threshold for the frequency of conserved nucleotides.
    - output_rbh2: Output filename for the RBH2 file.

    Packages used:
        - Bio.AlignIO
        - pandas
    """
    ali_block_counter = 1
    rbh2_entries = []
    aln_counter = 1
    total_conserved_regions = 0
    for multiple_alignment in AlignIO.parse(maf_file, "maf"):
        sequences, records, chrom_positions, strands, chroms = [], [], [], [], []

        for record in multiple_alignment:
            for species in species_list:
                if species in record.id:
                    sequences.append(str(record.seq).upper())
                    records.append(record.id)
                    chrom_positions.append(record.annotations['start'])

                    first_period_idx = record.id.find('.')
                    if first_period_idx != -1:
                        chrom_part = record.id[first_period_idx + 1:]  # Extract the scaffold part
                    else:
                        chrom_part = "unknown"  #  if no period is found
                    chroms.append(chrom_part)

                    strand = record.annotations['strand']
                    if strand == -1:
                        strand = '-'
                    elif strand == 1:
                        strand = '+'
                    strands.append(strand)
        aln_counter += 1
        if not sequences:
            continue

        matching_indices = NtCounter(sequences, threshold)
        range_indices = indices_to_ranges(matching_indices)

        for start_index, end_index in range_indices:
            r = lambda: random.randint(0,255) # This is here for random color generation
            rbh2_entry = {}
            rbh2_entry[f"rbh"] = f"block_{ali_block_counter}"
            rbh2_entry[f"gene_group"] = "not_assigned_yet"
            rbh2_entry[f"color"] = '#%02X%02X%02X' % (r(),r(),r())
            for i, record_id in enumerate(records):
                genomic_start = chrom_positions[i] + start_index
                genomic_end   = chrom_positions[i] + end_index
                rbh2_entry[f"{record_id}_scaf"]   = chroms[i]
                rbh2_entry[f"{record_id}_strand"] = strands[i]
                rbh2_entry[f"{record_id}_start"]  = genomic_start # make sure this is in bed format, note here whether it is or not
                rbh2_entry[f"{record_id}_stop"]   = genomic_end   # make sure this is in bed format, note here whether it is or not
            rbh2_entries.append(rbh2_entry)
            ali_block_counter += 1
            total_conserved_regions += len(range_indices)

    bed_output = output_bed if output_bed.endswith('.bed') else f"{output_bed}.bed"
    rbh2_output = bed_output.replace('.bed', '.rbh2')
    if rbh2_output.endswith('.rbh2'):
        rbh2_output = rbh2_output.replace('.rbh2', '')
    if not rbh2_output.endswith('.rbh2'):
        rbh2_output += '.rbh2'
    df = pd.DataFrame(rbh2_entries)
    df.to_csv(rbh2_output, sep="\t", index=False)

def main():
    args = get_args()
    maf_file = args.file
    Bp_Threshold    = args.threshold
    include_species = args.species
    output_filename = args.output

    process_alignments(maf_file, include_species, Bp_Threshold, output_filename)

if __name__ == '__main__':
    main()