#!/usr/bin/env python
"""

Program     : ConservFinder for .maf files.
Language    : Python
Author      : Yehor Tertyshnyk, Darrin T. Schultz
Email       : egor.tertyshnyk@gmail.com
Github      : https://github.com/Nepomorpha/ConservFinder
License     : ?--?
Citation    : -TBS-

Description : This program takes .maf file and output conserved sequences as well as their genomic range; relative to beginnning of alignment (1) and relative to position on chromosome (2). 

Usage
instructions: see -TBS-

"""
from collections import Counter
from Bio import AlignIO
import pybedtools

bed_intervals = [] # we need that for bed

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
    """
    for start, end in range_indices:
        conserved_sequence = sequences[0][start:end+1]
        print(f"Conserved sequence: {conserved_sequence} Relative range: {start}-{end}")
        for i, record_id in enumerate(records):
            genomic_start = Chrom_Position[i] + start + 1
            genomic_end = Chrom_Position[i] + end + 1
            print(f"{record_id} {genomic_start} {genomic_end}")
            bed_intervals.append([record_id, genomic_start, genomic_end, "0", ".", f"block_{ali_block_counter}"])
        print()

def main():
    maf_file = '/Users/egortertyshnyk/Desktop/Simakov_Group/Conserved_regions/MOLLUSC_Chr10_small.maf'
    include_species = ["Octopusvulgaris6645.OX597823.1", "Octopusbimaculoides37653.NC_068990.1", "Octopussinensis2607531.NC_043005.1"]
    Bp_Threshold = 0.6
    ali_block_counter = 1 # that is needed to track which maf ali is shown in bed

    for multiple_alignment in AlignIO.parse(maf_file, "maf"):
        print("\n--------------------------------New Alignment Block--------------------------------")
        sequences = []
        records = []
        Chrom_Position = []

        for record in multiple_alignment:
            if record.id in include_species:
                sequences.append(str(record.seq).upper())
                records.append(record.id)
                Chrom_Position.append(record.annotations['start'])

        if not sequences:
            continue

        matching_indices = NtCounter(sequences, Bp_Threshold)
        range_indices = indices_to_ranges(matching_indices)
        ranges_to_coordinates(range_indices, sequences, records, Chrom_Position, ali_block_counter)
        ali_block_counter += 1
    if bed_intervals:
        bedtool = pybedtools.BedTool(bed_intervals)
        bedtool.saveas('example.bed')
    else:
        print("No intervals to save.")

if __name__ == '__main__':
    main()