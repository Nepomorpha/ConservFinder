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

TODO:
  - Write a docstring for the whole script here: Use something like the template (DONE)
  - Delete the programs that you aren't using
    - If you want to have some code that you are testing or are working on, make a folder in the repo called "dev" and put files in there.
  - Rename this program to something that makes sense, like conservfinder.py, just use lowercase characters

"""
from collections import Counter
from Bio import AlignIO

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

maf_file = '/Users/egortertyshnyk/Desktop/Simakov_Group/Conserved_regions/MOLLUSC_Chr10_small.maf'
include_species = ["Octopusvulgaris6645.OX597823.1", "Octopusbimaculoides37653.NC_068990.1", "Octopussinensis2607531.NC_043005.1"]

Bp_Threshold = 0.5  # for __%

# Put this into a main function()
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
    print(matching_indices)

    range_indices = []
    if matching_indices:
        # put this into a function called indices_to_ranges()
        # For each function write a docstring
        start = end = matching_indices[0]
        for index in matching_indices[1:] + [None]:
            if index is not None and index == end + 1:
                end = index
            else:
                if end - start >= 4:
                    range_indices.append((start, end))
                if index is not None:
                    start = end = index

        # put this into a function called ranges_to_coordinates()
        # For each function write a docstring
        for start, end in range_indices:
            conserved_sequence = sequences[0][start:end+1]
            print(f"Conserved sequence: {conserved_sequence} Relative range: {start}-{end}")
            for i, record_id in enumerate(records):
                genomic_start = Chrom_Position[i] + start + 1
                genomic_end = Chrom_Position[i] + end + 1
                print(f"{record_id} {genomic_start} {genomic_end}")
            print()

# put the piece of code down here that calls the main function
# if __name__ == '__main__': ... et cetera
