#!/usr/bin/env python
"""
TODO:
  - Write a docstring for the whole script here: Use something like the template
  - Delete the programs that you aren't using
    - If you want to have some code that you are testing or are working on, make a folder in the repo called "dev" and put files in there.
  - Rename this program to something that makes sense, like conservfinder.py, just use lowercase characters

Do something like this - Delete the irrelevant parts and replace with your own info

Program  : odp
Language : snakemake
Date     : 2021-01-21
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
Support  : For issues or questions, please search if the topic has been discussed already
           on github and open a new issue if not: https://github.com/conchoecia/odp/issues
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.
Citation : If you use this software for your scientific publication, please cite:
           Schultz, DT; Haddock, SHD; Bredeson, JV; Green, RE; Simakov, O & Rokhsar, DS
           Ancient gene linkages support ctenophores as sister to other animals. Nature (2023).
           https://doi.org/10.1038/s41586-023-05936-6

Description:
  This program is part of the Oxford Dot Plot (odp) package on github.
  This script performs mutual-best protein diamond BLAST searches,
  then makes synteny plots of those results. Performs many other functions.
  See the documentation for more details.

Usage instructions:
  - See https://github.com/conchoecia/odp#getting-started
"""
from collections import Counter
from Bio import AlignIO

def NtCounter(sequences, threshold):
    """
    This is a docstring. The point of a docstring is to document what this piece of code does.
    It is for the person writing the code, for their own reference, and for people who use the code
    to do something.

    The purpose of this function:
        The purporse of this function is to find the conserved nucleotides in a multiple sequence alignment.

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
# delete trailing spaces
Bp_Threshold = 0.5  # for __% 

# Put this into a main function()
for multiple_alignment in AlignIO.parse(maf_file, "maf"):
    # delete trailing spaces
    print("\n--------------------------------New Alignment Block--------------------------------")      
    sequences = []
    records = []
    # delete trailing spaces
    Chrom_Position = [] 
    
    for record in multiple_alignment:
        if record.id in include_species:
            sequences.append(str(record.seq).upper())
            # delete trailing spaces
            records.append(record.id)  
            Chrom_Position.append(record.annotations['start'])

    if not sequences:
        continue

    # delete whitespaces
   
    matching_indices = NtCounter(sequences, Bp_Threshold)
    print(matching_indices)

    # delete whitespaces
   
    range_indices = []
    # delete trailng spaces
    if matching_indices:  
        # put this into a function called indices_to_ranges()
        # For each function write a docstring
        start = end = matching_indices[0]
        # delete trailing spaces
        for index in matching_indices[1:] + [None]: 
            if index is not None and index == end + 1:
                end = index
            else:
                # delete trailing spaces
                if end - start >= 4: 
                    range_indices.append((start, end))
                if index is not None:
                    start = end = index

        # put this into a function called ranges_to_coordinates()
        # For each function write a docstring
        # delete trailing spaces
        for start, end in range_indices:  
            conserved_sequence = sequences[0][start:end+1]
            print(f"Conserved sequence: {conserved_sequence} Relative range: {start}-{end}")
            for i, record_id in enumerate(records):
                # delete trailing spaces
                genomic_start = Chrom_Position[i] + start + 1 
                genomic_end = Chrom_Position[i] + end + 1
                print(f"{record_id} {genomic_start} {genomic_end}")
            # delete trailing spaces
            print()  

# put the piece of code down here that calls the main function
# if __name__ == '__main__': ... et cetera