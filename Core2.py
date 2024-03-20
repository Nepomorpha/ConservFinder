from collections import Counter
from Bio import AlignIO

def NtCounter(sequences, threshold):
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
        start = end = matching_indices[0]
        for index in matching_indices[1:] + [None]: 
            if index is not None and index == end + 1:
                end = index
            else:
                if end - start >= 4: 
                    range_indices.append((start, end))
                if index is not None:
                    start = end = index

        for start, end in range_indices:  
            conserved_sequence = sequences[0][start:end+1]
            print(f"Conserved sequence: {conserved_sequence} Relative range: {start}-{end}")
            for i, record_id in enumerate(records):
                genomic_start = Chrom_Position[i] + start + 1 
                genomic_end = Chrom_Position[i] + end + 1
                print(f"{record_id} {genomic_start} {genomic_end}")
            print()  
