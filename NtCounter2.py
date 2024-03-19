from collections import Counter

def NtCounter(sequences, threshold):

    for sn in range(len(sequences[0])):
        column = ''.join(seq[sn] for seq in sequences) 
        c = Counter(column)
        r = {}
        most_freq_seq = []
        print(column)
        for key in c:
            r[key] = c[key]/len(column)
        for nucleotide, frequency in r.items():
            if frequency >= threshold:
                most_freq_seq.append(nucleotide)
            else: 
                print("All below threshold!")
        print(r)
        print(most_freq_seq) 

# for test
threshold = 0.7
seq1 = "ATCATCATCATCGTT"
seq2 = "TCATCATCATCGTTT"
seq3 = "CATCATCATCGTTGA"
seq4 = "ATCATCATCGTTGAT"
sequences = [seq1,seq2,seq3,seq4]

# Call function 

NtCounter(sequences, threshold)