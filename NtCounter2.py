from collections import Counter

threshold = 0.7

seq1 = "ATCATCATCATCGTT"
seq2 = "TCATCATCATCGTTT"
seq3 = "CATCATCATCGTTGA"
seq4 = "ATCATCATCGTTGAT"

sequences = [seq1,seq2,seq3,seq4]

for sn in range(len(sequences[0])):
    column = seq1[sn] + seq2[sn] + seq3[sn] + seq4[sn]
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
    #print(2/len(column)) 

