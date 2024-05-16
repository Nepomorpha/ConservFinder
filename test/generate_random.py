import random

output_file = "artificial_file.maf"

with open(output_file, 'w') as f:
    f.write("##maf version=1\n\n") # MAF header
    total_length = 0
    for i in range(10):
        num_sequences = 5  # Number of sequences
        num_diff1 = random.randint(5, 10)  # Number of different base pairs
        num_conserv1 = random.randint(10, 15)  # Number of conserved base pairs
        num_diff2 = random.randint(5, 20)  # Number of different base pairs
        num_same1 = random.randint(1, 5)  # Number of same bases that should be not couted by ConservFinder
        num_diff3 = random.randint(10, 15)  # Number of different base pairs
        num_same2 = random.randint(10, 15)   # Number of same base pairs
        bases = ['A', 'C', 'G', 'T']  # Possible bases

        sequences = {}

        # This is to create empty sequences
        for i in range(num_sequences):
            sequences[f"sequence{i+1}"] = ""


        for i in range(num_diff1):
            for seq_key in sequences.keys():
                sequences[seq_key] += random.choice(bases)


        conserved_base1 = random.choice(bases)
        for seq_key in sequences.keys():
            sequences[seq_key] += conserved_base1 * num_conserv1


        for i in range(num_diff2):
            for seq_key in sequences.keys():
                sequences[seq_key] += random.choice(bases)


        same_base1 = random.choice(bases)
        for seq_key in sequences.keys():
            sequences[seq_key] += same_base1 * num_same1


        for i in range(num_diff3):
            for seq_key in sequences.keys():
                sequences[seq_key] += random.choice(bases)


        same_base2 = random.choice(bases)
        for seq_key in sequences.keys():
            sequences[seq_key] += same_base2 * num_same2

        # Print the sequences to verify
        for key, sequence in sequences.items():
            print(f"{key}: {sequence}")

        print(f"First {num_diff1} are different. (Position 0 to {num_diff1})")
        print(f"Next {num_conserv1} are conserved. (Position {num_diff1} to {num_diff1+num_conserv1})")
        print(f"Next {num_diff2} are different. (Position {num_diff1+num_conserv1} to {num_diff1+num_conserv1+num_diff2})")
        print(f"Next {num_same1} are same. (Position {num_diff1+num_conserv1+num_diff2} to {num_diff1+num_conserv1+num_diff2+num_same1})")
        print(f"Next {num_diff3} are different. (Position {num_diff1+num_conserv1+num_diff2+num_same1} to {num_diff1+num_conserv1+num_diff2+num_same1+num_diff3})")
        print(f"Next {num_same2} are same. (Position {num_diff1+num_conserv1+num_diff2+num_same1+num_diff3} to {num_diff1+num_conserv1+num_diff2+num_same1+num_diff3+num_same2})")
        print(f"Conserved sequences are from position: \n {num_diff1} to {num_diff1+num_conserv1} \n {num_diff1+num_conserv1+num_diff2} to {num_diff1+num_conserv1+num_diff2+num_same1} \n {num_diff1+num_conserv1+num_diff2+num_same1+num_diff3} to {num_diff1+num_conserv1+num_diff2+num_same1+num_diff3+num_same2}")

        full_seq_or_something_like_that = random.randint(1000, 5000)


        f.write("a\n")
        for key, sequence in sequences.items():
            f.write(f"s\t{key}\t{total_length}\t{len(sequence)}\t+\t {full_seq_or_something_like_that} \t{sequence}\n")
        f.write("\n")

        total_length += len(next(iter(sequences.values())))