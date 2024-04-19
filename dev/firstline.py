#!/usr/bin/env python

import argparse

def parse_args():
# Here I write some info about script 
    parser = argparse.ArgumentParser(
        description="Process input and output file paths, Output first line of input either in terminal or in output file.",

# Here I write what in- and outputs are. 
    parser.add_argument('--input', type=str, required=True, help='Path to the input file')  
    parser.add_argument('--output', type=str, help='Path to the output file:')  

    return parser.parse_args()

def process_files(input_path, output_path=None):
    with open(input_path, 'r') as file:
        first_line = file.readline().strip()
        message = f"Behold, the first line of the input file: {first_line}" # some random text, can and should be deleted later. 

    if output_path:
        with open(output_path, 'w') as file: # This part is for writting to output file IF provided. 
            file.write(message + "\n")
        print(f"Message written to {output_path}")  
    else:
        print(message) # This is when I didn't specify path to output.

def main():
    args = parse_args()
    process_files(args.input, args.output)  # Make sure these match the argument names defined in parse_args

if __name__ == "__main__":
    main()


# python firstline.py --input /Users/egortertyshnyk/Desktop/Simakov_Group/Conserved_regions/testFile.txt --output /Users/egortertyshnyk/Desktop/Simakov_Group/Conserved_regions/Output_of_testFile.txt
# Important, if you don't specify --output the first line will be printed just in terminal. 