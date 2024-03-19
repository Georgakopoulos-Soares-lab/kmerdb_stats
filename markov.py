import gzip
import pandas as pd
import json
import sys
import os

def get_identifier(filename):
    parts = filename.split('_')
    return '_'.join(parts[:2])

def calculate_probability(kmer, probabilities):
    prob = 1.0
    k = len(kmer)
    for i in range(len(kmer) - 1):
        pair = kmer[i:i+2]
        prob *= probabilities.get(pair, 0)  # Default to 0 
    return prob

def process_files(csv_filepath, gz_filepath, output_directory):
    csv_data = pd.read_csv(csv_filepath)
    transition_probabilities = {}

    for index, row in csv_data.iterrows():
        identifier = get_identifier(row['FileName'])
        transition_probabilities[identifier] = row[1:24].to_dict()

    gz_identifier = get_identifier(os.path.basename(gz_filepath))

    if gz_identifier not in transition_probabilities:
        print(f"Identifier {gz_identifier} not found in the CSV file.")
        return

    probabilities = transition_probabilities[gz_identifier]

    os.makedirs(output_directory, exist_ok=True)
    output_filename = os.path.join(output_directory, os.path.basename(gz_filepath).replace('.gz', '_probabilities.gz'))

    with gzip.open(gz_filepath, 'rt') as gz_file, gzip.open(output_filename, 'wt', compresslevel=9) as out_file:
        for kmer in gz_file:
            kmer = kmer.strip()  # Remove newline and any surrounding whitespace
            prob = calculate_probability(kmer, probabilities)
            # Format the probability with 4 decimal places in scientific notation
            formatted_prob = f"{prob:.4e}"
            out_file.write(f"{kmer} {formatted_prob}\n")

    print(f"Processed {gz_filepath}. Output written to {output_filename}.")

def main(bucket_id, output_dir):
    with open('/scratch/mpp5977/viruses_kmers.json', 'r') as json_file:
        file_paths_dict = json.load(json_file)

    if bucket_id not in file_paths_dict:
        print(f"Bucket ID {bucket_id} not found in JSON file.")
        return

    file_paths = file_paths_dict[bucket_id]

    csv_file = '/storage/group/izg5139/default/kmerdb_clean/merged_output.csv'  
    for gz_file in file_paths:
        process_files(csv_file, gz_file, output_dir)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python markov.py <bucket_id> <output_dir>")
        sys.exit(1)
    
    bucket_id = sys.argv[1]
    output_dir = sys.argv[2]
    main(bucket_id, output_dir)

