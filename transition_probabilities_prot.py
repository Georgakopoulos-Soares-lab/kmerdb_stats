import argparse
import json
from Bio import SeqIO
import os
import gzip
import csv

def update_kmers_counts(kmers_counts, distance):
    """Update k-mers counts based on the distance from the last character not in the amino acids list seen."""
    for k in range(3, 8):
        if distance >= k:
            kmers_counts[f"kmers{k}_pos"] += 1

def count_aa_transitions_protein(fasta_file):
    """Count transitions between amino acids and k-mers distances in a protein FASTA file, ignoring characters not in the specified list."""
    amino_acids = "GALMFWKQESPVICYHRNDT"
    transitions = {aa1 + aa2: 0 for aa1 in amino_acids for aa2 in amino_acids}
    aa_counts = {aa: 0 for aa in amino_acids}
    kmers_counts = {f"kmers{k}_pos": 0 for k in range(3, 8)}
    last_seen_not_in_list = -1

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq)):
            if seq[i] not in amino_acids:
                last_seen_not_in_list = i
            else:
                if i - last_seen_not_in_list >= 6:
                    update_kmers_counts(kmers_counts, i - last_seen_not_in_list)
                if i < len(seq) - 1:
                    next_aa = seq[i+1]
                    if seq[i] in amino_acids and next_aa in amino_acids:
                        transitions[seq[i] + next_aa] += 1
                        aa_counts[seq[i]] += 1

    return transitions, aa_counts, kmers_counts

def calculate_aa_transition_probabilities(transitions, aa_counts):
    """Calculate the transition probabilities for each amino acid transition."""
    probabilities = {}
    for transition, count in transitions.items():
        prev_aa = transition[0]
        probabilities[transition] = count / aa_counts[prev_aa] if aa_counts[prev_aa] > 0 else 0
    return probabilities

def process_fasta_file(fasta_file_path, output_dir):
    file_name = os.path.basename(fasta_file_path).replace(".fasta", "")
    output_file_name = os.path.join(output_dir, f"{file_name}_output.csv")

    with open(fasta_file_path, 'rt') as fasta_file:
        transitions, aa_counts, kmer_counts = count_aa_transitions_protein(fasta_file)
        probabilities = calculate_aa_transition_probabilities(transitions, aa_counts)

    with open(output_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Type', 'Key', 'Value'])

        for transition, probability in sorted(probabilities.items()):
            writer.writerow(['Transition Probability', transition, f"{probability:.4f}"])

        for kmer, count in sorted(kmer_counts.items()):
            writer.writerow(['K-mer Count', kmer, count])

    print(f"Processed {file_name}. Results saved to {output_file_name}")

def main(bucket_id, output_dir):
    stats_file = "/scratch/mpp5977/proteomes_scheduler.json"  
    with open(stats_file) as f:
        stats = json.load(f)

    bucket_id_str = str(bucket_id)
    
    if bucket_id_str in stats:
        files_to_process = stats[bucket_id_str]
    else:
        print(f"No files found for bucket ID {bucket_id}.")
        return

    for gz_file_path in files_to_process:
        process_fasta_file(gz_file_path, output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process protein FASTA files based on bucket ID.")
    parser.add_argument("bucket_id", type=int, help="Bucket ID to process files for.")
    parser.add_argument("--output_dir", default="/path/to/output_dir/", help="Directory to store output files.")  
    args = parser.parse_args()

    main(args.bucket_id, args.output_dir)
