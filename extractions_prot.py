from Bio import SeqIO
import gzip
import sys
import os

def extract_kmers(sequence, k, valid_chars):
    """
    Extracts and returns k-mers from a given sequence that only contain valid characters.
    """
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmer = str(sequence[i:i+k])
        if all(char in valid_chars for char in kmer):
            kmers.add(kmer)
    return kmers

def process_fasta(file_path, k, output_file):
    """
    Processes a gzipped FASTA file to extract k-mers that only contain specified amino acids.
    """
    AminoAcidsL = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    kmers_set = set()
    
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            kmers_set.update(extract_kmers(record.seq, k, AminoAcidsL))
    
    with open(output_file, 'w') as out_file:
        for kmer in sorted(kmers_set):
            out_file.write(f"{kmer}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_kmers_gz.py <length of k-mers> <path to gzipped FASTA file>")
        sys.exit(1)

    k = int(sys.argv[1])
    input_file = sys.argv[2]
    input_basename = os.path.basename(input_file)
    # Modify here to ensure the output file name format is as desired
    file_root = os.path.splitext(input_basename)[0]  # Remove .gz extension
    file_root = os.path.splitext(file_root)[0]  # Remove .fasta extension if present
    output_file = f"{file_root}_{k}mers.txt"
    
    process_fasta(input_file, k, output_file)
    print(f"Process completed. {k}-mers written to {output_file}")

