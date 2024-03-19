import itertools
import sys

def generate_kmers(k, amino_acids):
    kmers = itertools.product(amino_acids, repeat=k)
    
    filename = f"{k}mers_proteins.txt"
    with open(filename, 'w') as file:
        for kmer in kmers:
            file.write(''.join(kmer) + '\n')
    
    print(f"All possible {k}-mers have been generated and stored in {filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <length_of_kmer>")
    else:
        k = int(sys.argv[1])  # Convert argument to integer
        amino_acids = "GALMFWKQESPVICYHRNDT"
        generate_kmers(k, amino_acids)

