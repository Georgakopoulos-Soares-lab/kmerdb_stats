import sys

def find_nullomers(input_file):
    k = len(open(input_file, 'r').readline().strip())
    all_kmers_file = f"{k}mers_proteins.txt"
    
    existing_kmers = set()
    with open(input_file, 'r') as file:
        for line in file:
            existing_kmers.add(line.strip())

    try:
        with open(all_kmers_file, 'r') as file:
            all_kmers = {line.strip() for line in file}
    except FileNotFoundError:
        print(f"File {all_kmers_file} not found. Please make sure it exists in the current directory.")
        sys.exit(1)

    nullomers = all_kmers - existing_kmers
    
    
    output_file = f"{k}nullomers_protein.txt"
    with open(output_file, 'w') as outfile:
        for kmer in sorted(nullomers):
            outfile.write(kmer + '\n')
    
    print(f"Nullomers for length {k} have been generated and stored in {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("input_file for kmer extraction needed")
    else:
        input_file = sys.argv[1]  
        find_nullomers(input_file)

