# K-mer Extraction Tool - extractions_prot.py

This script extracts k-mers of a specified length from sequences in a zipped FASTA file. It filters k-mers to include only those composed of specified amino acids.

# Reference Proteom transition probability - transition_probabilities_prot.py

This script calculates the transition probabilities for an input reference proteome. For 20 different amino acids we have 380 transition probabilites.
The script also calculates the valid positions for kmer lengths of 3-7.


# merged_output.csv 

Is an example of the transition probabilities of all reference genomes in the kmerdb database.

# markov.py

Script that calculates the formation probability for all kmers within a kmer extraction file

## Installation

 Install the required Python packages:

```bash
pip install -r requirements.txt

