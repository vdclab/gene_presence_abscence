from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Get all the proteins of the database inside a dict-like structure
all_index_fasta = SeqIO.index(snakemake.input.protein_fasta, 'fasta')

# Read line by line the protein ftable hit without loading it in memory
with open(snakemake.input.list_all_prot, 'rt') as r_file :
    header = r_file.readline().rstrip().split() 
    index_proteinId = header.index('protein_id') 
    proteins_of_interest = [] 

    for line in r_file : 
        tmp_line = line.split() 
        proteins_of_interest.append(tmp_line[index_proteinId]) 

# Filtering protein table and saving
with open(snakemake.output.fasta, 'w') as w_file :
    for protein in proteins_of_interest :
        SeqIO.write(all_index_fasta[protein], w_file, 'fasta')