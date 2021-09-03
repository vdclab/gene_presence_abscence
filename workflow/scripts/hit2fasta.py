from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Get all the proteins of the database inside a dict-like structure
all_index_fasta = SeqIO.index(snakemake.input.protein_fasta, 'fasta')

# Read line by line the protein table hit without loading it in memory
with open(snakemake.input.list_all_prot, 'rt') as r_file :
    with open(snakemake.output.fasta, 'wt') as w_file :
        header = r_file.readline().rstrip().split() 
        index_proteinId = header.index('protein_id') 

        for line in r_file : 
            tmp_line = line.split() 
            proteins_of_interest = tmp_line[index_proteinId]
            SeqIO.write(all_index_fasta[proteins_of_interest], w_file, 'fasta')
        