from Bio import SeqIO
import sys, os

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

parser_fasta = SeqIO.parse(snakemake.input.fasta_protein, 'fasta')

with open(snakemake.output.tsv, 'wt') as w_file:
    header = ['protein_id', 'length', 'protein_description']
    w_file.write('{"\t".join(header)}\n')

    for prot in parser_fasta:
        length = len(prot.seq)

        w_file.write("{prot.id}\t{length}\t{prot.description}\n")

        