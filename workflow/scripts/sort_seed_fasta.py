from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

with open(snakemake.output.fasta_for_psiblast, 'w') as general_out_file :
    with open(snakemake.input.fasta_seed, 'rt') as handle:
        for seed in SeqIO.parse(handle, 'fasta'):
            if seed.id in snakemake.params.seed2psiblast or seed.id[:-2] in snakemake.params.seed2psiblast:
                # Write a file per seed for psiblast
                SeqIO.write(seed, general_out_file, 'fasta')
