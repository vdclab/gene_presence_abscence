import pandas as pd
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Opening blastout
blast_names = ['qacc', 'qlen', 'qseq','qstart', 'qend', 'sacc', 'slen', 'sseq', 'sstart', 'send','length',
               'pident', 'evalue', 'bitscore', 'qcovs','qcovhsp', 'ssciname', 'sblastname', 'stitle']

psiblast_result = pd.read_table(snakemake.input.psiblast,
                              comment='#',
                              names=blast_names
                              )

# Cleaning blastout
psiblast_result = psiblast_result[psiblast_result.qacc != 'Search has CONVERGED!']
psiblast_result.to_csv(snakemake.output.clean_blast, sep='\t', index=False)

# Getting the list of protein matches
all_sacc = psiblast_result.sacc.unique()

with open(snakemake.output.list_all_prot, 'w') as w_file :
    w_file.write('protein_id\n')

    for sacc in all_sacc :
        w_file.write(f'{sacc}\n')