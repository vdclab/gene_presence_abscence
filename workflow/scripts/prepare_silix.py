import pandas as pd
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Seeds
seed_dict = (
        pd.read_table(snakemake.input.seed_file)
        .set_index('protein_id')
        .length
        .to_dict()
)

# Protein lenght dict
protein_dict = (
        pd.read_table(snakemake.input.protein_file)
        .set_index('protein_id')
        .length
        .to_dict()
)

protein_dict.update(seed_dict)

# Opening blast_out and preparation
blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']

blast_out = pd.read_csv(snakemake.input.blast_out, sep='\t', header=None, names=blast_names)

# I think it is totally useless
# blast_out.drop_duplicates(inplace=True)

# start filtering blast out on e-value, coverage and percent identity
blast_out = blast_out[blast_out.evalue <= float(snakemake.wildcards.eval)].reset_index(drop = True)
blast_out = blast_out[blast_out.pident >= float(snakemake.wildcards.pid)].reset_index(drop = True)

# Calculating the coverage on the query
for index, row in blast_out.iterrows() :
    blast_out.at[index, 'coverage'] = row.length / protein_dict[row.qseqid]

blast_out = blast_out[blast_out.coverage >= float(snakemake.wildcards.coverage)].reset_index(drop = True)

# Write the blast_out without the column cov because could messup in silix
blast_out[blast_names].to_csv(snakemake.output[0], sep='\t', index=False, header=False)