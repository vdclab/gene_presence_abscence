import pandas as pd
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Seeds
seed_table = (
        pd.read_table(snakemake.input.seed_file)
        .set_index('protein_id')
)

# Protein lenght dict
protein_dict = (
        pd.read_table(snakemake.input.protein_file)
        .set_index('protein_id')
        .length
        .to_dict()
)

protein_dict.update(seed_table.length
                              .to_dict())

# Opening blast_out and preparation
blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']

blast_out = pd.read_csv(snakemake.input.blast_out, sep='\t', header=None, names=blast_names)

# I think it is totally useless
# blast_out.drop_duplicates(inplace=True)
max_eval = seed_table.evalue.max()
min_pident = seed_table.pident.min()
min_coverage = seed_table.coverage.min()

# start filtering blast out on e-value, coverage and percent identity
blast_out = blast_out[blast_out.evalue <= max_eval].reset_index(drop = True)
blast_out = blast_out[blast_out.pident >= min_pident].reset_index(drop = True)

# Calculating the coverage on the query
for index, row in blast_out.iterrows() :
    blast_out.at[index, 'coverage'] = row.length / protein_dict[row.qseqid]

blast_out = blast_out[blast_out.coverage >= min_coverage].reset_index(drop = True)

# Because blast_out could be really big, generate all the output at once
for file_out in snakemake.output :
    tmp_blast_out = blast_out.copy()
    constrains = file_out.split('--')[-1].split('.')[0]
    constrains_split = constrains.split('_')

    seed = constrains_split[0]
    evalue = float(constrains_split[2])
    coverage = float(constrains_split[4])
    pident = float(constrains_split[6])

    tmp_blast_out = tmp_blast_out[tmp_blast_out.evalue <= evalue].reset_index(drop = True)
    tmp_blast_out = tmp_blast_out[tmp_blast_out.pident >= pident].reset_index(drop = True)
    tmp_blast_out = tmp_blast_out[tmp_blast_out.coverage >= coverage].reset_index(drop = True)

    # Write the blast_out without the column cov because could messup in silix
    tmp_blast_out[blast_names].to_csv(file_out, sep='\t', index=False, header=False)