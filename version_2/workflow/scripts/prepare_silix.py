from Bio import SeqIO
import pandas as pd

### Amelioration purposes : Seems we could parse line per line the taxid file without loading it
# Preparing seeds
seed_table = pd.read_table(snakemake.input.seed_file)
seed_list = seed_table.protein_id.to_list()
seed_table.set_index('protein_id', inplace=True)

# Indexing fasta
fasta_db = SeqIO.index(snakemake.input.fasta, 'fasta')

# Opening blast_out and preparation
blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']

blast_out = pd.read_csv(snakemake.input.blast_out, sep='\t', header=None, names=blast_names)

# Don't know if necessery see later 
# blast_out['qseqid'] = blast_out.qseqid.apply(lambda x: x.split('.')[0])
# blast_out['sseqid'] = blast_out.sseqid.apply(lambda x: x.split('.')[0])

# I think it is totally useless
# blast_out.drop_duplicates(inplace=True)

# start filtering blast out on e-value, coverage and percent identity
blast_out = blast_out[blast_out.evalue <= float(snakemake.wildcards.eval)].reset_index(drop = True)
blast_out = blast_out[blast_out.pident >= float(snakemake.wildcards.pid)].reset_index(drop = True)

# Calculating the coverage on the query
for index, row in blast_out.iterrows() :
    blast_out.at[index, 'coverage'] = row.length / len(fasta_db[row.qseqid])

blast_out = blast_out[blast_out.coverage >= float(snakemake.wildcards.coverage)].reset_index(drop = True)

# Write the blast_out without the column cov because could messup in silix
blast_out[blast_names].to_csv(snakemake.output[0], sep='\t', index=False, header=False)