import pandas as pd
import sys, os

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

prot_perso = snakemake.input.perso
perso_df = pd.read_table(prot_perso)

prot_taxid = snakemake.input.taxid
taxid_df = pd.read_table(prot_taxid)

merge_df = pd.concat([taxid_df, perso_df])

merge_df.to_csv(snakemake.output.tsv, index=False, sep='\t')