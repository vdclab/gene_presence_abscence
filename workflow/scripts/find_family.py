import pandas as pd
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Open fnodes
fnodes = pd.read_table(snakemake.input.fnodes, names=["family", "protein_id"])

seed = snakemake.wildcards.seed

### Amelioration purposes : Seems we could parse line per line the taxid file without loading it
seed_df = pd.read_table(snakemake.input.seed_file)
seed_protname = seed_df.set_index("seed").protein_id.to_dict()[seed]

# Detection families WARNING verification about the name of the seed and the name of the protein in the fasta at the end may be different
gene_family = fnodes.loc[
    fnodes.protein_id.str.contains(seed_protname), "family"
].unique()[0]

# writing file with only family
fnodes = fnodes.loc[fnodes.family == gene_family]

# Creation of a columns 'seed' to identify quickly the seed and familiy
fnodes["seed"] = seed

fnodes.to_csv(snakemake.output[0], sep="\t", index=False)
