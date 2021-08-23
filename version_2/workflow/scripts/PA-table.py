import pandas as pd

# Seed preparing
seed_table = pd.read_table(str(input.seed_file))
seed_list = seed_table.seed.to_list()
seed_color_dict = seed_table.set_index('seed').color.to_dict()

# list of all proteins
all_proteins = pd.read_table(str(input.protein_table))

# fnodes opening
fnodes_files = [pd.read_table(fnodes_file) for fnodes_file in input.fnodes]

# concat all fnodes dataframe to one and add the protein information from the protein table and genome table
fam_id_table = pd.concat(fnodes_files)
fam_id_table = fam_id_table.merge(all_proteins, on='protein_id')

patab = pd.crosstab(index = fam_id_table['genome_id'], columns = fam_id_table['seed'])

# To add missing seed if not find
seed_missing = [seed for seed in seed_list if seed not in patab.columns]
patab.loc[:,seed_missing] = 0

patab = patab[seed_list].sort_values(by = seed_list, ascending = False).reset_index()

# Add the genome name to the table in case needed
patab = patab.merge(fam_id_table[['genome_id','genome_name']].drop_duplicates(), on='genome_id')

patab = patab.melt(id_vars=['genome_id', 'genome_name'], var_name='seed', value_name='PA')

# Put color instead of number
for index, row in patab.iterrows():
    # Use the fact that 0 == False in python to test if it's 1 or 0
    if row.PA :
        patab.at[index, 'color'] = seed_color_dict[row.seed]
    else :
        patab.at[index, 'color'] = '#FFFFFF' # White color

# save the table
patab.to_csv(str(output.final_table), sep='\t', index=False)
