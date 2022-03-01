import pandas as pd
import sys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

blast_out = snakemake.input.blast_out

seed_family = snakemake.input.fnodes

# Seeds
seed_table = (
        pd.read_table(snakemake.input.seed_file,
                     index_col=0,
                     usecols=['protein_id', 'length'],
                     dtype={'protein_id':'string',
                             'length':'int16'},
                     )
)

# Protein length dict
protein_dict = {}

with open(snakemake.input.protein_file, 'rt') as r_file :
    # Remove the first line because header
    header = r_file.readline().split('\t')
    protein_index = header.index('protein_id')
    length_index = header.index('length')

    # Read the rest of the line
    for line in r_file:
        # Remove any blank character at the end
        r_line = line.rstrip()
        # Split by the tabulation as it is a tsv
        split_line = r_line.split('\t')
        
        # Put the information in the protein dict
        protein_dict[split_line[protein_index]] = int(split_line[length_index])

# Update with the seeds
protein_dict.update(seed_table.length
                              .to_dict())

# Opening blast_out and preparation
blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']

# Header of the output
final_header = ['protein1', 'protein2', 'pident', 'evalue', 'coverage', 'fam']

# Read the fnodes flushed file
seed_fam = pd.read_table(seed_family, usecols=['protein_id', 'seed'], dtype={'protein_id':'string', 'seed':'category'})

# Get all the protein id in the seed family
fam_protein_name = seed_fam.protein_id.tolist()

# Get the name of the seed
seed_name = snakemake.wildcards.seed

# Get all the threshold used by the user
split_info = seed_family.split('.fnodes')[0].split('_')

seed_cov = snakemake.wildcards.coverage
seed_pid = snakemake.wildcards.pid
seed_eval = snakemake.wildcards.eval

# Read blast_out line by line
with open(blast_out, 'rt') as r_file :
    with open(snakemake.output[0], 'wt') as w_file:
        # Write the header in two times because format string need that
        header = '\t'.join(final_header)
        w_file.write(f"{header}\n")

        # Read the blast line by line
        for line in r_file:            
            # Split the line to be easier to handle            
            line_split = line.split()

            # Get the information wanted
            evalue_blast = line_split[10]
            qseqid = line_split[0]
            sseqid = line_split[1]

            # Look if both proteins are in the family
            if qseqid in fam_protein_name and sseqid in fam_protein_name:
                # Try to save calculation time
                pident_blast = float(line_split[2]) / 100
                length = float(line_split[7]) - float(line_split[6]) + 1
                coverage_blast = length / protein_dict[qseqid]
                # If exist put in the table because both are in the family
                line2write = f'{qseqid}\t{sseqid}\t{pident_blast}\t{evalue_blast}\t{coverage_blast}\tin_family_{seed_name}\n'
                w_file.write(line2write)

            # Look if one protein is in the family
            elif qseqid in fam_protein_name or sseqid in fam_protein_name:
                # Try to save calculation time
                pident_blast = float(line_split[2]) / 100
                length = float(line_split[7]) - float(line_split[6]) + 1
                coverage_blast = length / protein_dict[qseqid]
                # If exist put in the table because one of them is in the family
                line2write = f'{qseqid}\t{sseqid}\t{pident_blast}\t{evalue_blast}\t{coverage_blast}\tout_family_{seed_name}\n'
                w_file.write(line2write)
                
##########################################################################
