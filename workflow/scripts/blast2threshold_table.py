import pandas as pd
import sys

import utils_blast

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
    header = r_file.readline().rstrip().split('\t')
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

# Get the types of te columns for multiple HSPs dataframe
blast_dtypes = {'qseqid':'string',
                'sseqid':'string',
                'pident':'float64',
                'length':'int32',
                'mismatch':'int32',
                'gapopen':'int32',
                'qstart':'int32',
                'qend':'int32',
                'sstart':'int32',
                'send':'int32',
                'evalue':'float64',
                'bitscore':'float64'}

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

# seed_cov = snakemake.wildcards.coverage
# seed_pid = snakemake.wildcards.pid
# seed_eval = snakemake.wildcards.eval

# Read blast_out line by line
with open(snakemake.output[0], 'wt') as w_file:
    # Write the header in two times because format string need that
    header = '\t'.join(final_header)
    w_file.write(f"{header}\n")

    # Read the blast hsp by hsp
    for sub_blast in utils_blast.iterrator_on_blast_hsp(blast_out=blast_out) :
        # Get the number of hsps
        num_HSPs = len(sub_blast)

        if num_HSPs == 1:
            qseqid, sseqid, pident_blast, coverage_blast, evalue_blast, score = utils_blast.summarize_hit_only(split_line = sub_blast[0], 
                                                                                                               blast_header = blast_names,
                                                                                                               dict_protein = protein_dict,
                                                                                                               option_cov = snakemake.params.option_cov,
                                                                                                               option_pid = snakemake.params.option_pid)
        else:
            df_hsps = utils_blast.prepare_df_hsps(list_hsps = sub_blast,
                                                blast_dtypes = blast_dtypes, 
                                                blast_names = blast_names,
                                                HSPMIN = snakemake.params.minimum_length)

            sseqid = df_hsps.sseqid.values[0]
            qseqid = df_hsps.qseqid.values[0]

            delta_lg, coverage_blast, pident_blast, evalue_blast, score = utils_blast.summarize_hits(df_hsps = df_hsps, 
                                                                                                    length_query = protein_dict[qseqid], 
                                                                                                    length_subject = protein_dict[sseqid],
                                                                                                    option_cov = snakemake.params.option_cov, 
                                                                                                    option_pid = snakemake.params.option_pid)

        # Look if both proteins are in the family
        if qseqid in fam_protein_name and sseqid in fam_protein_name:
            # If exist put in the table because both are in the family
            line2write = f'{qseqid}\t{sseqid}\t{pident_blast}\t{evalue_blast}\t{coverage_blast}\tin_family_{seed_name}\n'
            w_file.write(line2write)

        # Look if one protein is in the family
        elif qseqid in fam_protein_name or sseqid in fam_protein_name:
            # If exist put in the table because one of them is in the family
            line2write = f'{qseqid}\t{sseqid}\t{pident_blast}\t{evalue_blast}\t{coverage_blast}\tout_family_{seed_name}\n'
            w_file.write(line2write)
                
##########################################################################
