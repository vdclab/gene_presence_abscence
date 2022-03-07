from posixpath import split
import pandas as pd
import sys
from common import utils_blast

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Seeds
seed_table = (
        pd.read_table(snakemake.input.seed_file)
        .set_index('protein_id')
)

# Protein lenght dict
# protein_dict = {}

# with open(snakemake.input.protein_file, 'rt') as r_file :
#     r_file.readline()
#     for line in r_file:
#         r_line = line.rstrip()
#         split_line = r_line.split('\t')
        
#         protein_dict[split_line[0]] = int(split_line[-1])

# protein_dict.update(seed_table.length
#                               .to_dict())

# Blast_out and preparation
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

# Calculating the max and min value of interest
max_eval = seed_table.evalue.max()
# min_pident = seed_table.pident.min()
# min_coverage = seed_table.coverage.min()

# Open all the output
list_open_output = [open(tsv_output, 'wt') for tsv_output in snakemake.output]
num_output = len(list_open_output)

# Read the blast hsp by hsp
for sub_blast in utils_blast.iterrator_on_blast_hsp(blast_out=snakemake.input.blast_out) :
    # Variable to get the lines 
    line = ""

    # Get the number of hsps
    num_HSPs = len(sub_blast)

    if num_HSPs == 1:
        evalue_blast = float(sub_blast[0][blast_names.index('evalue')])
        line = "\t".join(sub_blast[0]) + "\n"
    else:
        df_hsps = utils_blast.prepare_df_hsps(list_hsps = sub_blast,
                                            blast_dtypes = blast_dtypes, 
                                            blast_names = blast_names,
                                            HSPMIN = snakemake.params.minimum_length)

        evalue_blast = df_hsps.evalue.max()

        for index in df_hsps.index:
            line += "\t".join(sub_blast[index]) + "\n"                

    if evalue_blast <= max_eval :
    #if evalue_blast <= max_eval and pident_blast >= min_pident :
        # coverage_blast = length / protein_dict[qseqid]

        # Calculating the coverage on the query
        # if coverage_blast >= min_coverage :
        for index in range(num_output) :
            file_out_path = snakemake.output[index]

            constrains = file_out_path.split('--')[-1].split('.out')[0]

            constrains_split = constrains.split('_')

            # As i have the control about the end of the file name it is better to go this way in case "_" in the gene name
            evalue = float(constrains_split[-5])
            # coverage = float(constrains_split[-3])
            # pident = float(constrains_split[-1])

            # if evalue_blast <= evalue and pident_blast >= pident and coverage_blast >= coverage:
            if evalue_blast <= evalue:
                list_open_output[index].write(line)

for file_open in list_open_output :
    file_open.close()
