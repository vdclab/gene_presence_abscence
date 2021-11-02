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
protein_dict = {}

with open("snakemake.input.protein_file)", 'rt') as r_file :
    r_file.readline()
    for line in r_file:
        r_line = line.rstrip()
        split_line = r_line.split('\t')
        
        protein_dict[split_line[0]] = int(split_line[-1])

protein_dict.update(seed_table.length
                              .to_dict())

# Opening blast_out and preparation
blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']

# Calculating the max and min value of interest
max_eval = seed_table.evalue.max()
min_pident = seed_table.pident.min()
min_coverage = seed_table.coverage.min()

# Open all the output
list_open_output = [open(tsv_output, 'wt') for tsv_output in snakemake.output]
num_output = len(list_open_output)

# Read blast_out line by line
with open(snakemake.input.blast_out, 'rt') as r_file :
    for line in r_file:
        line_split = line.split()

        evalue_blast = float(line_split[10])
        pident_blast = float(line_split[2]) / 100
        length = float(line_split[7]) - float(line_split[6]) + 1
        qseqid = line_split[0]

        # start filtering blast out on e-value, coverage and percent identity
        if evalue_blast <= max_eval and pident_blast >= min_pident :
            coverage_blast = length / protein_dict[qseqid]

            # Calculating the coverage on the query
            if coverage_blast >= min_coverage :
                for index in range(num_output) :
                    file_out_path = snakemake.output[index]

                    constrains = file_out_path.split('--')[-1].split('.out')[0]

                    constrains_split = constrains.split('_')

                    evalue = float(constrains_split[2])
                    coverage = float(constrains_split[4])
                    pident = float(constrains_split[6]) / 100

                    if evalue_blast <= evalue and pident_blast >= pident and coverage_blast >= coverage:
                        list_open_output[index].write(line)

for file_open in list_open_output :
    file_open.close()
