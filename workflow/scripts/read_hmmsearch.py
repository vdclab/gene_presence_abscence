import sys, os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Columns for hmm out
column_names = ['protein_id',
                'protein_accession',
                'protein_length',
                'query',
                'query_accession',
                'query_length',
                'full_evalue',
                'full_score',
                'full_bias',
                'domain_#',
                'domain_of',
                'domain_c-evalue',
                'domain_c-evalue',
                'domain_score',
                'domain_bias',
                'query_from',
                'query_to',
                'ali_from',
                'ali_to',
                'env_from',
                'env_to',
                'target_description']

all_ids = []

with open(snakemake.input.hmm_out, 'r') as r_file:
    with open(snakemake.output.list_all_prot, 'w') as w_file:
        w_file.write('protein_id\n')
        for line in r_file:
            split_line = line.split()
            prot_id = split_line[0]
            if not prot_id == '#' and not prot_id in all_ids:
                all_ids.append(prot_id)
                w_file.write(f'{prot_id}\n')
