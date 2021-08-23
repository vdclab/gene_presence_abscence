##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
from textwrap import dedent
import glob

import pandas as pd
from snakemake.utils import validate

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################

def get_final_output():
    final_output = multiext(os.path.join(OUTPUT_FOLDER,'results','plots','gene_PA'), 
    	'.png', '.pdf')
    return final_output

##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

# path to taxonomic id to search seeds in (TSV format, columns: TaxId, NCBIGroups)
taxid = config['taxid']

# path to seeds sheet (TSV format, columns: seed, protein_id, ...)
seed_file = config['seed'] 

# Name your project
project_name = config['project_name']

# Result folder
OUTPUT_FOLDER =  os.path.join(config['output_folder'], project_name)

# Blast e-value thershold
e_val = config['e_val'] 

# Option for ncbi_genome_download
section = config['ndg_option']['section'] 

# Values for assembly_levels :
assembly_levels = config['ndg_option']['assembly_levels'] 

# Values for refseq_categories : 
refseq_categories = config['ndg_option']['refseq_categories']

# Values for groups : 
groups = config['ndg_option']['groups'] 

# Seepup option that create a reduce dataset using a psiblast step with the seed 
if config['speedup'] :
    speedup = os.path.join(OUTPUT_FOLDER, 'results', 
    				f'all_protein--eval_{e_val:.0e}.fasta')
else  :
    speedup = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 
    				'taxid_all_together.fasta')

##########################################################################
##########################################################################
##
##                                Main
##
##########################################################################
##########################################################################

# Validation of the seed file
seed_table = (
    pd.read_table(seed_file)
    .set_index("seed", drop=False)
)

validate(seed_table, schema="../schemas/seeds.schema.yaml")

# Validation of the taxid file
taxid_table = (
    pd.read_table(taxid)
    .set_index("TaxId", drop=False)
)

validate(taxid_table, schema="../schemas/taxid.schema.yaml")



##########################################################################
##########################################################################
##
##                                Config
##
##########################################################################
##########################################################################

# Name your project, take the name of seed by default
project_name = config['project_name']

# Blast e-value thershold, 0.000001 by default but can be changed in -C
e_val = config['e_val'] 

# Option for ncbi_genome_download
# Values for section : {refseq,genbank}
section = config['ndg_option']['section'] 

# Values for assembly_levels : ['all', 'complete', 'chromosome', 'scaffold', 'contig']
assembly_levels = config['ndg_option']['assembly_levels'] 

# Values for refseq_categories : {'reference', 'all'}
refseq_categories = config['ndg_option']['refseq_categories']

# Values for groups : ['all', 'archaea', 'bacteria', 'fungi', 'invertebrate', 'metagenomes', 'plant', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral']
groups = config['ndg_option']['groups'] 

OUTPUT_FOLDER =  os.path.join(config['output_folder'], project_name)

# Definition of the requirements for each seed
gene_constrains = [f'{row.seed}_evalue_{row.evalue:.0e}_cov_{row.coverage}_pid_{row.pident}'
                                                    for index, row in seed_table.iterrows()]

