##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
import pandas as pd
from snakemake.utils import validate, logger

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################

def get_final_output():
    """
    Generate final output name
    """
    final_output = multiext(os.path.join(OUTPUT_FOLDER,'results','plots','gene_PA'), 
    	'.png', '.pdf')
    return final_output

##########################################################################

def infer_gene_constrains(seed_table):
    '''
    Infer gene_constrains from default config value or table
    '''

    list_constrains = []

    for index, row in seed_table.iterrows() :
        if 'evalue' in seed_table.columns:
            tmp_evalue = row.evalue
        else :
            tmp_evalue = config['default_blast_option']['e_val']

        if 'coverage' in seed_table.columns:
            tmp_coverage = row.coverage
        else :
            tmp_coverage = config['default_blast_option']['cov']

        if 'pident' in seed_table.columns:
            tmp_pident = row.pident
        else :
            tmp_pident = config['default_blast_option']['pident']


        tmp_text = f'{row.seed}_evalue_{tmp_evalue:.0e}_cov_{tmp_coverage}_pid_{tmp_pident}'

        list_constrains.append(tmp_text)

    return list_constrains

##########################################################################

def infer_ngs_groups(taxid_df):
    '''
    Infer taxid ngs option if not in taxid
    '''

    if 'NCBIGroups' not in taxid_df:
        taxid_df['NCBIGroups'] = config['ndg_option']['groups'] 


    return taxid_df

##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

# path to seeds sheet (TSV format, columns: seed, protein_id, ...)
seed_file = config['seed'] 

# Validation of the seed file
seed_table = (
    pd.read_table(seed_file)
    .set_index("seed", drop=False)
)

validate(seed_table, schema="../schemas/seeds.schema.yaml")

# path to taxonomic id to search seeds in (TSV format, columns: TaxId, NCBIGroups)
taxid = config['taxid']

# Validation of the taxid file
taxid_table = (
    pd.read_table(taxid)
    .set_index("TaxId", drop=False)
)

validate(taxid_table, schema="../schemas/taxid.schema.yaml")

##########################################################################
##########################################################################
##
##                        Core configuration
##
##########################################################################
##########################################################################

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_basedir_short__"] = os.path.basename(workflow.basedir)
config["__workflow_workdir__"] = os.getcwd()

if workflow.config_args :
    tmp_config_arg = '" '.join(workflow.config_args).replace('=', '="')
    config["__config_args__"] = f' -C {tmp_config_arg}"'
else :
    config["__config_args__"] = ''

with open(os.path.join(workflow.basedir, '../VERSION'), 'rt') as version:
    config["__workflow_version__"] = version.readline()
    config["__workflow_version_link__"] = f'https://github.com/vdclab/gene_presence_abscence/tree/main/version_2'


##########################################################################
##########################################################################
##
##                           Options
##
##########################################################################
##########################################################################

# Name your project
project_name = config['project_name']

# Result folder
OUTPUT_FOLDER =  os.path.join(config['output_folder'], project_name)
# Adding to config for report
config['__output_folder__'] = os.path.abspath(OUTPUT_FOLDER)

# Psiblast default e-value thershold
e_val_psiblast = config['psiblast_e_val'] 

# Option for ncbi_genome_download
section = config['ndg_option']['section'] 

# Values for assembly_levels :
assembly_levels = config['ndg_option']['assembly_levels'] 

# Values for refseq_categories : 
refseq_categories = config['ndg_option']['refseq_categories']

# Values for groups : 
taxid_table = infer_ngs_groups(taxid_table)

# Seepup option that create a reduce dataset using a psiblast step with the seed 
if config['speedup'] :
    speedup = os.path.join(OUTPUT_FOLDER, 'databases', 'reduce_taxid',
                    f'all_protein--eval_{e_val_psiblast:.0e}.fasta')
else  :
    speedup = os.path.join(OUTPUT_FOLDER, 'databases', 'all_taxid', 
                    'taxid_all_together.fasta')

# Definition of the requirements for each seed
gene_constrains = infer_gene_constrains(seed_table)

