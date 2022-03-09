import pandas as pd
from ete3 import NCBITaxa 
import os, sys
import gzip
from Bio import SeqIO
from urllib.request import urlopen
import subprocess
import shlex
import ncbi_genome_download as ngd 
import resource
import numpy as np

##########################################################################

def get_cmdline_ndg(section, flat_output, file_formats, assembly_levels, 
                    output, metadata_table, taxids, groups, refseq_categories, parallel):
    '''
    Launch the software in cmdline instead of function
    '''
    
    python_version = sys.version.split('.')[0]
    ncbi_genome_download = ngd.__file__.replace('init', 'main')

    if flat_output :
        cmd_line = f"python{python_version} {ncbi_genome_download} -s {section} -F {file_formats}\
                                      -l {assembly_levels} --flat-output\
                                      -o {output} -p {parallel}\
                                      -m {metadata_table} -R {refseq_categories}\
                                      -t {','.join(taxids)} {groups}"
    else :
        cmd_line = f"python{python_version} {ncbi_genome_download} -s {section} -F {file_formats}\
                                      -l {assembly_levels}\
                                      -o {output} -p {parallel}\
                                      -m {metadata_table} -R {refseq_categories}\
                                      -t {','.join(taxids)} {groups}"

    subprocess.run(shlex.split(cmd_line))
    return

##########################################################################


def find_NCBIgroups(taxid, ncbi):
    """
    Infer taxid groups
    """

    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    names = [names[taxid] for taxid in lineage]

    protazoa = {'Cryptophyceae', 'Apusozoa', 'Amoebozoa', 'Haptista', 'Metamonada', 'Discoba', 'Sar'}

    if 'etagenomes' in ','.join(names):
        print(lineage)

    if 'metagenomes' in names :
        return 'metagenomes'
    elif 'Archaea' in names:
        return 'archaea'
    elif 'Bacteria' in names :
        return 'bacteria'
    elif 'Fungi' in names :
        return 'fungi'
    elif 'Viruses' in names:
        return 'viral'
    elif 'Viridiplantae' in names or 'Rhodophyta' in names:
        return 'plant'
    elif protazoa.intersection(set(names)) :
        return 'protazoa'
    elif 'Opisthokonta' in names:
        if 'Vertebrata' in names:
            if 'Mammalia' in names:
                return 'vertebrate_mammalian'
            else :
                return 'vertebrate_other'
        else :
            return 'invertebrate'
    else :
        return 'all'


##########################################################################


def infer_ngs_groups(df_taxid, ncbi):
    """
    Infer taxid ngs option if not in taxid
    """

    if snakemake.config["ndg_option"]["groups"] != 'all':
        if "NCBIGroups" not in df_taxid:
            df_taxid["NCBIGroups"] = snakemake.config["ndg_option"]["groups"].lower()
        elif df_taxid.NCBIGroups.isnull().any():
            df_taxid.loc[df_taxid.NCBIGroups.isnull(), "NCBIGroups"] = snakemake.config["ndg_option"]["groups"].lower()
        else :
            df_taxid["NCBIGroups"] = df_taxid["NCBIGroups"].str.lower()
    else:
        if "NCBIGroups" not in df_taxid:
            for index, row in df_taxid.iterrows():
                df_taxid.at[index, "NCBIGroups"] = find_NCBIgroups(row.TaxId, ncbi)
        elif df_taxid.NCBIGroups.isnull().any():
            sub_df = df_taxid.loc[df_taxid.NCBIGroups.isnull()]
            df_taxid = df_taxid.astype(str)
            for index, row in sub_df.iterrows():
                df_taxid.at[index, "NCBIGroups"] = find_NCBIgroups(row.TaxId, ncbi)

    return df_taxid


##########################################################################

def fetch_seq_assembly(mygroup, assembly_final_df, keyargs):
    """
    Function that will fetch the genome from the taxid in group and 
    concat the assembly that is created by ngs
    """

    keyargs['taxids'] = [str(taxid) for taxid in mygroup.TaxId.tolist()]

    #ngd.download(**keyargs)  
    get_cmdline_ndg(**keyargs)

    # Test if we download genomes else return assembly_final_df
    if os.path.isfile(snakemake.output.assembly_output):
        # Read the information about the assembly and concatenate with previous one
        tmp_assembly = pd.read_table(snakemake.output.assembly_output)

        # Remove the file
        os.remove(snakemake.output.assembly_output)

        return pd.concat([assembly_final_df, tmp_assembly])
    else :
        return assembly_final_df

##########################################################################

def main():
    '''
    Here I define a main function because the library multiprocess 
    bug sometime when it is not this format of script
    '''

    # The soft limit imposed by the current configuration
    soft_limit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]

    # If never done, will download the last version of NCBI Taxonomy dump and create the SQL database
    ncbi = NCBITaxa()  

    # To update the database 
    if snakemake.params.update_db :
        ncbi.update_taxonomy_database()  

    # If wanted the taxdump could be remove after the database is created
    if os.path.isfile('taxdump.tar.gz') :
       os.remove('taxdump.tar.gz')

    taxid_df = pd.read_table(snakemake.input[0])
    taxid_df = infer_ngs_groups(taxid_df, ncbi)

    new_taxid_df = pd.DataFrame()

    ### Amelioration purposes : Seems we could parse line per line the taxid file without loading it
    # For each item in the taxid file, check the lastest descendant taxa from the taxid
    # If a taxid doesn't exist, it is not added in the dataframe
    for index, row in taxid_df.iterrows():
        try :
            list_taxid = ncbi.get_descendant_taxa(row.TaxId)

            tmp_df = pd.DataFrame({'TaxId':list_taxid,
                                   'NCBIGroups': row.NCBIGroups})

            new_taxid_df = pd.concat([new_taxid_df, tmp_df])
        except ValueError :
            pass


    # Drop duplicate in case user have multiple time the same taxonomic id
    new_taxid_df.drop_duplicates('TaxId', inplace=True)        
    new_taxid_df.to_csv(snakemake.output.new_taxid_file, index=False, sep='\t')

    # Dictionnary of all the option of ncbi_genome_download
    keyargs = { "section": snakemake.params.section, 
                "file_formats": "protein-fasta", 
                'flat_output':True,  
                "output": snakemake.params.output_database_folder, 
                "parallel": snakemake.threads, 
                "metadata_table": snakemake.output.assembly_output,
                "assembly_levels":snakemake.params.assembly_levels,
                "refseq_categories":snakemake.params.refseq_categories}       

    # Create a dataframe that will contains all the assembly info
    assembly_final_df = pd.DataFrame()

    # If think it is better to search for all the same NCBI groups in the same time as the function can be parallized
    for NCBIgroup, group in new_taxid_df.groupby('NCBIGroups'):
        # Set the NCBI group
        keyargs['groups'] = NCBIgroup

        # Test if the number of TaxId is inferior to the soft limit imposed by the current configuration
        if soft_limit < group.shape[0]:
            number_of_file = soft_limit - 20
            batch_group = np.arange(group.shape[0])//number_of_file

            for index, batch in group.groupby(batch_group):
                assembly_final_df = fetch_seq_assembly(batch, assembly_final_df, keyargs)
        else :
            assembly_final_df = fetch_seq_assembly(group, assembly_final_df, keyargs)          

    # Remove last column because the file doesn't exist at the end
    assembly_final_df.iloc[:,:-1].to_csv(snakemake.output.assembly_output, index=False, sep='\t')

    # Now handeling the last step that is concatenate the fasta files downloaded
    columns_protein_file =  ['protein_id',
                             'protein_name',
                             'genome_name',
                             'genome_id',
                             'length']

    with open(snakemake.output.output_table_protein, 'wt') as protein_file:
        header_prot = '\t'.join(columns_protein_file)
        protein_file.write(f'{header_prot}\n')

        with open(snakemake.output.fasta_final, 'w') as fasta_file :
            for index, genome in assembly_final_df.iterrows() :
                # Test if file exists before
                if os.path.isfile(genome.local_filename):
                    with gzip.open(genome.local_filename, "rt") as handle:
                        parser = SeqIO.parse(handle, "fasta")

                        for protein in parser : 
                            description_split = protein.description.split(' [')

                            # To avoid duplicate name of proteins between close genomes
                            protein.id = f'{protein.id}--{genome.assembly_accession}'

                            protein_id   = protein.id
                            protein_name = ' '.join(description_split[0].split(' ')[1:])
                            genome_name  = genome.organism_name
                            genome_id    = genome.assembly_accession
                            length       = len(protein.seq)

                            line_prot = f'{protein_id}\t{protein_name}\t{genome_name}\t{genome_id}\t{length}\n'
                            protein_file.write(line_prot)

                            SeqIO.write(protein, fasta_file, 'fasta')

                    # Here to save space I decided to remove the file avec concatenation
                    os.remove(genome.local_filename)
    return
    
##########################################################################

if __name__ == "__main__":
    # Put error and out into the log file
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")
    main()
