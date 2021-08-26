import pandas as pd
from ete3 import NCBITaxa 
import ncbi_genome_download as ngd 
import os
import gzip
from Bio import SeqIO
from urllib.request import urlopen
import subprocess

##########################################################################

def execfile(script_py, global_vars):
    '''
    Alternative to python2 execfile
    '''

    if '://' not in script_py:
        script_py = script_py.replace(':/', '://')
        global_vars['__file__'] = script_py

    with urlopen(script_py) as f:
        script_txt = script_py.read()
        subprocess.run(['python3', '-'], text=True, input=script_txt)

##########################################################################

def main():
    '''
    Here I define a main function because the library multiprocess 
    bug sometime when it is not this format of script
    '''
    
    print('g = ', g['__file__'])
    print('globals()', globals()['__file__'])

    # If never done, will download the last version of NCBI Taxonomy dump and create the SQL database
    ncbi = NCBITaxa()  

    # To update the database 
    if snakemake.params.update_db :
        ncbi.update_taxonomy_database()  

    # If wanted the taxdump could be remove after the database is created
    if os.path.isfile('taxdump.tar.gz') :
       os.remove('taxdump.tar.gz')

    taxid_df = pd.read_table(snakemake.input[0])

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
        keyargs['groups'] = NCBIgroup
        keyargs['taxids'] = [str(taxid) for taxid in group.TaxId.tolist()]

        ngd.download(**keyargs)  

        # Read the information about the assembly and concatenate with previous one
        tmp_assembly = pd.read_table(snakemake.output.assembly_output)
        assembly_final_df = pd.concat([assembly_final_df, tmp_assembly])

    # Remove last column because the file doesn't exist at the end
    assembly_final_df.iloc[:,:-1].to_csv(snakemake.output.assembly_output, index=False, sep='\t')

    # Now handeling the last step that is concatenate the fasta files downloaded
    df_proteins = pd.DataFrame(columns = ['protein_id',
                                          'protein_name',
                                          'genome_name',
                                          'genome_id',
                                          'length'])

    with open(snakemake.output.fasta_final, 'w') as w_file :
        for index, genome in assembly_final_df.iterrows() :
            with gzip.open(genome.local_filename, "rt") as handle:
                parser = SeqIO.parse(handle, "fasta")

                for protein in parser : 
                    description_split = protein.description.split(' [')

                    # To avoid duplicate name of proteins between close genomes
                    protein.id = f'{protein.id}--{genome.assembly_accession}'

                    df_proteins.at[-1, 'protein_id']   = protein.id
                    df_proteins.at[-1, 'protein_name'] = ' '.join(description_split[0].split(' ')[1:])
                    df_proteins.at[-1, 'genome_name']  = description_split[1].split(']')[0]
                    df_proteins.at[-1, 'genome_id']    = genome.assembly_accession
                    df_proteins.at[-1, 'length']       = len(protein.seq)

                    df_proteins.index += 1

                    SeqIO.write(protein, w_file, 'fasta')



            # Here to save space I decided to remove the file avec concatenation
            os.remove(genome.local_filename)

    df_proteins.to_csv(snakemake.output.output_table_protein, sep='\t', index=False)

##########################################################################

if __name__ == "__main__":
    g = globals().copy()
    print(globals()["__file__"])

    if '/https' in g["__file__"] :
        g["__file__"] = g["__file__"].split('https')[-1]
        g["__file__"] = f'https{g["__file__"]}'
        execfile(g["__file__"], g)
    else :
        main()
