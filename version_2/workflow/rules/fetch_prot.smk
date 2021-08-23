# Rules to fetch proteins fasta files from NCBI

##########################################################################
##########################################################################

rule fetch_proteins_database:
    """
    Create the protein database with all the taxid given by the user. 
    Also providing a tab-delimited file with genome metadata.

    Input : str
    -----
        list of taxid in column, no header, from the rule plit_taxid_file
        
    Output : str
    ------
        table of the information collected on the proteins, without header.
        format: protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
    """

    input:
        taxid
    output:
        fasta_final = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 'taxid_all_together.fasta'),
        assembly_output = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 'summary_assembly_taxid.tsv'),
        new_taxid_file = os.path.join(OUTPUT_FOLDER, 'taxid_checked.txt'),
        output_table_protein = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 'protein_table.tsv')
    log:
        "logs/fetch_proteins_database.log",  
    params:
        output_database_folder = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid'),
        section = section,
        assembly_levels = assembly_levels,
        refseq_categories = refseq_categories,
    resources: 
        cpus=5,
    threads :
        5
    conda:
        "../envs/fetch_NCBI.yaml"
    script :
        "../scripts/fetch_proteome.py"

##########################################################################
##########################################################################

rule fetch_fasta_from_seed:
    """
    from the seed table and fetch the fasta of the seed. Then they are writen in the output file.

    Input : str
    -----
        the seed file input.
        table without header in the format : 
            name | protein id | e-value | percentage of identity | coverage | color
        
    Outputs : str
    -------   
        multifasta output of the seed sequences
    """
    input :
        seed_file
    output:
        os.path.join(OUTPUT_FOLDER, 'seeds.fasta')
    log:
        "logs/fetch_fasta_from_seed.log",  
    conda:
        "../envs/biopython.yaml"
    script :
        "../scripts/fetch_seed.py"


##########################################################################
##########################################################################
