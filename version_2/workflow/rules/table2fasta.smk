# Rule to make fasta from table

##########################################################################
##########################################################################

rule make_fasta:
    """
    Create a fasta file from the psiblast results and the result of the protein information in the rule cat_proteins_info
    
    Inputs
    ------
    protein_table : str
        final table of protein information from the rule cat_proteins_info, without header.
        format: protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
    list_all_prot : str
        list of all protein identifications gathered in the psiBLAST in column
        
    Outputs 
    -------
    fasta : str
        multifasta file of all the unique protein ids.
    reduced_protein_table : str
        final table of protein information with removed duplicates, without header.
        format: protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
    """
    input:
        protein_fasta = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 'taxid_all_together.fasta'),
        list_all_prot = os.path.join(OUTPUT_FOLDER, 'processing_files',f'list_all_protein--eval_{e_val}.tsv')
    output:
        fasta = os.path.join(OUTPUT_FOLDER, 'results', f'all_protein--eval_{e_val}.fasta'),
    log:
        "logs/make_fasta.log",        
    conda:
        "../envs/biopython.yaml"
    script :
        "../scripts/hit2fasta.py"


##########################################################################
##########################################################################
