# Rule to make fasta from table

##########################################################################
##########################################################################

rule make_fasta:
    input:
        protein_fasta = os.path.join(OUTPUT_FOLDER, 'databases', 'all_taxid', 'taxid_all_together.fasta'),
        list_all_prot = os.path.join(OUTPUT_FOLDER, 'processing_files', 'psiblast', f'list_all_protein--eval_{e_val_psiblast:.0e}.tsv')
    output:
        fasta = os.path.join(OUTPUT_FOLDER, 'databases', 'reduce_taxid', f'all_protein--eval_{e_val_psiblast:.0e}.fasta'),
    log:
        os.path.join(OUTPUT_FOLDER, 'logs', 'format_table', "make_fasta.log"),        
    conda:
        "../envs/biopython.yaml"
    script :
        "../scripts/hit2fasta.py"


##########################################################################
##########################################################################
