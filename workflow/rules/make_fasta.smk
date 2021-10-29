# Rule to make fasta from table

##########################################################################
##########################################################################


rule make_fasta:
    input:
        protein_fasta=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "taxid_all_together.fasta"
        ),
        list_all_prot=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "psiblast",
            f"list_all_protein--eval_{e_val_psiblast:.0e}.tsv",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "reduce_taxid",
            f"all_protein--eval_{e_val_psiblast:.0e}.fasta",
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "format_table", "make_fasta.log"),
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/hit2fasta.py"


##########################################################################
##########################################################################


rule extract_protein:
    input: 
        PAtab=os.path.join(OUTPUT_FOLDER, "results", "patab_melt.tsv"),
        fasta_protein=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "taxid_all_together.fasta"
        ),        
    output:
        expand(
            os.path.join(OUTPUT_FOLDER, "results", "fasta", "{seed}.fasta"),
            seed = seed_table.seed,
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fetch_proteins",
            "extract_protein.log",
        ),
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/extract_protein.py"


##########################################################################
##########################################################################


rule merge_fasta:
    input: 
        taxid_fasta=speedup,
        seed_fasta=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.fasta"),        
    output:
        fasta_for_blast=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fetch_proteins",
            "extract_protein.log",
        ),
    shell:
        """cat '{input.taxid_fasta}' '{input.seed_fasta}' > '{output.fasta_for_blast}'"""


##########################################################################
##########################################################################
