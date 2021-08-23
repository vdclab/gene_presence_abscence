# Rule for formating tables

##########################################################################
##########################################################################

rule read_psiblast:
    """
    Read the psiBLAST, remove unwanted lines ane extract the list of matches
    
    Inputs
    ------
    psiblast : str
        blast out format in tabulation, no header, from the psiblast rule
        format : query accession | query length | query sequence | query start position | querry end position |
                subject accession | subject length | subject sequence| subject start position | subject end position | 
                length of alignment | percentage of identity | e-value | bitscore | querry coverage
            
    Outputs
    -------
    clean_blast : str
        cleaned blast out format in tabulation, no header
        format : query accession | query length | query sequence | query start position | querry end position |
                subject accession | subject length | subject sequence| subject start position | subject end position | 
                length of alignment | percentage of identity | e-value | bitscore | querry coverage
    list_all_prot : str
        list of all potein identifications gathered in the psiBLAST in column
    """

    input:
        psiblast = os.path.join(OUTPUT_FOLDER, 'processing_files',f'psiblast--eval_{e_val}_raw.out')
    output:
        clean_blast = os.path.join(OUTPUT_FOLDER, 'processing_files',f'psiblast--eval_{e_val}_cleaned.out'),
        list_all_prot = os.path.join(OUTPUT_FOLDER, 'processing_files',f'list_all_protein--eval_{e_val}.csv')
    log:
        "logs/read_psiblast.log",
    conda:
        "../envs/pandas.yaml"
    script :
        "../scripts/format_psiblast.py"


##########################################################################
##########################################################################

rule prepare_for_silix:
    """
    Filter the blast results from the rule blast with the threshold specified for each seed in the seed file.
    Filters include the identity score, coverage and e-value.
    Create one new filtered blast result for each seed.
    
    Inputs
    ------
    fasta : str
        multifasta of proteins with seed from the rule blast
    blast_out : str
        blast output from the rule blast
        format: query id | subject id | percentage of identity | length of match  | mismatch | gapopen |
                query start position | query end position | subject start position | subject end position |
                e-value | bitscore
        
    Output : list of str
    ------
        list of blast output filtered for each seed.
        format: query id | subject id | percentage of identity | length of match  | mismatch | gapopen |
                query start position | query end position | subject start position | subject end position |
                e-value | bitscore
    
    Params
    ------
    new_dir : str
        work directory
    project_name : str
        project name determined by the user
    """

    input:
        fasta = os.path.join(OUTPUT_FOLDER, 'database', 'reduce_taxid', 'all_protein_with_seeds.fasta'),
        blast_out = os.path.join(OUTPUT_FOLDER, 'processing_files', 'blastp--blast_evalue_1e-2.out'),
    output:
        os.path.join(OUTPUT_FOLDER, 'processing_files', 'blast_out_per_gene', 'filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.out'),
    log:
        "logs/{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.prepare_for_silix.log",
    group:
        'independent_seed'
    conda:
        "../envs/biopython_pandas.yaml"
    script :
        "../scripts/prepare_silix.py"

##########################################################################
##########################################################################

rule find_family:
    """
    Find the group of each seed in each individual seed and record it
    
    Input
    -----
    fnodes : str
        fnodes file, table of protein id and family number, without headers from the rule silix.
        format: family | protein id
        
    Output : str
    ------
        updated fnodes with only the family of the seed.
        format: family | protein id | seed
    """

    input:
        fnodes = os.path.join(OUTPUT_FOLDER, 'processing_files','blast_out_per_gene', 'filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.fnodes'),
        seed_file = seed_file
    output:
        os.path.join(OUTPUT_FOLDER, 'processing_files','blast_out_per_gene', 'filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.fnodes.flushed')
    log:
        "logs/{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.find_family.log",
    group:
        'independent_seed'
    conda:
        "../envs/pandas.yaml"
    script :
        "../scripts/find_family.py"

##########################################################################
##########################################################################

rule make_table:
    """
    Check the presence of protein similar to the seed in each taxid and create a table of presence abscence
    This table is then plotted in a colored table.
    
    Inputs
    ------
    protein_table : str
        final table of protein information from the rule cat_proteins_info, without header.
        format: protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
    fnode : str
        concatenated fnodes with each seed family from 
        format: family | protein id | seed
            
    Outputs
    -------
    final_table : str
        presence/abscence table, with header. Each line is a genome, each column is a seed.
        format: genome id | genome name | seed 1 | seed 2 .. seed x
    pdf : list of str
        plots in pdf of the final table centered on one seed
    png : list of str
        plots in png of the final table centered on one seed
    """

    input:
        seed_file = seed_file,
        protein_table = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 'protein_table.tsv'),
        assembly_table = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 'summary_assembly_taxid.txt'),
        fnodes = expand(os.path.join(OUTPUT_FOLDER, 
                                    'processing_files','blast_out_per_gene', 
                                    'filtered_blast--{gene_constrains}.fnodes.flushed'),
                                        gene_constrains=gene_constrains)
    output:
        final_table = os.path.join(OUTPUT_FOLDER, 'results', 'patab.tsv'),
    log:
        "logs/make_table.log",
    conda:
        "../envs/pandas.yaml"
    script :
        "../scripts/PA-table.py"

##########################################################################
##########################################################################
