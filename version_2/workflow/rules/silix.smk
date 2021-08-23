# Rule for silix 

##########################################################################
##########################################################################

rule silix:
    """
    Uses Silix to create a network of protein and give a file of the protein segregated in groups.
    If the blast output file is empty, just create an empty file
    
    Inputs
    ------
    blast_out : str
        blast output filtered for a specific seed from the rule prepare_for_silix.
        format: query id | subject id | percentage of identity | length of match  | mismatch | gapopen |
                query start position | query end position | subject start position | subject end position |
                e-value | bitscore
    fasta : str
        multifasta of proteins with seed from the rule blast
        
    Output : str
    ------
        fnodes file, table of protein id and family number, without headers.
        format: family | protein id
        
    Params
    ------
    silix_version : str
        version of silix to use
    """

    input:
        fasta = os.path.join(OUTPUT_FOLDER, 'database', 'reduce_taxid', 'all_protein_with_seeds.fasta'),
        blast_out = os.path.join(OUTPUT_FOLDER, 'processing_files','blast_out_per_gene', 'filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.out'),
    output:
        os.path.join(OUTPUT_FOLDER, 'processing_files','blast_out_per_gene', 'filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.fnodes')
    log:
        "logs/blast/{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.silix.log",
    group:
        'independent_seed'
    conda :
        '../envs/silix.yaml'
    envmodules:
        "silix/1.2.11"
    shell:
        """
        if [ -s {input.blast_out} ]
        then   
            sh -c 'silix "{input.fasta}" "{input.blast_out}" -f "{wildcards.seed}" -i 0.05 -r 0.05 > "{output}"'
        else
            touch '{output}'
        fi
        """

##########################################################################
##########################################################################
