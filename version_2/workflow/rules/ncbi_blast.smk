# Module containing all the ncbi-blast related rules

##########################################################################
##########################################################################

rule psiblast:
    """
    Use the sequences of the seeds to make a psiBLAST against all the taxid
    
    Inputs
    ------
    seed : str
        the seed multifasta file input from rule fetch_fasta_from_seed
    taxid : str
        list of taxid in columns, no header
        
    Output
    ------
    blast out : str
        blast out format in tabulation, no header
        format : query accession | query length | query sequence | query start position | querry end position |
                subject accession | subject length | subject sequence| subject start position | subject end position | 
                length of alignment | percentage of identity | e-value | bitscore | querry coverage
                
    Params
    ------
    e_val : int
        e-value threshold for psi-blast chosen by user
    blast_version : str
        blast version to use
    """

    input:
        seed = os.path.join(OUTPUT_FOLDER, 'seeds.fasta'),
        taxid_db = os.path.join(OUTPUT_FOLDER, 'database', 'all_taxid', 'taxid_all_together.fasta')
    output:
        os.path.join(OUTPUT_FOLDER, 'processing_files', f'psiblast--eval_{e_val}_raw.out')
    log:
        "logs/blast/psiblast.log",
    params:
        e_val = e_val,
    resources: 
        cpus=5, mem_mb='10gb', time_min='5-0'    
    envmodules:
        "ncbi_blast/2.10.1"
    threads:
        5
    shell:
        '''
        makeblastdb -dbtype prot -in '{input.taxid_db}' -parse_seqids

        psiblast -query '{input.seed}' -db '{input.taxid_db}' -evalue {e_val} \
                 -outfmt '7 qacc qlen qseq qstart qend sacc slen sseq sstart send length pident evalue bitscore qcovs' \
                 -num_threads {threads} -num_iterations 3 -out '{output}'
        '''

##########################################################################
##########################################################################

rule blast:
    """
    blast all versus all of the fasta of all protein generated in the rule make_fasta
    
    Inputs
    ------
    prot_sequence : str
        multifasta file of all the unique protein ids from the rule make_fasta
    seed_fasta : str
        multifasta file of all the seeds from the rule fetch_fasta_from_seed
        
    Outputs
    -------
    blast_out : str
        output format of blast
        format: query id | subject id | percentage of identity | length of match  | mismatch | gapopen |
                query start position | query end position | subject start position | subject end position |
                e-value | bitscore
    fasta_for_blast : str
        concatenation of the 2 input multifasta files
        
    Params
    ------
    blast_version : str
        version of blast
    """

    input:
         taxid_fasta = speedup,
         seed_fasta = os.path.join(OUTPUT_FOLDER, 'seeds.fasta')
    output:
        blast_out = os.path.join(OUTPUT_FOLDER, 'processing_files', 'blastp--blast_evalue_1e-2.out'),
        fasta_for_blast = os.path.join(OUTPUT_FOLDER, 'database', 'reduce_taxid', 'all_protein_with_seeds.fasta')
    log:
        "logs/blast/blast.log",
    resources: 
        cpus=5, mem_mb='10gb', time_min='5-0' 
    envmodules:
        "ncbi_blast/2.10.1"  
    threads:
        5
    shell:
         """
         cat '{input.taxid_fasta}' '{input.seed_fasta}' > '{output.fasta_for_blast}'

         makeblastdb -dbtype prot -in '{output.fasta_for_blast}'

         blastp -query '{output.fasta_for_blast}' -db '{output.fasta_for_blast}' -evalue 0.01 \
                -outfmt 6 -out '{output.blast_out}' -num_threads {threads} -num_alignments 25000
         """