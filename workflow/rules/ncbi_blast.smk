# Module containing all the ncbi-blast related rules

##########################################################################
##########################################################################


rule psiblast:
    input:
        seed=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.fasta"),
        taxid_db=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "taxid_all_together.fasta"
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "psiblast",
            f"psiblast--eval_{e_val_psiblast:.0e}_raw.out",
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "psiblast", "psiblast.log"),
    params:
        e_value=e_val_psiblast,
    resources:
        cpus=5,
        mem_mb=10000,
        time_min=7200,
    conda:
        "../envs/blast.yaml"
    envmodules:
        "ncbi_blast/2.10.1",
    threads: 5
    shell:
        """
        makeblastdb -dbtype prot -in '{input.taxid_db}' -parse_seqids

        psiblast -query '{input.seed}' -db '{input.taxid_db}' -evalue {params.e_value} \
                 -outfmt '7 qacc qlen qseq qstart qend sacc slen sseq sstart send length pident evalue bitscore qcovs' \
                 -num_threads {threads} -num_iterations 3 -out '{output}'

        rm {input.taxid_db}.*
        """


##########################################################################
##########################################################################


rule blast:
    input:
        taxid_fasta=speedup,
        seed_fasta=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.fasta"),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER, "processing_files", "blast", "blastp--blast_evalue_1e-2.out"
        ),
        fasta_for_blast=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "blast", "blast.log"),
    resources:
        cpus=5,
        mem_mb=10000,
        time_min=7200,
    conda:
        "../envs/blast.yaml"
    envmodules:
        "ncbi_blast/2.10.1",
    threads: 5
    shell:
        """
        cat '{input.taxid_fasta}' '{input.seed_fasta}' > '{output.fasta_for_blast}'

        makeblastdb -dbtype prot -in '{output.fasta_for_blast}'

        blastp -query '{output.fasta_for_blast}' -db '{output.fasta_for_blast}' -evalue 0.01 \
               -outfmt 6 -out '{output.blast_out}' -num_threads {threads} -num_alignments 25000

        rm {output.fasta_for_blast}.*
        """
