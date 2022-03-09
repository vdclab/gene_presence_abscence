# Rules to deal with or create HMM profiles


##########################################################################
##########################################################################


rule hmm_merger:
    input:
        hmm=expand(
            os.path.join(hmm_folder,"{hmm_profiles}"),
            hmm_profiles=hmm_profiles
        ),
        seeds=start_seed_file,
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "HMM",
            f"merged_HMM.hmm",
        ),
    log:
        os.path.join(OUTPUT_FOLDER,"logs","hmm","hmm_merger.log"),
    conda:
        "../envs/hmm.yaml"
    shell:
        """
        cat {input.hmm} > {output}
        """


##########################################################################
##########################################################################


rule hmmsearch:
    input:
        seed_hmm=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "HMM",
            f"merged_HMM.hmm",
        ),
        taxid_db=os.path.join(
                    OUTPUT_FOLDER, "databases", "all_taxid", "taxid_all_together.fasta"
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "HMM",
            f"HMM--eval_{e_val_HMM:.0e}.out",
        ),
    params:
        hmm_type = hmm_type,
        e_value = e_val_HMM,
    resources:
        cpus=5,
        mem_mb=10000,
        time=7200,
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "hmm", "hmmsearch.log"),
    threads:
        5
    conda:
        "../envs/hmm.yaml"
    envmodules:
        "hmmer/3.3.2",
    shell:
        """
        hmmsearch --domtblout {output} {params.hmm_type} {params.e_value} --cpu {threads} \
            {input.seed_hmm} {input.taxid_db} > NOT_IN_LOG_VAR
        """


##########################################################################
##########################################################################
