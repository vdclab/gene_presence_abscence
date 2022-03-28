# Rules to deal with or create HMM profiles

##########################################################################
##########################################################################


rule hmmsearch:
    input:
        hmm=list_hmm,
        taxid_db=merge_db,
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmsearch",
            f"HMM--eval_{e_val_HMM:.0e}.out",
        ),
    params:
        hmm_type=hmm_type,
        e_value=e_val_HMM,
        seed_hmm=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "hmm_db",
            "merged_HMM.hmm",
        ),
    resources:
        cpus=5,
        mem_mb=10000,
        time=7200,
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "hmmer", "hmmsearch.log"),
    threads:
        5
    conda:
        "../envs/hmmer.yaml"
    envmodules:
        "hmmer/3.3.2",
    shell:
        """
        cat {input.hmm} > '{params.seed_hmm}' 2> '{log}'

        hmmsearch --domtblout {output} {params.hmm_type} {params.e_value} --cpu {threads} \
            {params.seed_hmm} {input.taxid_db} 2>> '{log}'
        """


##########################################################################
##########################################################################