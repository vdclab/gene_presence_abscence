# Rule for pluggin threshold

##########################################################################
##########################################################################


rule blast2threshold_table:
    input:
        seed_file=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "new_seeds.tsv"),
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "blastp--blast_evalue_1e-2.out"
        ),
        fnodes=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "silix",
            "modif",
            "filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.fnodes.flushed",
        ),
        protein_file=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "protein_table.tsv"
        ),     
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "analysis_thresholds",
            "table_hits_family--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.tsv"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "analysis_thresholds",
            "blast2threshold_table--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.log"
        ),
    resources:
        mem_mb=8000,
    conda:
        "../envs/pandas_plots.yaml"
    script:
        "../scripts/blast2threshold_table.py"


##########################################################################
##########################################################################


rule report_threshold:
    input:
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "analysis_thresholds",
                "table_hits_family--{gene_constrains}.tsv"
            ),
            gene_constrains=gene_constrains,
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "analysis_thresholds",
            "report_figure_thresholds.html",
        ),
    params:
        css=workflow.source_path("../report/threshold_report.css"),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "analysis_thresholds",
            "analysis_thresholds_fig.log"
        ),
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/blast2threshold_fig.py"


##########################################################################
##########################################################################
