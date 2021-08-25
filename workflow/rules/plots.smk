# Rules for plotting the results

##########################################################################
##########################################################################


rule plots:
    input:
        final_table=os.path.join(OUTPUT_FOLDER, "results", "patab_melt.tsv"),
    output:
        png=report(
            os.path.join(OUTPUT_FOLDER, "results", "plots", "gene_PA.png"),
            caption="../report/PA_plot.rst",
            category="Plots",
        ),
        pdf=report(
            os.path.join(OUTPUT_FOLDER, "results", "plots", "gene_PA.pdf"),
            caption="../report/PA_plot.rst",
            category="Plots",
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "plots", "plots.log"),
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/plot-PA.py"


##########################################################################
##########################################################################
