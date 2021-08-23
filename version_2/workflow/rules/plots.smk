# Rules for plotting the results

##########################################################################
##########################################################################

rule plots :
    input :
        final_table = os.path.join(OUTPUT_FOLDER, 'results', 'patab.tsv')
    output :
        multiext(os.path.join(OUTPUT_FOLDER,'results','plots', 'gene_PA'), 
                                                         '.png', '.pdf')
    log:
        "logs/blast/plots.log",
    conda:
        "../envs/plots.yaml"
    script :
        "../scripts/plot-PA.py"

##########################################################################
##########################################################################
