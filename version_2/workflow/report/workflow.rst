gene_presence_abscence workflow results
=======================================

gene-presence-abscence_ is a Snakemake workflow for ???
. This analysis is based on commit version {{ snakemake.config["__workflow_commit__"] }}_.

The analysis can be rerun with the following command:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --use-conda
{% else %}
   snakemake -j 1 --use-conda -s {{ snakemake.config["__workflow_basedir__"] }}/Snakefile
{% endif %}

.. note::

   Since the workflow is still work in progress, make sure to first
   run the commands with the `--dry-run` (`-n`) flag to make sure you
   don't inadvertedly have to regenerate large parts of the results.
   Many workflow dependencies are complex and in particular when
   running smaller parts of the workflow, unexpected things may
   happen.


Workflow summary
----------------

The workflow runs the following steps:

1. Find and Expand taxonomic id to create database
2. Fetching user seeds
3. Psiblast the seeds againts DB to speedup the blast all VS all (optional)
4. Compare all the pairs of sequences using blast
5. Using independent threshold for each seed
6. Highlighting the different seeds on the genomes


Data organization
=================

.. code-block:: text

   {{ snakemake.config["project_name"] }}/            <- top-level project folder
   │
   ├── logs                     <- Collection of log outputs, e.g. from cluster managers
   │
   ├── reports                  <- Generated analyses and articles as html, pdf and more.
   |   ├── qc                   <- QC reports, including multiqc.html
   │   └── figures              <- Graphics for use in reports.
   │
   └── results                  <- Final results for sharing with collaborators, typically derived from analysis sets
       ├── databases            <- All the fasta downloaded or created during the analysis
       ├── processing_files     <- All tables created during analysis
       └── figures              <- Plots and table on which the plot are created



Analysis overview
=================

The analyses can basically be divided in two parts: `Raw data
analysis`_ and `Analysis sets`_.

Raw data analysis
------------------

.. figure:: {{ Rule_basedir }}/report/dag_first_steps.svg
   :width: 30%
   :align: center

   Schematic overview of the first steps.

The raw data analysis generates silix results that serve as a
starting point for subsequent analyses.

Analysis sets
--------------

Once fnodes data has been generated it can be further analyzed.

For each seeds the workflow analyses in which family it belong.
Find the proteins in the same family. Highlight the family on
the genome


General analyses
=================


Workflow graph
==============


.. _gene-presence-abscence: https://github.com/vdclab/gene_presence_abscence
.. _{{ snakemake.config["__workflow_commit__"] }}: {{ snakemake.config["__workflow_commit_link__"] }}
