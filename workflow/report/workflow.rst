============================
Gene presence-absence report
============================

gene-presence-abscence_ is a Snakemake workflow for ???
. This analysis is based on commit version {{ snakemake.config["__workflow_version__"] }}_.

The analysis can be rerun with the following command:

Local:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --use-conda {{ snakemake.config["__config_args__"] }}
{% else %}
   snakemake -j 1 --use-conda {{ snakemake.config["__config_args__"] }} -s {{ snakemake.config["__workflow_basedir_short__"] }}/Snakefile
{% endif %}

On SLURM cluster:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake --profile profile/slurm -j 1 --use-conda {{ snakemake.config["__config_args__"] }}
{% else %}
   snakemake --profile profile/slurm -j 1 --use-conda {{ snakemake.config["__config_args__"] }} -s {{ snakemake.config["__workflow_basedir_short__"] }}/Snakefile
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
-----------------

.. code-block:: text

   {{ snakemake.config["project_name"] }}/                          <- top-level project folder
   │
   │
   ├── report.html                           <- This report file      
   │
   ├── logs                                  <- Collection of log outputs, e.g. from cluster managers
   │
   ├── databases                             <- Generated analysis database related files
   │   ├── all_taxid                         
   │   │   ├─ protein_table.tsv              <- Table with the informations about the proteins of the downloaded taxid
   │   │   ├─ summary_assembly_taxid.tsv     <- Table with the informations about the downloaded genome from NCBI
   │   │   ├─ taxid_all_together.fasta       <- Fasta file of the downloaded taxid
   │   │   └─ taxid_checked.txt              <- List of the downloaded taxid
   │   │
   │   ├── merge_fasta                       
   │   │   └─ taxid_checked.txt              <- Fasta with the concatenation of the genome of interest and seeds
   │   │
   │   └── seeds                             
   │       ├─ seeds.fasta                    <- Fasta file of the seeds
   │       └─ new_seeds.tsv                  <- Table with the informations about the seeds
   │
   └── results                               <- Final results for sharing with collaborators, typically derived from analysis sets
       ├── patab_melt.tsv                    <- Table with the information of presence-absence one information by line
       ├── patab_table.tsv                   <- Table with the information of presence absence with genome in index and seeds in columns and proteins Id in the cell
       └── plots                             <- Plots and table on which the plot are created



Analysis overview
-----------------

The analyses can basically be divided in two parts: `Raw data
analysis`_ and `Analysis sets`_.

Raw data analysis
*****************

The raw data analysis generates silix results that serve as a
starting point for subsequent analyses.

Analysis sets
*************

Once fnodes data has been generated it can be further analyzed.

For each seeds the workflow analyses in which family it belong.
Find the proteins in the same family. Highlight the family on
the genome


General results
---------------

Figure
******

.. figure:: gene_PA.png_
   :width: 30%
   :align: center


Table format
************

.. csv-table:: Presence Absence in genomes
   :file: {{ snakemake.config["__output_folder__"] }}/results/patab_table.tsv
   :width: 20%
   :delim: tab
   :align: center

Workflow graph
--------------


.. _gene-presence-abscence: https://github.com/vdclab/gene_presence_abscence
.. _{{ snakemake.config["__workflow_version__"] }}: {{ snakemake.config["__workflow_version_link__"] }}
