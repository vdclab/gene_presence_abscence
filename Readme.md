# Snakemake workflow: sORTholog

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.14.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/vdclab/sORTholog/workflows/Tests/badge.svg?branch=main)](https://github.com/vdclab/sORTholog/actions?query=branch%3Amain+workflow%3ATests)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Aim

This software will produce a table visualization of the presence or absence of encoded proteins in a given set of genomes. Both, the proteins and set of proteomes are required as input by the user.

The primary output will be the table of sORTholog protein data filled using the NCBI Protein Database, as well as PDF and PNG files of the table visualization.

## Table of content

- [Snakemake workflow: sORTholog](#snakemake-workflow-sortholog)
  - [Aim](#aim)
  - [Table of content](#table-of-content)
  - [Installation and Running](#installation-and-running)
    - [Step 1: install Snakemake and Snakedeploy](#step-1-install-snakemake-and-snakedeploy)
    - [Step 2: deploy workflow](#step-2-deploy-workflow)
    - [Step 3: configure workflow](#step-3-configure-workflow)
      - [Taxid file](#taxid-file)
      - [Seed file](#seed-file)
      - [Config file](#config-file)
    - [Step 4: run the workflow](#step-4-run-the-workflow)
    - [Step 5: generate report](#step-5-generate-report)
  - [Important Notes for running on the UFL cluster](#important-notes-for-running-on-the-ufl-cluster)
    - [Step C1: log on the cluster [UFL users]](#step-c1-log-on-the-cluster-ufl-users)
    - [Step C2: load conda](#step-c2-load-conda)
    - [Step C3: configure sbatch jobs (Need to be changed if you are not from VDCLab)](#step-c3-configure-sbatch-jobs-need-to-be-changed-if-you-are-not-from-vdclab)
    - [Step C4: run the workflow](#step-c4-run-the-workflow)
    - [Additional information](#additional-information)
  - [Walk-Through and File Production](#walk-through-and-file-production)
    - [Pipeline in image](#pipeline-in-image)
    - [Rule desciptions](#rule-desciptions)
      - [Rule 1: fetch_fasta_from_seed](#rule-1-fetch_fasta_from_seed)
      - [Rule 2: fetch_proteins_database](#rule-2-fetch_proteins_database)
      - [Rule 2.1.1: make_seed_psiblast (optional)](#rule-211-make_seed_psiblast-optional)
      - [Rule 2.1.2: psiblast (optional)](#rule-212-psiblast-optional)
      - [Rule 2.1.3: read_psiblast (optional)](#rule-213-read_psiblast-optional)
      - [Rule 2.2.1: hmmsearch (optional)](#rule-221-hmmsearch-optional)
      - [Rule 2.2.2: read_hmmsearch (optional)](#rule-222-read_hmmsearch-optional)
      - [Rule 2.3: merge_hmmsearch_psiblast (optional)](#rule-23-merge_hmmsearch_psiblast-optional)
      - [Rule 2.4: make_fasta (optional)](#rule-24-make_fasta-optional)
      - [Rule 3: merge_fasta](#rule-3-merge_fasta)
      - [Rule 4: blast](#rule-4-blast)
      - [Rule 5: prepare_for_silix](#rule-5-prepare_for_silix)
      - [Rule 6: silix](#rule-6-silix)
      - [Rule 7: find_family](#rule-7-find_family)
      - [Rule 8: make_PA_table](#rule-8-make_pa_table)
      - [Rule 9: plots](#rule-9-plots)
    - [Additional rules](#additional-rules)
      - [Rule A1: report_threshold](#rule-a1-report_threshold)
      - [Rule A2: extract_protein](#rule-a2-extract_protein)
      - [Rule A3: quick_plots](#rule-a3-quick_plots)
      - [Rule A4: clean](#rule-a4-clean)

## Installation and Running

### Step 1: install Snakemake and Snakedeploy

This pipeline uses the Snakemake workflow, that is only available on unix environment. Any linuxOS or MacOS should be able to run sORTholog. Windows users could use the workflow using the [Ubuntu virtual environment](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview).

Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).

**Make sure that mamba is properly configure in your PATH**

To install Mamba by conda, run

```shell
conda install mamba -n base -c conda-forge
```

Given that Mamba is installed, run 

```shell
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
```

to install both Snakemake and Snakedeploy in an isolated environment. 

**Notes**

For all following commands (step 2 to 5) ensure that this environment is activated via 

```shell
conda activate snakemake
```
:warning: **This command has to be used every single time you want to use sORTholog to be able to activate the conda environment that contains snakemake.**

### Step 2: deploy workflow

Given that Snakemake and Snakedeploy are installed and available (see [Step 1: install Snakemake and Snakedeploy](#step-1-install-snakemake-and-snakedeploy)), the workflow can be deployed as follows.

First, create an appropriate working directory, in which sORTholog will be deployed, on your system in the place of your choice as follow (note to change the path and file name to the one you want to create):

```shell
mkdir path/to/sORTholog
```
NB: [More on the mkdir command.](https://en.wikipedia.org/wiki/Mkdir)

Then go to your sORTholog working directory as follow:

```shell
cd path/to/sORTholog
```

NB: [More on the cd command.](https://en.wikipedia.org/wiki/Cd_(command))

Second, run (if your are in the sORTholog directory)

```shell
snakedeploy deploy-workflow https://github.com/vdclab/sORTholog . --tag 0.4.7
```

else, run 

```shell
snakedeploy deploy-workflow https://github.com/vdclab/sORTholog path/to/sORTholog --tag 0.4.7
```

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder.

[▲ Back to table of content ▲](#table-of-content)

### Step 3: configure workflow

The three files described in this section ([Taxid file](#taxid-file), [Seed file](#seed-file), and [Config file](#config-file)) can be edited with a text editor. We recommend [VScode](https://code.visualstudio.com/) or [SublimeText](https://www.sublimetext.com/) if you do not have already a favorite text editor.

#### Taxid file

The taxid file is a table, in the format of a tabulated text file (e.g. .txt, .tsv). This table should contain 2 columns: `TaxId` and `NCBIGroups`. 
The `TaxId` columns should be filled with the TaxId of the genome you want to analyze. You can use any clade TaxId, sORTholog will take care of downloading all the genome of this clade.
The `NCBIGroups` can be left blanked, but knowing this group would make the download of genomes faster. Here a list of accepted terms: `all`, `archaea`, `bacteria`, `fungi`, `invertebrate`, `metagenomes`, `plant`, `protozoa`, `vertebrate_mammalian`, `vertebrate_other`, `viral`. This will allow to look in the right folder in the FTP of NCBI to download the proteome. By choosing `all`, it will look in all folders and slow down the process.
Ideally your taxid file should be in your sORTholog working directory or in the `config` folder. 
Example provided in the file [doc/dummy_taxids.tsv](https://github.com/vdclab/sORTholog/blob/main/doc/dummy_taxids.tsv) folder in the GitHub page.

#### Seed file

The seed file contains the list of proteins you want to identify in your genomes. It is a table, in the format of a  tabulated text file (e.g. .txt, .tsv). Ideally your seed file should be in your `workdir` folder or in the `config` folder. 
Here is a list of collumns for this file:
- <u>seed</u>: Mandatory, name of the protein as you want it to appear in the final figure.
- <u>protein_id</u>: Mandatory, either the NCBI protein id of a protein you are sure have the function you are searching for.
- <u>hmm</u>: name of the hmm profile file associated with the protein family, including the extension, often `.hmm`. (no default, if not mentioned, psiBLAST will be triggered instead)
- <u>evalue</u>: Number between 0 and 1, e-value threshold of blast.
- <u>pident</u>: Number between 0 and 1, percentage of identity (expressed in frequency) threshold of blast.
- <u>cov</u>: Number between 0 and 1, coverage threshold (expressed in frequency) for blast results. The coverage is based on the coverage of the query.
- <u>color</u>: hexadecimal code of the color you want the positive results to be in the final figure for this seed.
Example provided in the file [doc/dummy_seeds.tsv](https://github.com/vdclab/sORTholog/blob/main/doc/dummy_seeds.tsv) folder in the GitHub page.

**Common mistakes**

- Columns in the seed files should be separated by tabulation, introduction of spaces or other hidden characters might disrupt the reading of the file.
- Do not introduce empty lines in the seed file.
- Make sure to have the correct extension for your hmm profile.

#### Config file

This file is in the `config` folder, and is named `config.yaml`. You can edit this file following these instructions (also in comments in the `config.yaml` file)

**General settings**

- <u>seed</u>: absolute or relative path and name of your seed file.
- <u>taxid</u>: absolute or relative path and name of your taxid file.
- <u>output_folder</u>: Name of the folder in which you want to write the project folder. You can use absolute path or relative path.
- <u>project_name</u>: Name of your project. This name will be used as a folder name to put all the files produced by sORTholog.

NB: For windows user, unix uses forward slashes `/`. Make sure to use the latter. 

**Options to download proteomes**

- ndg_option:
  - <u>section</u>: `refseq` (default) or `genbank`. Database from which you are pulling the proteomes. 
  - <u>assembly_levels</u>: `all` (default), `complete`, `chromosome`, `scaffold` or `contig`. Type of assembly you want.
  - <u>refseq_categories</u>: `reference`, `all` (default). Only refseq reference genomes or all of refseq.
  - <u>groups</u>: `all` (default), `archaea`, `bacteria`, `fungi`, `invertebrate`, `metagenomes`, `plant`, `protozoa`, `vertebrate_mammalian`, `vertebrate_other`, `viral`. Set this group if unspecified in the `NCBITaxId` columns of the [Taxid file](#taxid-file).
- <u>update_db</u>: `True` or `False` (default). Update the Taxonomy dump, will increase run time.

**Options to add your personal proteome**

- <u>perso_database</u>: path to personal proteome database. It should consist of a multifasta file with all the proteins you want to add to the search. 
- <u>perso_annotation</u>: table in tabulated text format that contains the information of the annotation of the fasta file in `perso_database`, columns: protein_id, genome_id[, genome_name]

NB: For windows user, unix uses forward slashes `/`. Make sure to use the latter. 

**Speed up options**

- <u>speedup</u>: `True` or `False`. Define if you want to use psi-blast or hmmsearch to filter your proteome before running a blast all vs all. It is recommended on very large amount of genomes.

**Psiblast options**

  - <u>psiblast_e_val</u>: Number between 0 and 1, default `0.01`. E-value used to accept a result in psi-blast. It is recommended to use a high e-value to gather super family related protein rather than being too stringent in this step. 
  - <u>iteration</u>: Number between 0 and 5, default `5`. Number of iteration of the psi-blast. The more iterations, the more you gather phylogenetically distant proteins.

**HMMsearch options**

- <u>hmm_profiles</u>: path and name of your folder containing all the hmm profiles mentioned in your seed file.
- <u>e_val</u>: Number between 0 and 1, default `0.0001`. E-value default threshold if left empty in the seed file.
- <u>focus</u>: `domain` or `full` (default: full). Analyse the results over the full size of the query or based on domain detection.

**Analysis options**

**Default blast options**

  - <u>filter</u>: `e_value`, `score` or `both` (default: e_value). Se up if you want to filter your blast results by e-value, bit score or both. 
  - <u>e_val</u>: Number between 0 and 1, default `0.0001`. E-value default threshold if left empty in the seed file.
  - <u>pid</u>: Number between 0 and 1, default `0.35`. Percentage of identity (expressed in frequency) default threshold if left empty in the seed file. 
  - <u>cov</u>: Number between 0 and 1, default `0.8`. Coverage of the query (expressed in frequency) default threshold if left empty in the seed file.
  
**Silix options**

  - <u>cov_min</u>: `mean`, `subject`, `query`, `shortest` or `longest` (default: mean). Source or the length divider to calculate the coverage (from query, subject or an average of both)
  - <u>pid_min</u>: `mean`, `subject`, `query`, `shortest`, `longest` or `HSP` (default: mean). Source or the length divider to calculate the percentage of identity (from query, subject or an average of both)
  - <u>length_min</u>: positive number, default 100. Mininimum length to accept partial sequences in families
  
**Plot settings**

default_values_plot:
- <u>color</u>: hexadecimal color, (default: `#131516`)
- <u>colored_border</u>: `True` or `False`. Turn the color of the border a darker shade inferred from the background.

Here a comparison of the two behaviours:

<p align="center">
  <img src="doc/colored_border_option.png?raw=true">
</p>    

- <u>round_border</u>: `True` or `False`. Turn border to roundish shape.

Here a comparison of the two behaviours:

<p align="center">
  <img src="doc/round_border_option.png?raw=true">
</p>

**Only plot table options**

- <u>PAtab_table</u>: Path of the table you want to use to transform it into the plot figure. The table has to be a tabulated text file of a table that has the seeds in columns, the genomes in lines. (The option is for the use of the rule `A3: quick_plots` only)

[▲ Back to table of content ▲](#table-of-content)

### Step 4: run the workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda, run Snakemake with 

```shell
snakemake --cores 1 --use-conda 
```

Here the number of used for the workflow is set to 1 but the number could be increase if needed.

If you want to run the workflow with the additional step to speedup the analysis that contains a big dataset you can either:
- Change in the `config.yaml` the parameter `speedup` and change the value to `True`
- Run the workflow adding `-C speedup=True` to the command line as follow

```shell
snakemake --cores 1 --use-conda -C speedup=True
```

You will find the table of presence of absence `patab_table.tsv`, the melt version of this table `patab_melt.tsv` and a folder `plots` containing the figure representing the table.

[▲ Back to table of content ▲](#table-of-content)

### Step 5: generate report

After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser via 

```shell
snakemake --report report.zip --report-stylesheet config/report.css
```
The resulting report.zip file can be passed on to collaborators, provided as a supplementary file in publications.

[▲ Back to table of content ▲](#table-of-content)

## Important Notes for running on the UFL cluster

### Step C1: log on the cluster [UFL users]

1) To work on the HiPerGator server, user must first open a terminal. 

User must access the server by entering the following command:  

```
ssh (GatorLink ID)@hpg.rc.ufl.edu
```

User will be required to enter their GatorLink account password, then asked to be verified by dual push.
      
2) User will require access to files on the server. To access these files via command line, use the following: 

```
cd /blue/(your group)/(GatorLink ID)
```

You can access to the files of this folder on your computer by adding a server. The address for Windows is as follow: `\\exasmb.rc.ufl.edu`, for Unix (Linux or MacOS) the address: `smb://exasmb.rc.ufl.edu`.  
More info [here](https://help.rc.ufl.edu/doc/Samba_Access).

### Step C2: load conda

Before starting the installation of the program, you will need to load `conda`; this required software is listed as follows:

```
module load conda
```

:warning: **This step will be necessary each time you want to use sORTholog on HiPergator.**

After that the step are the same as in [Step 1: install Snakemake and Snakedeploy](#step-1-install-snakemake-and-snakedeploy), [Step 2: deploy workflow](#step-2-deploy-workflow), and [Step 3: configure workflow](#step-3-configure-workflow).

[▲ Back to table of content ▲](#table-of-content)

### Step C3: configure sbatch jobs (Need to be changed if you are not from VDCLab)

Before running the workflow, another config file need to be assessed to configure the sbatch jobs on HiPerGator. This file is `config/slurm/cluster-config.yaml`

Here are the options:

- <ins>account:</ins> account linked to your HiPerGator. (insert your account here)
- <u>qos:</u> qos account name for HiPerGator. (insert your account here)
- <u>time:</u> time per rule in minutes.
- <u>nodes:</u> number of nodes to use per rules.
- <u>ntasks:</u> number of simultaneous tasks.
- <u>cpus-per-task:</u> number of CPU per tasks on HiPerGator.
- <u>mem:</u> memory allocated per rule if not indicated elsewhere.
- <u>output:</u> destination of the output log.
- <u>error:</u> destination of the error log.
- <u>job-name:</u> name of the job on HiPerGator.

Only the account and qos options need to be changed, the other values can be left as defaults.

[▲ Back to table of content ▲](#table-of-content)

### Step C4: run the workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda and running on the cluster, run Snakemake with 

```shell
snakemake -j 5 --use-conda --profile config/slurm
```

As previously, if you want to run the workflow with the additional step to speedup the analysis that contains a big dataset you can either:
- Change in the `config.yaml` the parameter `speedup` and change the value to `True`
- Run the workflow adding `-C speedup=True` to the command line as follow:

```shell
snakemake -j 5 --use-conda --profile config/slurm -C speedup=True 
```

As previously you can also generate a report following [Step 5: Generate report](#step-5-generate-report).

[▲ Back to table of content ▲](#table-of-content)

### Additional information

If you want to run it as a job,

1. Make sure that the pipeline is deploy where you are and that you change the `config/config.yaml` accordingly to your needs
2. Create a file name, for exemple `my_project.sh`, and copy and paste the following, changing your e-mail address and other relevant parameters:

```bash
#!/bin/sh
#SBATCH --mail-user=ufid@ufl.edu              # Where to send mail
#SBATCH --mail-type=END                       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --job-name=sbatch_snakemake           # Name of the job in slurm 
#SBATCH --cpus-per-task=1                     # Use 1 core
#SBATCH --mem=4G                              # Memory allocated to run the workflow
#SBATCH --time=30-0:00:00                      # Time for the workflow to run, format DD-HH:MM:SS
#SBATCH --nodes=1                             # Number of node to start the workflow
#SBATCH --output=logs/sbatch_snakemake-%j.out # Output of the sbatch job
#SBATCH --error=logs/sbatch-snakemake-%j.err  # Error log of the sbatch job, including the snakemake workflow notifications

echo "Date start $(date)"

# activate a specific conda environment, if you so choose
module load conda
conda activate snakemake 

# make things fail on errors
set -o nounset
set -o errexit
set -x

### run your commands here!
snakemake -j 5 --use-conda --profile config/slurm
snakemake --report report.zip --report-stylesheet config/report.css

echo "Date end $(date)"
```

3. run the following command

```bash
sbatch my_project.sh
```

[▲ Back to table of content ▲](#table-of-content)

## Walk-Through and File Production

This pipeline consists of 9/12 steps called rules that take input files and create output files. Here is a description of the pipeline.

1. As snakemake is set up, there is a last rule, called `all`, that serves to call the last output files and make sure they were created.

2. A folder containing your work will be created:

```
   [project_name]/                           <- top-level project folder (your project_name)
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
   ├── analysis_thresholds                   <- (Optional) Generate by the rule report_thresholds
   │   ├── tables                         
   │   │   └─ table--seed.tsv                <- (Optional) Table with the information of each pair of hits in for the seed family (one per seed)
   │   │
   │   └── report_figure_thresholds.html     <- (Optional) HTML report with the plots made by report_thresholds                  
   │
   └── results                               <- Final results for sharing with collaborators, typically derived from analysis sets
       ├── fasta                             <- (Optional) folder with the fasta file of the orthologs of the seeds
       ├── patab_melt.tsv                    <- Table with the information of sORTholog one information by line
       ├── patab_table.tsv                   <- Table with the information of presence absence with genome in index and seeds in columns and proteins Id in the cell
       └── plots                             <- Plots and table on which the plot are created

```

3. Before starting the pipeline, your taxids will be checked and updated for the lowest descendant levels. This will create a file of checked taxids. It will skip this step if this file already exists.

4. When restarting the pipeline, the software will check if you made any changes in the seed file before running. If changes have been made, it will run what is necessary, else nothing will happen.

[▲ Back to table of content ▲](#table-of-content)

### Pipeline in image 

**Only blast and silix behavior**

<p align="center">
  <img src="doc/dummy_dag.png?raw=true" height="400">
</p>

**Speedup behavior**

<p align="center">
  <img src="doc/dummy_dag_speedup.png?raw=true" height="500">
</p>

[▲ Back to table of content ▲](#table-of-content)

### Rule desciptions

#### Rule 1: fetch_fasta_from_seed

Rule description: Fetch the fasta protein sequence of the seed from the seed table. Then they are writen in the output file.

```
· input file: 
    - seed file: type = str
                 format = seed name | protein id | hmm | e-value | percentage of identity | coverage | color
                 description: seed file input by user

· output files: 
    - multi fasta file: type = str
	                    description = multifasta output of the seed protein sequences 
    - new seed file: type = str
                     columns = seed name | protein id | hmm | e-value | percentage of identity | coverage | color      
                     description = Same file as input but the 'protein id' of each file is changed to match the fasta file                
```

#### Rule 2: fetch_proteins_database

Rule description: Fetch the information for each protein of each genome in the taxid list. That includes: the protein ncbi id, sequence, length and annotation, as well as in which genome it is found. Information for the genome include genome ncbi id, name, taxid and if complete or partial.
 
```
· input file: 
    - taxid file: type = str
                  description = list of taxid in column, second additional column for teh Taxid NCBI group

  
· output files: 
    - multi fasta file: type = str
                        description = multifasta output of the genome proteomes in one file 
    - table of genomes informations: type = str
                                     format = assembly_accession | bioproject | biosample | wgs_master | excluded_from_refseq | refseq_category | relation_to_type_material | taxid | species_taxid | organism_name | infraspecific_name | isolate | version_status | assembly_level | release_type | genome_rep | seq_rel_date | asm_name | submitter | gbrs_paired_asm | paired_asm_comp | ftp_path
                                     description = table of the information collected from the NCBI Assembly.
    - table of protein informations: type = str
                                     format = protein id | protein name | genome name | genome id | length
                                     description = table of the information collected on the proteins.
    - new taxid file: type = str
                      description = list of taxid downloaded in column
                      
```

#### Rule 2.1.1: make_seed_psiblast (optional)

Rule description: filter the seed fasta to keep only the seed that do not have a hmm profile associated with them in the seed initial seed file

```
· input file:
    - fasta seed: type = str
                  description = multifasta of the seed sequences
                  
· output file:    
    - fasta: type = str
             description = multifasta of the seed sequences that do not have a hmm profile associated with 
             
```

#### Rule 2.1.2: psiblast (optional)

Rule description: Use the seed protein sequences to make a psi-BLAST against all the taxid proteome. 
 
```
· input files:
    - seed fasta: type = str
                  description = multifasta of the seed sequences that do not have a hmm profile associated with
    - taxid fasta: type = str 
                   description = multifasta of the genome sequences
        
· output file:  
    - blast out: type = str
                 format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                 description = blast out format 6 in tabulation, no header
```

#### Rule 2.1.3: read_psiblast (optional)

Rule description: Read the psiBLAST output, remove unwanted lines and extract the list of protein matches. 

```
· input files: 
    - psiblast output: type = str
                       format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                       description = blast out format 6 in tabulation, no header
               
· output files:
    - list all prot: type = str
                     description = list of all potein identifications gathered in the psiBLAST in column
```

#### Rule 2.2.1: hmmsearch (optional)

Rule description: Use hmm profile(s) provided to fetch similar proteins in the taxid proteome using hmmsearch from hmmer 

```
· input file: 
    - hmm: type = list of str
           description = list of hmm profile provided in the hmm folder and indictaed in the seed table
    - taxid db: type = str
                description = multifasta of the taxid proteome
               
· output file:
    - domtblout: type = str
                 format = protein id | protein accession | protein length | query | query accession | query length | full evalue | full score | full bias | domain number | domain of | domain c-evalue | domain c-evalue | domain score | domain_bias | query from | query to | ali from | ali to | env from | env to | target description
                 description = domtblout output format of hmmsearch
    - seed hmm: type = str
                description = merge file of the hmm profile provided

```

#### Rule 2.2.2: read_hmmsearch (optional)

Rule description: Read the hmmsearch output and extract the list of protein matches. 

```
· input file: 
    - hmm out: type = str
               format = protein id | protein accession | protein length | query | query accession | query length | full evalue | full score | full bias | domain number | domain of | domain c-evalue | domain c-evalue | domain score | domain_bias | query from | query to | ali from | ali to | env from | env to | target description
               description = domtblout output format of hmmsearch
               
· output file:
    - list all prot: type = str
                     description = list of all potein identifications gathered in the hmmsearch in column
```

#### Rule 2.3: merge_hmmsearch_psiblast (optional)

Rule description: merge the list of proteines from both the psiBLAT and hmmsearch results. 

```
· input files: 
    - list hmmsearch: type = str
                      description = list of all potein identifications gathered in the hmmsearch in column
    - list psiblast: type = str
                     description = list of all potein identifications gathered in the psiBLAST in column             
               
· output file:
    - list all prot: type = str
                     description = list of all potein identifications gathered in both the psiblast and the hmmsearch in column
```

#### Rule 2.4: make_fasta (optional)

Rule description: Create a fasta file from the psiblast results and the result of the protein information in the rule cat_proteins_info. 

```
· input files:
    - protein table: type = str
                     format = protein id | protein name | genome name | genome id | length
                     description = table of protein information 
    - list all prot: type = str
                     description = list of all potein identifications gathered in both the psiblast and the hmmsearch in column
  
· output file: 
    - multifasta: type = str
                  description = multifasta file of all the unique protein ids.
```

#### Rule 3: merge_fasta

Rule description: Merge the proteome of interest with the seeds 

```
· input files:
    - taxid fasta: type = str
                     description =   all proteome multifasta or multifasta file of all the unique protein ids if speedup true
    - seed fasta: type = str
                  description = multifasta file of all the seeds from the rule fetch_fasta_from_seed
  
· output file:
    - fasta for blast: type = str
                       description = concatenation of the 2 input multifasta files
```

#### Rule 4: blast

Rule description: Blast all versus all of the fasta of all proteins. 

```
· input file:
    - fasta for blast: type = str
                       description = concatenation of the the proteome taxid and seeds multifasta files
  
· output file:
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description = output format 6 of blast, no header
```

#### Rule 5: prepare_for_silix

Rule description: Filter the blast results from the rule blast with the threshold specified for each seed in the seed file. Filters include the identity score, coverage and e-value, decided by the user. Create one new filtered blast result for each seed.

```
· input files: 
    - seed file: type = str
                 format = seed name | protein id | hmm | e-value | percentage of identity | coverage | color
                 description: seed file input by user           
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description =  blast output 6 from the rule blast  
  
· output files: 
    - blast output to send to silix: type = list of str
                                     format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                                     description = list of blast output filtered for each seed.
                                   
```

#### Rule 6: silix

Rule description: Uses Silix to create a network of protein and give a file of the protein segregated in groups.  If the blast output file is empty, just create an empty file. 

```
· input files:
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description = blast output filtered for a specific seed from the rule prepare_for_silix.
    - fasta: type = str
             description = multifasta of proteins with seed from the rule merge fasta
  
· output file:
    - fnodes: type = str 
              format = family | protein id
              description = fnodes file, table of protein id and family number, without headers.        
```

#### Rule 7: find_family
 
Rule description: Find the group of each seed in each individual seed and record it. 

```
· input files:
    - fnodes: type = str
              format = family | protein id
              description = fnodes file, table of protein id and family number, without headers from the rule silix.
    - seed file: type = str
                 format = seed name | protein id | hmm | e-value | percentage of identity | coverage | color
                 description: seed file input by user  
  
· output file:
    - updadted fnodes: type = str
                       format = family | protein id | seed
                       description = updated fnodes with only the family of the seed.
```

#### Rule 8: make_PA_table

Rule description: Check the presence of protein similar to the seed in each taxid and create a table of presence abscence. This table will be plotted in a colored table.

```
· input files:
    - seed file: type = str
                 format = seed name | protein id | hmm | e-value | percentage of identity | coverage | color
                 description: seed file input by user  
    - protein table: type = str
                     format =  protein id | protein name | genome name | genome id | length
                     description = final table of protein information
    - updadted fnodes: type = str
                       format = family | protein id | seed
                       description =  updated fnodes with only the family of the seed
  
· output files:
    - final table melt: type = str
                        format = genome_id | genome_name | seed | PA | color | protein_id
                        description = presence/abscence table, with header. Each line is a protein information
    - final table crosstab: type = str
                            format = row = genome_id, columns = seed 1 | seed 2 .. seed x
                            description = presence/abscence table, with header. Each line is a genome, each column is a seed.                   
```

#### Rule 9: plots

Rule description: Plot the table from `make_PA_table` in a colored table.

```
· input file:
    - final table melt : type = str
                         format = genome_id | genome_name | seed | PA | color | protein_id
                         description = presence/abscence table, with header. Each line is a protein information
  
· output files:
    - pdf: type = list of str
           description = plots in pdf of the final table centered on one seed 
    - png: type = list of str
           description = plots in png of the final table centered on one seed
```

### Additional rules

#### Rule A1: report_threshold

Rule description: Pair of rules that works together. From the blast output, create a report containing 2 scatter plots: one comparing the coverage and percentage identity of the hits and the other adding the evalue to the comparaison.

To run the rules do

```bash
snakemake report_threshold --cores 1 --use-conda
```

*Sub rule A1a: blast2threshold_table*

```
· input files:
    - seed file: type = str
                 format = seed name | protein id | hmm | e-value | percentage of identity | coverage | color
                 description: seed file input by user  
    - protein table: type = str
                     format =  protein id | protein name | genome name | genome id | length
                     description = final table of protein information.
    - updadted fnodes: type = str
                       format = family | protein id | seed
                       description =  updated fnodes with only the family of the seed
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description = output format of blast                                 
  
· output file:
    - table_thresholds : type = str
                         format = protein1 | protein 2 | pident | coverage | evalue | length
                         description = table with the columns from the blast output of only the pairs containing at least one member of the seed family.
```

*Sub rule A1a: report_threshold*

```
· input file:
    - tables_thresholds : type = list of str
                          format = protein1 | protein 2 | pident | coverage | evalue | length
                          description = table with the columns from the blast output of only the pairs containing at least one member of the seed family.                               
  
· output file:
    - report_html : type = str
                    description = report in html that contains the figure made using the coverage/percentage identity/evalue of all hits for each seeds family

```

#### Rule A2: extract_protein

Rule description: Create a fasta for each orthologs found (one fasta file per seed)

To run the rule do

```bash
snakemake extract_protein --use_conda --cores 1
```

or if on the cluster, do

```bash
snakemake extract_protein --use_conda --cores 1 --profile config/slurm
```

```
· input files:
    - final table melt : type = str
                         format = genome_id | genome_name | seed | PA | color | protein_id
                         description = presence/abscence table, with header. Each line is a protein information
    - prot sequence: type = str
                     description =  multifasta file of all the unique protein ids                         
  
· output file:
    - fasta files : type = list of str
                    description = name of all the fasta output with the format [name of the seed].fasta
```

#### Rule A3: quick_plots

Rule description: Plot a presence absence table that you may have edited in pdf and png.

To run the rule do

```bash
snakemake quick_plots -C PAtab_table=path/to/pa_table.tsv --use_conda --cores 1
```
or if on the cluster, do

```bash
snakemake quick_plots -C PAtab_table=path/to/pa_table.tsv --use_conda --cores 1 --profile config/slurm
```

Example of the `presence/absence table` in the [doc](https://github.com/vdclab/sORTholog/blob/main/doc/dummy_PAtab.tsv) folder in the GitHub page

```
· input files:
    - PA table : type = str
                 format = genome_id in the first column and the next columns are annotation columns
                 description = presence/abscence table, with header.
  
· output files:
    - pdf: type = list of str
           description = plots in pdf of the final table centered on one seed 
    - png: type = list of str
           description = plots in png of the final table centered on one seed
```

#### Rule A4: clean

Rule description: Remove the folder database, logs and processing_files to only keep results

To run the rule do

```bash
snakemake clean --cores 1
```

```
· input files:
    - database, logs and processing_files directory
```

[▲ Back to table of content ▲](#table-of-content)
