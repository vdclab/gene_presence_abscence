# Snakemake workflow: sORTholog

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.14.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/vdclab/sORTholog/workflows/Tests/badge.svg?branch=main)](https://github.com/vdclab/sORTholog/actions?query=branch%3Amain+workflow%3ATests)
[![License (AGPL version 3)](https://img.shields.io/badge/license-GNU%20AGPL%20version%203-green.svg)](COPYING)

## Aim

This software will produce a table visualization of the presence and absence of encoded proteins in a given set of genomes. Both, the proteins and set of genomes are required as input by the user.

The primary output will be the table of sORTholog protein data filled using the NCBI Protein Database, as well as PDF and PNG files of the table visualization.

## Installation

### Step 1: install Snakemake and Snakedeploy

Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).

To install Mamba by conda, run

```shell
conda install mamba -n base -c conda-forge
```

Given that Mamba is installed, run 

```shell
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
```

to install both Snakemake and Snakedeploy in an isolated environment. 

#### Notes 

For all following commands (step 2 to 5) ensure that this environment is activated via 

```shell
conda activate snakemake
```

### Step 2: deploy workflow

 Given that Snakemake and Snakedeploy are installed and available (see Step 1), the workflow can be deployed as follows.

First, create an appropriate project working directory on your system in the place of your choice as follow (note to change the path and file name to the one you want to create): : 

```shell
mkdir path/to/project-workdir
```

Then go to your project working directory as follow:

```shell
cd path/to/project-workdir
```

In all following steps, we will assume that you are inside of that directory.

Second, run 

```shell
snakedeploy deploy-workflow https://github.com/vdclab/sORTholog . --tag 0.1.2
```

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder.

### Step 3: configure workflow

#### General settings

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.  

Seed and taxonomy sheet
- Add seeds to `config/samples.tsv`. For each protein, the columns `seed`, and `protein_id` have to be defined. The protein_id is the protein id define by [NCBI](https://www.ncbi.nlm.nih.gov/). To include other relevant variables such as color, threshold (e-value, percentage identities and coverage), add a new column to the sheet. Exemple of the `seed file` is present in the `config` folder if needed or in the [doc](https://github.com/vdclab/sORTholog/blob/main/doc/dummy_seeds.tsv) folder in the GitHub page
- Add taxonomic id list to `config/taxid.tsv`. For each genome, the column `TaxId` has to be defined. The TaxId is define by [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy). To include other relevant variables such as taxonimical groups (`NCBIGroups`), add a new column to the sheet. The groups that could be used are : `all`, `archaea`, `bacteria`, `fungi`, `invertebrate`, `metagenomes`, `plant`, `protozoa`, `vertebrate_mammalian`, `vertebrate_other` or `viral`.
- You can also configure other parameter such as `project_name`, `output_folder` in the `config/config.yaml` file. Exemple of the `tanonomy file` is present in the `config` folder if needed or in the [doc](https://github.com/vdclab/sORTholog/blob/main/doc/dummy_taxids.tsv) folder in the GitHub page

Missing values can be specified by empty columns or by writing `NA`.

#### Useful options

##### Choose the border shade of the final figure

The shade of the border can be changed using the `colored_border` options in the `config/config.yaml` file:
- if set to `True` the shade of the border is infered from the color set for the seed in the `seed file`
- if set to `False` the color will be black (Default setting)

Here a comparison of the two behaviours:

<p align="center">
  <img src="doc/colored_border_option.png?raw=true">
</p>

##### Choose the border shape of the final figure

The shape of the border can be changed using the `round_border` option in the `config/config.yaml` file:
- if set to `True` the shape of the border will be round
- if set to `False` the shape of the border will be straight (Default setting)

Here a comparison of the two behaviours:

<p align="center">
  <img src="doc/round_border_option.png?raw=true">
</p>

### Step 4: run the workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda, run Snakemake with 

```shell
snakemake --cores 1 --use-conda 
```

If you want to run the workflow with the additional step to speedup the analysis that contains a big dataset you can either:
- Change in the `config.yaml` the parameter `speedup` and change the value to `True`
- Run the workflow adding `-C speedup=True` to the command line as follow

```shell
snakemake --cores 1 --use-conda -C speedup=True
```

### Step 5: generate report

After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser via 

```shell
snakemake --report report.zip --report-stylesheet config/report.css
```
The resulting report.zip file can be passed on to collaborators, provided as a supplementary file in publications.

### Side notes 

The efficacy of snakemake is dependent upon pre-existing output files, but, without supply of these files by the user, the rule will be triggered.

If you want to re-run a part of the pipeline without running previous rule you can specify the rule you want to run with the option `allowed-rules`

```bash
  --allowed-rules ALLOWED_RULES [ALLOWED_RULES ...]
                        Only consider given rules. If omitted, all rules in Snakefile are used. Note that this is intended primarily for internal use and may lead to unexpected results otherwise. (default: None)
```

For exemple if I want to re-run from silix to all

```bash
snakemake --cores 1 --use-conda --allowed-rules silix find_family make_table plots all
```

## Important Notes for running on the cluster

### Step C1: log on the cluster [UFL users]

1) To work on the HiPerGator server, user must first open a terminal (for unix) or command prompt (for windows).  
User must access the server by entering the following command:  
```
ssh (GatorLink ID)@hpg.rc.ufl.edu
```
User will be required to enter their GatorLink account password.
      
2) User will require access to files on the server. To access these files via command line, use the following:  
```
cd /blue/(your group)/(GatorLink ID)
```

For exemple `group` is `lagard` for the VDCLab

You can access to the files of this folder on your computer by adding a server. The address for Windows is as follow: `\\exasmb.rc.ufl.edu`, for Unix (Linux or MacOS) the address: `smb://exasmb.rc.ufl.edu`.  
For more info : https://help.rc.ufl.edu/doc/Samba_Access

### Step C2: load conda

Before starting the installion of the program, you will need to load `conda`; this required software is listed as follows:

```
module load conda
```

After that the step are the same as in `Step 1: install Snakemake and Snakedeploy`, `Step 2: deploy workflow`, `Step 3: configure workflow`.

### Step C3: run the workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda and running on the cluster, run Snakemake with 

```shell
snakemake --cores 5 --use-conda --profile config/slurm
```

As previously, if you want to run the workflow with the additional step to speedup the analysis that contains a big dataset you can either:
- Change in the `config.yaml` the parameter `speedup` and change the value to `True`
- Run the workflow adding `-C speedup=True` to the command line as follow

```shell
snakemake --cores 1 --use-conda --profile config/slurm -C speedup=True 
```

As previously you can also generate a report following `Step 5: Generate report`.

### Additionnal information

If you want to run it as a job,

1. make sure that the pipeline is deploy where you are and that you change the `config/config.yaml` accordingly to your needs
2. create a file name for exemple `my_project.sh` and copy paste the following

```bash
#!/bin/bash -login
#SBATCH -J sbatch_snakemake 
#SBATCH -t 7-0:00:00
#SBATCH -N 1
#SBATCH --output $PWD/sbatch_snakemake-%j.out
#SBATCH --error $PWD/sbatch-snakemake-%j.err

# activate conda in general
source /home/$(whoami)/.bashrc # if you have the conda init setting

# activate a specific conda environment, if you so choose
conda activate snakemake 

# make things fail on errors
set -o nounset
set -o errexit
set -x

### run your commands here!
snakemake --cores 5 --use-conda --profile config/slurm
snakemake --report report.zip --report-stylesheet config/report.css
```

3. run the following command

```bash
sbatch my_project.sh
```

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

### Pipeline in image 

#### Normal behavior

<p align="center">
  <img src="doc/dummy_dag.png?raw=true" height="400">
</p>

#### Speedup behavior

<p align="center">
  <img src="doc/dummy_dag_speedup.png?raw=true" height="500">
</p>

### Rule 1: fetch_fasta_from_seed

Rule description: Fetch the fasta of the seed from the seed table. Then they are writen in the output file.

```
· params file: 
    - seed file: type = str
                 columns = seed name | protein id | e-value | percentage of identity | coverage | color

· output file: 
    - multi fasta file: type = str
	                    description = multifasta output of the seed sequences 
    - new seed file: type = str
                     columns = seed name | protein id | e-value | percentage of identity | coverage | color      
                     description = Same file as input but the 'protein id' of each file is changed to match the fasta file                
```
        
### Rule 2: fetch_proteins_database

Rule description: Fetch the information for each protein of each genome in the taxid list. That includes: the protein ncbi id, sequence, length and annotation, as well as in which genome is found. Information for the genome include genome ncbi id, name, taxid and if complete or partial.
 
```

· input file: 
    - taxid file: type = str
                  description = list of taxid in column

  
· output file: 
    - multi fasta file: type = str
                        description = multifasta output of the genomes proteomes in one file 
    - table of genomes informations: type = str
                                     format = assembly_accession | bioproject | biosample | wgs_master | excluded_from_refseq | refseq_category | relation_to_type_material | taxid | species_taxid | organism_name | infraspecific_name | isolate | version_status | assembly_level | release_type | genome_rep | seq_rel_date | asm_name | submitter | gbrs_paired_asm | paired_asm_comp | ftp_path
                                     description = table of the information collected from the NCBI Assembly.
    - table of protein informations: type = str
                                     format = protein id | protein name | genome name | genome id | length
                                     description = table of the information collected on the proteins.
    - new taxid file: type = str
                      description = list of taxid downloaded in column

```
 
### Rule 2.1: psiblast (optional)

Rule description: Use the sequences of the seeds to make a psiBLAST against all the taxid. 
 
```
· input files:
    - seed fasta: type = str
                  description = multifasta of the seeds sequence
    - taxid fasta: type = str 
                   description = multifasta of the genomes sequences
        
· output file (aka blast out):  
    - blast out: type = str
                 format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                 description = blast out format 6 in tabulation, no header
```

### Rule 2.2: read_psiblast (optional)

Rule description: Read the psiBLAST, remove unwanted lines ane extract the list of matches. 

```
· input file : 
    - psiblast output: type = str
                       format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                       description = blast out format in tabulation, no header
               
· output files:
    - clean blast: type = str
                   format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                   description = cleaned blast out format
    - list all prot: type = str
                     description = list of all potein identifications gathered in the psiBLAST in column
```

### Rule 2.3: make_fasta (optional)

Rule description: Create a fasta file from the psiblast results and the result of the protein information in the rule cat_proteins_info. 

```

· input files:
    - protein table: type = str
                     format = protein id | protein name | genome name | genome id | length
                     description = table of protein information 
    - list all prot: type = str
                     description = list of all protein identifications gathered in the psiBLAST in column
  
· output file: 
    - multifasta: type = str
                  description = multifasta file of all the unique protein ids.

```

### Rule 3: merge_fasta

Rule description: Merge the proteome of interest with the seeds 

```

· input files:
    - prot sequence: type = str
                     description =  multifasta file of all the unique protein ids
    - seed fasta: type = str
                  description = multifasta file of all the seeds from the rule fetch_fasta_from_seed
  
· output files:
    - fasta for blast: type = str
                       description = concatenation of the 2 input multifasta files
```



### Rule 4: blast

Rule description: Blast all versus all of the fasta of all proteins. 

```

· input files:
    - fasta for blast: type = str
                       description = concatenation of the the proteome taxid and seeds multifasta files
  
· output files:
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description = output format of blast
```



### Rule 5: prepare_for_silix

Rule description: Filter the blast results from the rule blast with the threshold specified for each seed in the seed file. Filters include the identity score, coverage and e-value, decided by the user. Create one new filtered blast result for each seed.

```

· input file: 
    - protein table: type = str
                     format = protein id | protein name | genome name | genome id | length
                     description = table of protein information 
    - seed file: type = str
                 columns = seed name | protein id | e-value | percentage of identity | coverage | color             
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description =  blast output from the rule blast  
  
  - output file: 
    - blast output to send to silix: type = list of str
                                     format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                                     description = list of blast output filtered for each seed.

```


### Rule 6: silix

Rule description: Uses Silix to create a network of protein and give a file of the protein segregated in groups.  If the blast output file is empty, just create an empty file. 

```

· input files:
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description = blast output filtered for a specific seed from the rule prepare_for_silix.
    - fasta: type = str
             description = multifasta of proteins with seed from the rule blast
     
· params:
    - silix version: type = str
                     description =  version of silix to use
  
· output file:
    - fnodes: type = str 
              format = family | protein id
              description = fnodes file, table of protein id and family number, without headers.

```


### Rule 7: find_family
 
Rule description: Find the group of each seed in each individual seed and record it. 

```

· input file:
    - fnodes: type = str
              format = family | protein id
              description = fnodes file, table of protein id and family number, without headers from the rule silix.
  
· output file:
    - updadted fnodes: type = str
                       format = family | protein id | seed
                       description = updated fnodes with only the family of the seed.

```


### Rule 8: make_table

Rule description: Check the presence of protein similar to the seed in each taxid and create a table of presence abscence. This table is then plotted in a colored table.

```

· input files:
    - seed file: type = str
                 columns = seed name | protein id | e-value | percentage of identity | coverage | color 
    - protein table: type = str
                     format =  protein id | protein name | genome name | genome id | length
                     description = final table of protein information.
    - updadted fnodes: type = str
              format = family | protein id | seed
              description =  updated fnodes with only the family of the seed
    - table of genomes informations: type = str
                                     format = assembly_accession | bioproject | biosample | wgs_master | excluded_from_refseq | refseq_category | relation_to_type_material | taxid | species_taxid | organism_name | infraspecific_name | isolate | version_status | assembly_level | release_type | genome_rep | seq_rel_date | asm_name | submitter | gbrs_paired_asm | paired_asm_comp | ftp_path
                                     description = table of the information collected from the NCBI Assembly.
  
· output files:
    - final table melt : type = str
                         format = genome_id | genome_name | seed | PA | color | protein_id
                         description = presence/abscence table, with header. Each line is a protein information
    - final table crosstab : type = str
                             format = row = genome_id, columns = seed 1 | seed 2 .. seed x
                             description = presence/abscence table, with header. Each line is a genome, each column is a seed.                   
```



### Rule 9: plots

Rule description: The table from `make_table` is then plotted in a colored table.

```

· input files:
    - final table melt : type = str
                         format = genome_id | genome_name | seed | PA | color | protein_id
                         description = presence/abscence table, with header. Each line is a protein information
  
· output files:
    - pdf: type = list of str
           description = plots in pdf of the final table centered on one seed 
    - png: type = list of str
           description = plots in png of the final table centered on one seed
```

### Additionnal rules

#### Rule A1: report_threshold

Rule description: Pair of rules that works together. From the blast output, create a report containing 2 scatter plots: one comparing the coverage and percentage identity of the hits and the other adding the evalue to the comparaison.

To run the rules do

```bash
snakemake report_threshold --cores 1 --use-conda
```

##### Sub rule A1a: blast2threshold_table
```

· input files:
    - seed file: type = str
                 columns = seed name | protein id | e-value | percentage of identity | coverage | color 
    - protein table: type = str
                     format =  protein id | protein name | genome name | genome id | length
                     description = final table of protein information.
    - updadted fnodes: type = str
              format = family | protein id | seed
              description =  updated fnodes with only the family of the seed
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description = output format of blast                                 
  
· output files:
    - table_thresholds : type = str
           format = protein1 | protein 2 | pident | coverage | evalue | length
           description = table with the columns from the blast output of only the pairs containing at least one member of the seed family.
```

##### Sub rule A1a: report_threshold
```

· input files:
    - tables_thresholds : type = list of str
           format = protein1 | protein 2 | pident | coverage | evalue | length
           description = table with the columns from the blast output of only the pairs containing at least one member of the seed family.                               
  
· output files:
    - report_html : type = str
           description = report in html that contains the figure made using the coverage/percentage identity/evalue of all hits for each seeds family
```

#### Rule A2: extract_protein

Rule description: Create a fasta for each orthologs found (one fasta file per seed)

To run the rule do

```bash
snakemake extract_protein --use_conda --cores 1
```

```

· input files:
    - final table melt : type = str
                         format = genome_id | genome_name | seed | PA | color | protein_id
                         description = presence/abscence table, with header. Each line is a protein information
    - prot sequence: type = str
                     description =  multifasta file of all the unique protein ids                         
  
· output files:
    - fasta files : type = list of str
           description = name of all the fasta output with the format [name of the seed].fasta
```

#### Rule A3: quick_plots

Rule description: Plot a presence absence table that you may have edited in pdf and png.

To run the rule do

```bash
snakemake quick_plots -C PAtab_table=path/to/pa_table.tsv --use_conda --cores 1
```

Exemple of the `presence/absence table` in the [doc](https://github.com/vdclab/sORTholog/blob/main/doc/dummy_PAtab.tsv) folder in the GitHub page

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
