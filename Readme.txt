Presence Absence 1.0


#######
# Aim #
#######

This software will produce a table visualization of the presence and absence of encoded proteins in a given set of genomes. Both, the proteins and set of genomes are required as input by the user.

The primary output will be the table of presence-absence protein data filled using the NCBI Protein Database, as well as PDF and PNG files of the table visualization.


################
# Requirements #
################

To run this software you will need 4 files: the Snakefile, the user's seed file (format?), the user's file containing a list of NCBI Genome IDs, and a cluster.json file for appropriate connection to the appropriate dedicated clusters for running of the program.

You will also need these programs to be installed in the server (already done) or locally, on your computer:
- Snakemake 5.17.0
- python 2.8
- biopython 1.78	(python package)
- pandas 1.1.5		(python package)
- matplotlib 3.3.3	(python package)
- BLAST 2.10.1
- SILIX 1.2.11

the 4 files required should be as followed:

* Snakefile:
	do not edit it
	
* User seed file:
	User may select any name, but the file must be a tabular text table of either the .txt or .csv file format.
	This table should contain 3 columns, without headers: 
		chosen name for your protein of interest | ncbi protein ids | color in hexadecimal for image file (#000000 for black)
		
* User list of NCBI Genome IDs:
	User may select any name, but the file must be a tabular text table of either the .txt or .csv file format.
	This table should contain 1 column containing all NCBI Genome IDs; columns should be without names/headers.

* cluster.json file:
	This file contains the following:
		{
			"__default__":{
				"c" : "1",
				"qos" : "lagard",
				"time" : "72:00:00", 
				"account" : "lagard", 
				"mail-type" : "END",
				"mail-user": "(1)",
				"mem": "15gb"
			},
			"blast":{
				"c" : "5"
			}
		}
		
	The user will find this file ready to be edited where the user may enter their own e-mail in the place of the number 1, as it is shown in line 3 of the script.
	"time"= maximum amount of time allowed on the server
	"mem"= maximum memory allocated fo your computation
	Both variables, "time" and "mem", may also be edited by the user.


#######
# Use #
#######

Before starting the program, you will need to load other required software; this required software is listed as follows:
	module load snakemake/5.17.0
	module load python/3.8
	
The first will load the pipeline software, snakemake, that is required to read the Snakefile.
The second will load python v3.8 (release), which is the programmatic interpreter for snakemake.

To start the program, provide the following in the command line:
	python3.8 -m (path-1)/snakemake --cluster-config (path-2)/cluster.json --cluster "sbatch -c {cluster.c} --qos={cluster.qos} --time={cluster.time} --account={cluster.account} --mail-type={cluster.mail-type} --mail-user={cluster.mail-user} --mem={cluster.mem}"  -j 5 -d (path-3) -C project_name=(1) ncbi_id_list=(2) gene_tab=(3) (4)

	The following variables of the above script may be edited:
	(path-1)	path where user's snakefile is located. Ideally, this can be skipped if your are working within the same address that the file is stored (using the command cd)
	(path-2)	path where user's cluster.json file is located. Ideally, this file is located at the same place as your Snakefile, or with the user's Genome ID list and Seed files.
				If the cluster.json file is present in your current location, leave this blank.
	(path-3)	path where user's protein and genome data are. Always specify the full path to avoid problems.
	(1)			name of user's project provided in single or double quotes. Avoid spaces, replace any spaces with an underscore.
	(2)			name of user's file containing your list of NCBI Genome IDs, with the extension .txt or .csv (e.g. "file_name.txt" or 'file_name.txt').
	(3)			name of user's file containing the ID of user's protein of interest (seed file), with the extension .txt or .csv, in between single or double quotes.
	(4)			In addition, some optional inputs may be added and edited by the user:
					eval=(float)	Used to change the e-value threshold for BLAST. Can be on the 0.000001 (decimal) format or 10**-6 (exponential) format.
						DEFAULT = 1E-6
					id=(int)		Used to change the percentage of identity over the alignment threshold for SILIX. It should be an integer between 1 and 100.
						DEFAULT = 35
					cov=(int)		Used to change the coverage of the alignment threshold for SILIX. It should be an integer between 1 and 100.
						DEFAULT = 80
	
	*!* It is recommended that the user copy and paste this complete line and description of its components into user's notebook, changing the inputs as desired separately prior to copy-pasting into the command line for the relevant server.


	*******************
	* Important Notes *
	*******************
	
	1) To work on the HiPerGator server, user must first open a terminal (for unix) or command prompt (for windows).
	   User must access the server by entering the following command: ssh hpg2.rc.ufl.edu -l (GatorLink ID)
	   User will be required to enter their GatorLink account password. No double tap verification needed.
	   
	2) User will require access to files on the server. To access these files via command line, use the following: cd /blue/lagard
	   From here, if it has not already be done, user may create a working file to be associated with user using: mkdir (dir_name)
	   Then you can access it by using: cd /blue/lagard/(dir_name)
	   You can access these file on your computer by adding a server. The adresse is as follow: \\exasmb.rc.ufl.edu
	   For more info : *insert link*
	   
	3) The efficacy of snakemake is dependent upon pre-existing output files, but, without supply of these files by the user, the rule will not be triggered.
	   To trigger a rule again you can delete the output files of the rule, or change there names.
	   Alternatively, you can force a rule by adding it to the end of the user-entered command (in command line): --force (rule_name)

#################################
# Walk-Through and File Production #
#################################

This pipeline consists of 6 steps called rules that take input files and create output files. Here is a description of the pipeline.
*!* As snakemake is set up, there is a 7th rule, called all, that serves to call the last output files and make sure they were created.

	*****************************
	* Rule 1 : sequence_fetcher *
	*****************************
	
	*input files:
		- ncbi_id_list: user's list of NCBI Genome IDs [input ncbi_id_list]
		- seed: user's table containing the protein of interest [input gene_tab]
	
	*description:
		This rule will download the GenBank using the user-provided Genome ID list, recording the genome ID and name, as well as the protein ID and sequence of all the CDS of each genome.
		Then if the user's protein of interest has been detected across any of the genomes called, they will be downloaded from NCBI and the protein sequence recorded.
		All sequences recorded will be printed in a single FASTA file.
		*!* this step's length will increase linearly with the number of genomes contained in user's provided Genome ID list.
	
	*output files:
		- genome_prot_table: 	type = tabulated table with headers
								name = (input project_name)_all_protein.csv 
								format = genome_id | genome name | protein_id | sequence
								note = the seed not present in ypur dataset will have their genome id and name set to "Seed", not to be changed.
		- fasta_file:	type = multifasta
						name = (input project_name)_all_protein.fasta
						format = 	>protein_id
									sequence
						
						
	******************
	* Rule 2 : blast *
	******************
	
	*input file: sequence_fetcher fasta_file
	
	*decsription:
		Blast proteins created in rule 1, all against all, with a e_value threshold of 1E-3.
		This step lenght is exponentially increasing with the number of proteins.
	
	*output file:	type = blast out 6 file, tabulated table format, no headers
					name = blastall_(input project_name).out
					format = qseqid | sseqid | pident | length | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore
	
	*miscelenious files created: all file requiered to be created for a blast, a.k.a. the database files:
				(input project_name)_all_protein.fasta.pdb
				(input project_name)_all_protein.fasta.phr
				(input project_name)_all_protein.fasta.pin
				(input project_name)_all_protein.fasta.pot
				(input project_name)_all_protein.fasta.psq
				(input project_name)_all_protein.fasta.ptf
				(input project_name)_all_protein.fasta.pto
	
	
	*************************
	* Rule 3 : filter_blast *
	*************************
	
	*input file: blast output
	
	*input parameter: input eval
	
	*decsription:
		Filter the blast out file to remove all match higher than the e-value threshold (1E-6 if unchanged).
		Create a new folder for the rest of the analysis : (project_name)_eval(e_val)_id(id)_cov(cov), rfered as (new_dir) alter in this Readme.
	
	*output files:	type = blast out 6 file, tabulated table format, with headers
					name = (new_dir[see description])/blastall_(input project_name)_eval(input eval).out
	
	******************
	* Rule 4 : silix *
	******************
	
	*input files:
		- fasta: sequence_fetcher fasta_file output
		- blast_out: filter_blast output
		
	*input parameters:
		- id: input id
		- cov: input cov
		- eval: input eval
		- new_dir: new_dir defined in filter_blast
	
	*decsription:
		Uses graph theory to cluster proteins based on the input covergage and percentage of identity of the alignments.
		Each cluster is attributed a family number, and this info is recorded in a fnode file fo each protein.
		For each family a multifasta will be created with the family number containing the sequence of all the proteins of this family.
	
	*output files:	type = fnodes, tabulated table, no headers
					name = (new_dir)/(project_name)_eval(e_val)_id(id)_cov(cov).fnodes
					format = family | protein_id
	
	*miscelenious files created: one multifasta (name = cluster(family number)) per protein family in the directory (new_dir)/family_eval(e_val)_id(id)_cov(cov)
	
	***********************
	* Rule 5 : make_table *
	***********************
	
	*input files:
		- fnodes: silix output
		- seeds: your table of protein of interest [input gene_tab]
		- genome_prot_table: sequence_fetcher genome_prot_table output
	
	*decsription:
		Locate in wich silix family the given seeds are, and gather these families, using the fnodes file and the seeds file.
		Identify what protein is linked to what genome, using the genome_prot_table file.
		Complete a table of presence/absence based on these results, with the seed sin column, and the genomes in lines. 
		If a seed is not detected, the column will not be created.
	
	*output files:
		- gene_table:	type = tabulated table, with headers
						name = (new_dir)/gene_table_(project_name)_eval(e_val)_id(id)_cov(cov).csv
						format = genome_id | genome_name | seed_name_1 | seed_name_2 | ... | seed_name_n
						note = the seed columns are filled with ncbi protein ids
		- seed_table:	type = tabulated table, with headers
						name = (new_dir)/familied_seeds_table_(project_name)_eval(e_val)_id(id)_cov(cov).csv
						format = seed_name | seed_ncbi_id | color | silix_family
	
	*****************
	* Rule 6 : plot *
	*****************
	
	*input files:
		- gene_table: make_table gene_table output
		- seed_table: make_table seed_table output
	
	*decsription:
		This rule will make a graphic representation of the gene_table input, after removing the genome_id column.
		The graphic will be colored squared proteins are present, white square if absent.
		The color will be assigned by the seed_table input.
	
	*output files:
		- pdf:	type = pdf file
				name = (new_dir)/gene_table_(project_name)_eval(e_val)_id(id)_cov(cov).pdf
				format = graphic image of the input gene_table without the genome_id column
		- png:	type = png file
				name = (new_dir)/gene_table_(project_name)_eval(e_val)_id(id)_cov(cov).png
				format = graphic image of the input gene_table without the genome_id column
