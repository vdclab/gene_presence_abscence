# Presence-Abscence2.0

## Aim

This software will produce a table visualization of the presence and absence of encoded proteins in a given set of genomes. Both, the proteins and set of genomes are required as input by the user.

The primary output will be the table of presence-absence protein data filled using the NCBI Protein Database, as well as PDF and PNG files of the table visualization.

## Requirements 

To run this software you will need 4 files: the Snakefile (conatining the code, do not edit), the user's seed file (format explained bellow), the user's file containing a list of NCBI Genome IDs (format explained bellow), and a cluster.json file (format explained bellow, to be edited) for appropriate connection to the appropriate dedicated clusters for running of the program.

You will also need these programs to be installed in the server (already done) or locally, on your computer:

### Software 
- python == 3.8
- Snakemake == 6.4.1
- BLAST == 2.10.1
- SILIX == 1.2.11

### Python libraries
- biopython == 1.78
- pandas == 1.1.5	
- matplotlib == 3.3.3
- ete3 == 3.1.2
- ncbi-genome-download == 0.3.0

### Requiered files 

The 4 files required should be as followed:

```
- Snakefile:
  code file
	do not edit it
	
- User seed file:
  By default, the pipeline will recognize this file with the name 'seeds.txt'. This may be changed.
	User may select any name, but the file must be a tabular text table of either the .txt or .csv file format.
	This table should contain 3 columns, without headers in order: 
	  * chosen name for your protein of interest
    * ncbi protein id
    * BLAST e-value threshold for this protein
    * Minimum percentage of identity threshold for this protein
    * Minimum coverage thershold for this protein
    * color in hexadecimal for image file (#000000 for black)
		
- User list of NCBI Genome IDs:
	User may select any name, but the file must be a tabular text table of either the .txt or .csv file format.
	This table should contain 1 column containing all NCBI Genome IDs; columns should be without names/headers.

- cluster.json file:
	This file contains the following:
		{
			"__default__":{
				"c" : "1",
				"qos" : "lagard",
				"time" : "D-HH:MM:SS", 
				"account" : "lagard", 
				"mail-type" : "END",
				"mail-user": "(1)",
				"mem-per-cpu": "15gb"
			},
			"psiblast":{
		    "c" : "5"
	    },
	    "fetch_proteins_info_and_blastall":{
		    "c" : "5"
			}
		}
		
	The user will find this file ready to be edited where the user may enter their own e-mail in the place of the number 1, as it is shown in line 3 of the script.
	"time"= maximum amount of time allowed on the server
	"mem"= maximum memory allocated fo your computation
	Both variables, "time" and "mem", may also be edited by the user.
```

These 4 files should be in the same folder, hereafter your woking folder

## Use 

Before starting the program, you will need to load other required software; this required software is listed as follows:

```
module load python/3.8 snakemake/6.4.1
```

The first will load the pipeline software, snakemake, that is required to read the Snakefile.
The second will load python3.8 (release), which is the programmatic interpreter for snakemake.

Before running the program, you need to be in your working folder. To do so use this command line:

```
cd /blue/lagard/(GatorLink ID)/(woking folder)
```

To start the program, provide the following in the command line:

```
python3.8 -m snakemake --cluster-config cluster.json --cluster "sbatch -c {cluster.c} --qos={cluster.qos} --time={cluster.time} --account={cluster.account} --mail-type={cluster.mail-type} --mail-user={cluster.mail-user} --mem-per-cpu={cluster.mem-per-cpu}"  -j 5 -d (working path) -C taxid=(1) (2)
```
If you want to prevent the software from sending you e-mail, you can remove the --mail-user option

The following variables of the above script may be edited

```			
(working path)	path where user's protein and genome data are. Always specify the full path to avoid problems.
	(1)	name of user's file containing your list of NCBI TaxID, with the extension .txt or .csv (e.g. "file_name.txt" or 'file_name.txt').
  	(2)	In addition, some optional inputs may be added and edited by the user:
     		 project_name=(str) name of user's project provided in single or double quotes. Avoid spaces, replace any spaces with an underscore.
          		DEFAULT = seed_file name
      		seed=(str) name of user's file containing the ID of user's protein of interest (seed file), with the extension .txt or .csv, in between single or double quotes.
          		DEFAULT = 'seeds.txt'
		eval=(float)	Used to change the e-value threshold for PSI-BLAST. Can be on the 0.000001 (decimal) format or 10**-6 (exponential) format.
			DEFAULT = 1E-6
      		blast=(str) BLAST version to use.
         		 DEFAULT = '2.10.1'
      		silix=(str) Silix version to use.
          		DEFAULT = '1.2.11'
```

**/!\ It is recommended that the user copy and paste this complete line and description of its components into user's notebook, changing the inputs as desired separately prior to copy-pasting into the command line for the relevant server.**

### Important Notes 
	
1) To work on the HiPerGator server, user must first open a terminal (for unix) or command prompt (for windows).  
User must access the server by entering the following command:  
```
ssh (GatorLink ID)@hpg.rc.ufl.edu
```
User will be required to enter their GatorLink account password. **No double tap verification needed.**  
	  
2) User will require access to files on the server. To access these files via command line, use the following:  
```
cd /blue/lagard/(GatorLink ID)
```
From here, if it has not already be done, user may create a working file to be associated with user using: mkdir (dir_name)  
Then you can access it by using:  
```
cd (dir_name)
```
You can access these file on your computer by adding a server. The adresse is as follow: `\\exasmb.rc.ufl.edu`  
For more info : https://help.rc.ufl.edu/doc/Samba_Access
	   
3) The efficacy of snakemake is dependent upon pre-existing output files, but, without supply of these files by the user, the rule will be triggered.

To trigger a rule again you can delete the output files of the rule, or change there names.

Alternatively, you can force a rule by adding it to the end of the user-entered command (in command line): `--force (rule_name)` or `-f (rule_name)`.


## Walk-Through and File Production

This pipeline consists of 12 steps called rules that take input files and create output files. Here is a description of the pipeline.

1. As snakemake is set up, there is a 13th rule, called `all`, that serves to call the last output files and make sure they were created.

2. A folder containing your work will be created with two folder inside, if not already existing: `processing_file` and `results`. You will find any result files in the `results` folder.

3. Before starting the pipeline, your taxids will be checked and updated for the lowest levels. This will create a file of checked taxids. It will skip this step if this file already exists.

4. When restarting the pipeline, the software will check if you made any changes in the seed file before running. If changes have been made, it will run what is necessary, else nothing will happen.

### Pipeline in image 

#### Normal behavior

![](doc/dummy_dag.png "Normal behavior" =250x250)

#### Speedup behavior

![](doc/dummy_dag_speedup.png "Speedup behavior" =250x250)

### Rule 1 : fetch_fasta_from_seed

Rule description:  
- Fetch the fasta of the seed from the seed table. Then they are writen in the output file.  

```
· input file: 
    - seed file: type = str
                              format = name | protein id | e-value | percentage of identity | coverage | color
                              note = no header on this file.

· output file: 
    - multi fasta file: type = str
	                    description = multifasta output of the seed sequences 
```        
        
       
### Rule 2 : psiblast

Rule description: Use the sequences of the seeds to make a psiBLAST against all the taxid. 
 
```
· input files:
    - seed: type = str
            description = multifasta of the seed sequences
    - taxid: type = str 
             description = list of taxid in columns, no header
        
· output file (aka blast out):  
    - blast out: type = str
                 format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                 description = blast out format in tabulation, no header
```        
        
        
### Rule 3 : read_psiblast

Rule description: Read the psiBLAST, remove unwanted lines ane extract the list of matches. 

```
· input file : 
    - psiblast output: type = str
                       format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                       description = blast out format in tabulation, no header
    
· params:
    - e-val: type = int
             description =  e-value threshold chosen by user
    - blast version: type = str
                     description = blast version to use  
        
        
· output files:
    - clean blast: type = str
                   format = query accession | query length | query sequence | query start position | querry end position | subject accession | subject length | subject sequence| subject start position | subject end position | length of alignment | percentage of identity | e-value | bitscore | querry coverage
                   description = cleaned blast out format in tabulation, no header
    - list all prot: type = str
                     description = list of all potein identifications gathered in the psiBLAST in column
```

        
        
### Rule 4 : split_taxid_file

Rule description: Split the taxid file into 3. 

```

· input file (aka ): 
    - taxid file: type = str
                  description = verified taxid file, list of taxid in column, no header
        
· params:
    - nb_per_file: type = int
                   description = number of taxid per file, detrmined by the fetch_splitter function.

· output file:
    - taxid files: type = list of str
                   description = three files of list of taxid in column, no header, determined by the fetch_splitter function.
```       
        
        
### Rule 5 : fetch_proteins_info

Rule description: Fetch the information for each protein of each genome in the taxid list. That includes: the protein ncbi id, sequence, length and annotation, as well as in which genome is found. Information for the genome include genome ncbi id, name, taxid and if complete or partial.  
 
```

· input file: 
    - taxid file: type = str
                  description = list of taxid in column, no header, from the rule plit_taxid_file (one file)

  
· output file: 
    - table of protein informations: type = str
                                     format = protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
                                     description = table of the information collected on the proteins, without header.
  
```
  
  
### Rule 6 : cat_proteins_info

Rule description: Concatenate the different table of protein info created in the rule fetch_proteins_info. Then remove all file created in the rules split_taxid_file and fetch_proteins_info. 

```

· input file: 
    - all protein information files: type = list of str
                                     format = protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
                                     description = list of tables of the information collected on the proteins, without header.
        
· output file: 
    - concatenated all protein information file: type = str
                                                 format = protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
                                                 description = final table of protein information, without header.

```


### Rule 7 : make_fasta

Rule description: Create a fasta file from the psiblast results and the result of the protein information in the rule cat_proteins_info. 

```

· input files:
    - protein table: type = str
                     format = protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
                     description = final table of protein information from the rule cat_proteins_info, without header.
    - list all prot: type = str
                     description = list of all protein identifications gathered in the psiBLAST in column
  
· output file: 
    - multifasta: type = str
                  description = multifasta file of all the unique protein ids.
    - reduced protein table: type = str
                             format = protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
                             description = final table of protein information with removed duplicates, without header.

```


### Rule 8 : blast

Rule description: Blast all versus all of the fasta of all protein generated in the rule make_fasta. 

```

· input files:
    - prot sequence: type = str
                     description =  multifasta file of all the unique protein ids from the rule make_fasta
    - seed fasta: type = str
                  description = multifasta file of all the seeds from the rule fetch_fasta_from_seed
     
· params:
    - blast version: type = str
                     description = blast version to use 
  
· output files:
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description = output format of blast
    - fasta for blast: type = str
                       description = concatenation of the 2 input multifasta files
```



### Rule 9 : prepare_for_silix

Rule description: Filter the blast results from the rule blast with the threshold specified for each seed in the seed file. Filters include the identity score, coverage and e-value, decided by the user. Create one new filtered blast result for each seed.  

```

· input file: 
    - fasta: type = str
             description = multifasta of proteins with seed from the rule blast
    - blast out: type = str
                 format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                 description =  blast output from the rule blast  
  
· params:
    - new dir: type = str
               description = work directory
    - project name: type = str
                    description = project name determined by the user
  
  - output file: 
    - blast output to send to silix: type = list of str
                                     format = query id | subject id | percentage of identity | length of match  | mismatch | gapopen | query start position | query end position | subject start position | subject end position | e-value | bitscore
                                     description = list of blast output filtered for each seed.

```


### Rule 10 : silix

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


### Rule 11 : find_family
 
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


### Rule 12 : make_table

Rule description: Check the presence of protein similar to the seed in each taxid and create a table of presence abscence. This table is then plotted in a colored table.

```

· input files:
    - protein table: type = str
                     format = protein id | protein name | genome name | genome status | genome id | taxid | length | sequence
                     description = final table of protein information from the rule cat_proteins_info, without header.
    - fnodes: type = str
              format = family | protein id | seed
              description =  concatenated fnodes with each seed family
  
· output files:
    - final table: type = str
                   format = genome id | genome name | seed 1 | seed 2 .. seed x
                   description = presence/abscence table, with header. Each line is a genome, each column is a seed.
    - pdf: type = list of str
           description = plots in pdf of the final table centered on one seed 
    - png: type = list of str
           description = plots in png of the final table centered on one seed
```
