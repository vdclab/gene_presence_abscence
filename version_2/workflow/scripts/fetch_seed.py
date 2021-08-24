from Bio import Entrez

Entrez.tool = 'draw presence/abscence v2'
Entrez.email = 'decrecylab@gmail.com'

# Read line by line the seed file without loading it in memory
with open(snakemake.input[0], 'rt') as r_file :
    header = r_file.readline().split() 
    index_proteinId = header.index('protein_id') 
    id_list = [] 
    
    for line in r_file : 
        tmp_line = line.split() 
        id_list.append(tmp_line[index_proteinId]) 

# getting seed sequences and writing the fasta file
with Entrez.efetch(db='protein', id=id_list, rettype='fasta', retmode='text') as handle:
    with open(snakemake.output.fasta_seed, 'w') as out_file:
        out_file.write(handle.read())

# Getting lenght of the seeds proteins and adding to new table seeds
seed_dict = SeqIO.index(snakemake.output.fasta_seed, index)

# Parsing of the table to add the length to the table
with open(snakemake.input[0], 'rt') as r_file :
    with open(snakemake.output.new_seed_file, 'wt') as wt
        header = r_file.readline().split() 
        index_proteinId = header.index('protein_id') 
        
        # Adding length to header
        new_line = f"{line.rstrip()}\tlength\n"
        w_file.write(new_line)

        for line in r_file : 
            tmp_line = line.split() 

            protein_id = [protein_id for protein_id in seed_dict.keys() 
                                if tmp_line[index_proteinId] in protein_id][0]

            length = len(seed_dict[protein_id])

            new_line = f"{line.rstrip()}\t{length}\n"
            w_file.write(new_line)



