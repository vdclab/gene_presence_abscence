from Bio import Entrez

Entrez.tool = 'draw presence/abscence v2'
Entrez.email = 'decrecylab@gmail.com'

# Read line by line the seed file without loading it in memory
with open(snakemake.input[0], 'rt') as r_file :
    header = r_file.readline().rstrip().split() 
    index_proteinId = header.index('protein_id') 
    id_list = [] 
    
    for line in r_file : 
        tmp_line = line.split() 
        id_list.append(tmp_line[index_proteinId]) 

# getting seed sequences and writing the fasta file
with Entrez.efetch(db='protein', id=id_list, rettype='fasta', retmode='text') as handle:
    with open(snakemake.output[0], 'w') as out_file:
        out_file.write(handle.read())