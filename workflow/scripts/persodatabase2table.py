from Bio import SeqIO
import sys, os

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

parser_fasta = SeqIO.parse(snakemake.input.fasta_protein, 'fasta')

annotations_dict = {}

if annotations:
    annotations_dict = pd.read_table(annotations).set_index('protein_id').genome_id.to_dict()

with open(snakemake.output.tsv, 'wt') as w_file:
    header = ['protein_id', 'genome_id', 'length', 'protein_description']
    header = "\t".join(header)
    w_file.write(f'{header}\n')

    for prot in parser_fasta:
        length = len(prot.seq)

        # Security to remove protein id from description
        description = prot.description.replace(f"{prot.id} ", "")
        description = prot.description.replace(f"{prot.id}", "")

        if annotations_dict:
            genome_id = annotations_dict[prot.id]
            protein_id = f"{prot.id}--{genome_id}"

        else:
            id_split = prot.id.split('--')
            genome_id = id_split[1]
            protein_id = prot.id
        
        w_file.write(f"{protein_id}\t{genome_id}\t{length}\t{description}\n")
