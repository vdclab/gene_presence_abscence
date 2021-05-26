# snakemake/5.17.0
# usage: snakemake --cluster-config cluster.json --cluster "sbatch -c {cluster.c} --qos={cluster.qos} --time={cluster.time} --account={cluster.account} --mail-type={cluster.mail-type} --mail-user={cluster.mail-user}"  -j 5 -d "/blue/lagard/ghutinet/fpreQ1_test" -C fasta_name="unculturephages_and_que.fasta" gene_tab="seeds.txt" id=0.1 seeds_in=False 
# -j : number of max core to use
# -d : directory of the files
# -C : configuration of values bellow:

##########################################################################
##########################################################################
##
##                                Options
##
##########################################################################
##########################################################################

#name your project
project_name = config['project_name']

#name of the polyfasta containing your proteins
ncbi_id_list = config['ncbi_id_list']

#name of the file containing the gene and color table
gene_tab = config['gene_tab']

#Blast e-value thershold, 0.000001 by default but can be changed in -C
e_val = config['e_val'] if 'e_val' in config else 0.000001

#Blast percentage of identity threshold, 35% by default but can be changed in -C
id = config['id'] if 'id' in config else 0.35

#Blast percentage of coverage threshold, 80% by default but can be changed in -C
cov = config['cov'] if 'cov' in config else 0.8

##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os

##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################


new_dir = f"{project_name}_eval{e_val}_id{id}_cov{cov}"

if not os.path.exists(new_dir):
    os.mkdir(new_dir)

##########################################################################
##########################################################################
##
##                                Rules
##
##########################################################################
##########################################################################

rule all:
    input:
         pdf=f"{new_dir}/gene_table_{project_name}_eval{e_val}_id{id}_cov{cov}.pdf",
         png=f"{new_dir}/gene_table_{project_name}_eval{e_val}_id{id}_cov{cov}.png"

##########################################################################
##########################################################################

rule sequence_fetcher:
    input:
        ncbi_id_list = ncbi_id_list,
        seeds = gene_tab
    output:
        genome_prot_table = f'{project_name}_all_protein.csv',
        fasta_file = f'{project_name}_all_protein.fasta'
    run:
        from Bio import Entrez
        import pandas as pd
        import ssl

        Entrez.tool = 'draw presence/abscence'


        def _error_get_record(list_id, number):
            ncbi_records = []
            len_list = len(list_id)

            for new_start in range(0, len_list, 10 ** number):
                new_end = new_start + 10 ** number if new_start + 10 ** number <= len_list else len_list
                tab = (5 - number) * "\t"
                print(f'{tab}{new_start} {new_end} {len_list}')
                ncbi_records += get_record(list_id=list_id[new_start:new_end], number=number)

            return ncbi_records


        def get_record(list_id, number):
            try:
                with Entrez.efetch(db='nucleotide', id=list_id, rettype='gbwithparts', retmode='xml') as ncbi_handle:
                    ncbi_records = Entrez.read(ncbi_handle)
                del ncbi_handle

            except ssl.SSLError:
                print(ssl.SSLError)
                ncbi_records = _error_get_record(list_id, number-1)

            except MemoryError:
                print(f'MemoryError, number ={number}')
                ncbi_records = _error_get_record(list_id, number-1)

            return ncbi_records


        def break_record(ncbi_records, final_table):
            for ncbi_record in ncbi_records:

                for feature in ncbi_record['GBSeq_feature-table']:

                    if feature['GBFeature_key'] == 'CDS':
                        feature_quals = {}

                        for qual in feature['GBFeature_quals']:
                            try:
                                feature_quals[qual['GBQualifier_name']] = qual['GBQualifier_value']
                            except KeyError:
                                continue

                        if 'protein_id' in feature_quals:
                            if not feature_quals['protein_id'] in final_table.prot_id.to_list():
                                extracted_data = {'genome_id': [ncbi_record['GBSeq_accession-version']],
                                                  'genome_name': [ncbi_record['GBSeq_organism']],
                                                  'prot_id': [feature_quals['protein_id']],
                                                  'sequence': []
                                                  }

                                try:
                                    extracted_data['sequence'] = [feature_quals['translation']]
                                except KeyError:
                                    with Entrez.efetch(db='protein',
                                                       id=feature_quals['protein_id'],
                                                       rettype='gb', retmode='xml'
                                                       ) as new_handle:
                                        extracted_data['sequence'] = [Entrez.read(new_handle)[0]['GBSeq_sequence']]
                                    del new_handle

                                final_table = final_table.append(pd.DataFrame(extracted_data), ignore_index=True)

            return final_table


        id_table = pd.read_csv(str(input.ncbi_id_list), header=None)
        id_list = []
        removed = []

        for id in id_table[0].to_list():
            try:
                Entrez.read(Entrez.esummary(db='nucleotide',id=id))
                id_list += [id]
            except:
                removed += [id]

        if len(removed) > 0:
            pd.DataFrame(removed).to_csv('ncbi_id_removed_from_analysis.txt', header=False, index=False)

        records = []
        genome_prot_table = pd.DataFrame({'genome_id': [],
                                          'genome_name': [],
                                          'prot_id': [],
                                          'sequence': []
                                          })
        power = 3

        if id_table.shape[0] > 10 ** power:

            for start in range(0, id_table.shape[0], 10 ** power):
                end = start + 10 ** power if start + 10 ** power < id_table.shape[0] else id_table.shape[0]
                print(start, end, id_table.shape[0])
                records = get_record(list_id=id_list[start:end], number=power-1)
                genome_prot_table = break_record(records, genome_prot_table)
                del records

        else:
            print(0, id_table.shape[0], id_table.shape[0])
            records += get_record(list_id=id_list, number=4)
            genome_prot_table = break_record(records, genome_prot_table)
            del records

        seeds = pd.read_table(str(input.seeds),names=["gene", "protein_id", "color"])
        with Entrez.efetch(db='protein',id=seeds.protein_id.to_list(),rettype='gb',retmode='xml') as seed_handle:
            for record in Entrez.read(seed_handle):
                if not record['GBSeq_locus'] in genome_prot_table.prot_id.to_list():
                    genome_prot_table = genome_prot_table.append({'genome_id': 'Seed',
                                                                  'genome_name': 'Seed',
                                                                  'prot_id': record['GBSeq_locus'],
                                                                  'sequence': record['GBSeq_sequence']
                                                                  }, ignore_index=True)

        del seeds
        del seed_handle

        genome_prot_table.to_csv(str(output.genome_prot_table), index=False, sep='\t')

        with open(str(output.fasta_file), 'w') as fasta_file:
            for row in genome_prot_table.iterrows():
                row = row[1]
                fasta_file.write(f'>{row["prot_id"]}\n{row["sequence"]}\n')

        print(f'{genome_prot_table.shape[0]} proteins recorded')
        del genome_prot_table
        
##########################################################################
##########################################################################

rule blast:
    input:
         f'{project_name}_all_protein.fasta'
    output:
          f"blastall_{project_name}.out"
    threads: 5
    shell:
         """
         module load ncbi_blast/2.10.1
         makeblastdb -dbtype prot -in {input}
         blastp -query {input} -db {input} -evalue 0.001 -outfmt 6 -out {output} -num_threads 5 -num_alignments 25000
         """
##########################################################################
##########################################################################

rule filter_blast:
    input:
        f"blastall_{project_name}.out"
    output:
        f"{new_dir}/blastall_{project_name}_eval{e_val}.out"
    params:
        e_val = e_val
    run:
        import pandas as pd


        header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
                  'evalue', 'bitscore']
        blast_out = pd.read_table(str(input), names=header)
        blast_out_filtered = blast_out[blast_out.evalue <= params.e_val].reset_index(drop=True)
        blast_out_filtered.to_csv(str(output), index=False, sep='\t')

##########################################################################
##########################################################################

rule silix:
    input:
         fasta = f'{project_name}_all_protein.fasta',
         blast_out = f"{new_dir}/blastall_{project_name}_eval{e_val}.out"
    output:
        f"{new_dir}/{project_name}_eval{e_val}_id{id}_cov{cov}.fnodes"
    params:
          id = id,
          cov = cov,
          e_val = e_val,
          new_dir = new_dir
    shell:
         """
         module load silix/1.2.11
         sh -c 'silix {input.fasta} {input.blast_out} -f FAM -i {params.id} -r {params.cov} >{output}'
         mkdir -p {new_dir}/family_eval{params.e_val}_id{params.id}_cov{params.cov}
         sh -c 'silix-split -p {new_dir}/family_eval{params.e_val}_id{params.id}_cov{params.cov}/cluster -n 2 {input.fasta} {output}' 
         """

##########################################################################
##########################################################################

rule make_table:
    input:
        fnodes = f"{new_dir}/{project_name}_eval{e_val}_id{id}_cov{cov}.fnodes",
        seeds = gene_tab,
        genome_prot_table = f'{project_name}_all_protein.csv'
    output:
          gene_table = f"{new_dir}/gene_table_{project_name}_eval{e_val}_id{id}_cov{cov}.csv",
          seed_table = f"{new_dir}/famillied_seeds_{project_name}_eval{e_val}_id{id}_cov{cov}.txt"
    run:
        import pandas as pd

        #take list of headers of the seeds
        seed = pd.read_table(str(input.seeds), names=["gene", "protein_id", "color"])

        #open fnodes file result of silix
        fam_id_table = pd.read_table(str(input.fnodes), names=["family", "protein_id"])

        #import the table of protein id and genome id
        gp_table = pd.read_table(str(input.genome_prot_table), sep='\t')


        #sort by family number
        fam_id_table = fam_id_table.sort_values("family", ascending=True)

        #detection of the seeds families
        seed_fam = pd.merge(seed,
                            fam_id_table[fam_id_table.protein_id.isin(seed.protein_id)],
                            how='left',
                            on='protein_id'
                            )
        seed_fam.to_csv(str(output.seed_table), index=False, sep='\t')
        seeded_fit = fam_id_table[(fam_id_table.family.isin(seed_fam.family))]

        # Dict to find gene name from family Id
        dict_fam = seed_fam.set_index('family').gene.to_dict()
        seed_list = seed.gene.tolist()

        #create final table
        patab = pd.DataFrame(columns=['genome_id', 'genome_name'] + seed_list)

        patab['genome_id'] = gp_table.genome_id
        patab['genome_name'] = gp_table.genome_name
        patab.drop_duplicates(inplace=True)
        patab.set_index('genome_id', inplace=True)
        print(patab)
        patab.drop(index='Seed', inplace=True)
        gp_table.set_index('prot_id', inplace=True)

        for index, line in seeded_fit.iterrows():
            protein = line.protein_id
            genome = gp_table.genome_id[protein]
            if genome == 'Seed':
                continue
            family = line.family

            if patab.isna().loc[genome, dict_fam[family]]:
                patab.at[genome, dict_fam[family]] = protein

            else:
                patab.at[genome, dict_fam[family]] = "{before}, {addition}".format(
                    before=patab.loc[genome, dict_fam[family]], addition=protein)
                print(patab.loc[genome, dict_fam[family]])

        patab.dropna(how="all", axis=1, inplace=True)
        sorted_patab = patab.notnull().astype('int').sort_values(by=patab.columns.tolist(), ascending=False)
        patab = patab.loc[sorted_patab.index].reset_index()

        #print the table
        patab.to_csv(str(output.gene_table), index=False, sep='\t')

##########################################################################
##########################################################################

rule plot:
    input:
         gene_table = f"{new_dir}/gene_table_{project_name}_eval{e_val}_id{id}_cov{cov}.csv",
         seed_table = f"{new_dir}/famillied_seeds_{project_name}_eval{e_val}_id{id}_cov{cov}.txt"
    output:
          pdf = f"{new_dir}/gene_table_{project_name}_eval{e_val}_id{id}_cov{cov}.pdf",
          png = f"{new_dir}/gene_table_{project_name}_eval{e_val}_id{id}_cov{cov}.png"
    run:
        import matplotlib.pyplot as plt
        import pandas as pd


        font = {'family': 'DejaVu Sans', 'weight': 'light', 'size': 10, }
        plt.rc('font', **font)
        plt.rcParams['text.color'] = 'black'
        plt.rcParams['svg.fonttype'] = 'none'  # Editable SVG text
        patab = pd.read_table(str(input.gene_table)).fillna('None').set_index('genome_id')[::-1]
        #patab.drop(labels='genome_id', axis=1, inplace=True)
        seed_table = pd.read_table(str(input.seed_table))
        tab_to_draw = patab.reset_index().drop(labels='genome_name', axis=1)\
            .melt(id_vars='genome_id').rename(columns={'variable': 'gene', 'value': 'PA'})
        dict_color = seed_table.set_index("gene").color.to_dict()
        tab_to_draw['color'] = tab_to_draw.apply(lambda x: dict_color[x.gene] if x.PA != 'None' else 'white', axis=1)
        tab_to_draw['number'] = tab_to_draw.apply(lambda x: len(x.PA.split(',')) if x.PA != 'None' else 0, axis=1)
        tab_to_draw['x_pos'] = tab_to_draw.gene.apply(lambda x: patab.columns.tolist().index(x)-1)
        tab_to_draw['y_pos'] = tab_to_draw.genome_id.apply(lambda x: patab.index.tolist().index(x))

        fig, ax = plt.subplots(1, 1, figsize=((patab.shape[1]) / 3,
                                              patab.shape[0] / 3),
                               gridspec_kw={'width_ratios': [patab.shape[1]]})

        label_format = {'color': 'black', 'fontweight': 'bold',
                        'fontsize': 12}

        for _, row in tab_to_draw.iterrows():
            ax.plot(row.x_pos, row.y_pos, linestyle="None", marker="s",
                    markersize=15, mfc=row.color, mec='black', markeredgewidth=1)

            if row.number > 1:
                ax.text(row.x_pos, row.y_pos, str(row.number), fontsize=11, color='white', ha='center', va='center',
                        fontweight='heavy')

        plt.yticks(range(patab.shape[0]), patab.genome_name.tolist(), **label_format)
        plt.xticks(range(patab.shape[1]-1), tab_to_draw.gene.unique().tolist(), **label_format)

        ax.tick_params(axis='both', which='both', length=0)  # No tick markers
        ax.set_ylabel('')  # No ylabel
        ax.xaxis.tick_top()  # xticklabels on top
        ax.xaxis.set_label_position('top')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha='center')  # Rotate x labels

        for position in ['top', 'bottom', 'left', 'right']:
            ax.spines[position].set_visible(False)  # Remove border

        plt.xlim(-0.5, patab.shape[1] - 0.5)
        plt.ylim(-0.5, patab.shape[0] - 0.5)
        plt.savefig(str(output.png), bbox_inches="tight", dpi=300)
        plt.savefig(str(output.pdf), bbox_inches="tight", dpi=300)

##########################################################################
##########################################################################
