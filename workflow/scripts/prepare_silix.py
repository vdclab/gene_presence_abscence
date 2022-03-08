import pandas as pd
import sys
import numpy as np 

##########################################################################
##########################################################################

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################
# Try to load it from file but not good now because can't access to the file
# if deployed worlflow by snakedeploy
##########################################################################

def iterrator_on_blast_hsp(blast_out):
    """Iterate over blast output. It considers that the output
    is in outfmt 6 and that all the hsp should be one after the
    other

    Args:
        blast_out (string): Path to a blast tabulated output

    Yields:
        pandas.DataFrame: A pandas Dataframe of the two proteins and their hsps
    """

    current_pair = ''
    list_line = []

    with open(blast_out) as r_file:

        for line in r_file:
            split_line = line.rstrip().split('\t')

            if (split_line[0], split_line[1]) == current_pair:
                list_line.append(split_line)
            elif current_pair == '':
                current_pair = (split_line[0], split_line[1])
                list_line = [split_line]
            else:
                yield list_line
                current_pair = (split_line[0], split_line[1])
                list_line = [split_line]

        yield list_line

##########################################################################

def get_orientation(sstart, send):
    """Returns the orientation between two points.

    Args:
        sstart (int): Start of the sequence
        send (int): End of the sequence

    Returns:
        int: 1 or -1 depending on the orientation
    """
    
    return ((send - sstart) / abs((send - sstart))).astype(int)

##########################################################################

def update_values(old_values, sens):
    """Update sequence start/end with depending of the orientation.

    Args:
        old_values (int): Position in the sequence
        sens (int): 1 or -1 depending on the orientation

    Returns:
        int: Position in the sequence
    """

    return old_values * sens

##########################################################################

def length_aln_on_sequence(start, end):
    """Return the length of the sequence.

    Args:
        start (int): Start of the sequence
        end (int): End of the sequence

    Returns:
        int: Length of the protein
    """
    return end - start + 1

##########################################################################

def calculate_fraction(delta, lgHSP, bitscore, pid, pos):
    """Calculate the new score, identities and positives of the hsp2.

    Args:
        delta (int): Difference of size between old and new hsp2 length
        lgHSP (int): size of the hsp before trimming
        bitscore (float): Score of the alignment
        pid (int): Number of identical in the alignment 
        pos (int): Number of positives in the alignment

    Returns:
        float, int, int, int: The updated values of the score, identities, positives and length
    """

    # Calculate new score, id and positive
    # Calculation: initial_value * (franction of length that has been preserved)
    fraction = 1 - (delta / lgHSP)
    
    new_score = np.floor(bitscore * fraction)
    new_id = np.floor(pid * fraction)
    new_pos = np.floor(pos * fraction)
    
    # Calculate new length
    new_length = lgHSP - delta
    
    # Set expect value to -1 : this value should not be used after
    # having changed HSPs boundaries
    # new_evalue = -1
    
    return new_score, new_id, new_pos, new_length

##########################################################################

def remove_overlap_query(hsp1, hsp2):
    """Remove overlap with the given hsp2 in the query.

    Args:
        hsp1 (pandas.Series): Blast values of the given hsp (best)
        hsp2 (pandas.Series): Blast values of the given hsp (questioning)

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    # Calculate where is the overlap and remove the overlapping part
    if hsp2.qstart < hsp1.qstart:
        new_qstart = hsp2.qstart
        new_qend = hsp1.qstart -1
        delta = hsp2.qend - new_qend
        new_sstart = hsp2.sstart
        new_send = hsp2.send - delta
        
    elif hsp2.qend > hsp1.qend:
        new_qstart = hsp1.qend + 1
        new_qend = hsp2.qend
        delta = new_qstart - hsp2.qstart 
        new_sstart = hsp2.sstart + delta
        new_send = hsp2.send
        
    new_score, new_id, new_pos, new_length = calculate_fraction(delta=delta, 
                                                                lgHSP=hsp2.lgHSP,
                                                                bitscore=hsp2.bitscore, 
                                                                pid=hsp2.id,
                                                                pos=hsp2.pos)
        
    return {'bitscore':new_score, 'length':new_length,
            'qend':new_qend, 'qstart':new_qstart, 'send':new_send,
            'sstart':new_sstart, 'pos':new_pos, 'id':new_id}

##########################################################################

def remove_overlap_subject(hsp1, hsp2):
    """Remove overlap with the given hsp2 in the subject.

    Args:
        hsp1 (pandas.Series): Blast values of the given hsp (best)
        hsp2 (pandas.Series): Blast values of the given hsp (questioning)

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    # Calculate where is the overlap and remove the overlapping part
    if hsp2.sstart < hsp1.sstart:
        new_sstart = hsp2.sstart
        new_send = hsp1.sstart -1
        delta = hsp2.send - new_send
        new_qstart = hsp2.qstart
        new_qend = hsp2.qend - delta
        
    elif hsp2.send > hsp1.send:
        new_sstart = hsp1.send + 1
        new_send = hsp2.send
        delta = new_sstart - hsp2.sstart 
        new_qstart = hsp2.qstart + delta
        new_qend = hsp2.qend
        
    new_score, new_id, new_pos, new_length = calculate_fraction(delta=delta, 
                                                                lgHSP=hsp2.lgHSP,
                                                                bitscore=hsp2.bitscore, 
                                                                pid=hsp2.id,
                                                                pos=hsp2.pos)
        
    return {'bitscore':new_score, 'lgHSP':new_length,
            'qend':new_qend, 'qstart':new_qstart, 'send':new_send,
            'sstart':new_sstart, 'pos':new_pos, 'id':new_id}

##########################################################################

def checkHSPS(hsp1, hsp2, HSPMIN=100):
    """compare two HSPS in the blast output

    Args:
        hsp1 (pandas.Series): Blast values of the given hsp (best)
        hsp2 (pandas.Series): Blast values of the given hsp (questioning)
        HSPMIN (int, optional): Minumum length of the HSP. Defaults to 100.

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    dict_update = {'stat':0}
    
    # Check if the hsp2 is in a different orientation than hsp1
    if hsp1.sens != hsp2.sens:
        # print(f'orientation wrong: {hsp1.sens} != {hsp2.sens}')
        return dict_update
    
    # Check is hsp2 inside hsp1 for the query sequence
    if hsp1.qstart >= hsp2.qstart and hsp1.qend <= hsp2.qend:
        # print(f'hsp2 inside hsp1 for query: {hsp1.qstart} >= {hsp2.qstart} and {hsp1.qend} <= {hsp2.qend}')
        return dict_update
    
    # Check is hsp1 inside hsp2 for the query sequence
    elif hsp1.qstart <= hsp2.qstart and hsp1.qend >= hsp2.qend:
        # print(f'hsp1 inside hsp2 for query: {hsp1.qstart} <= {hsp2.qstart} and {hsp1.qend} >= {hsp2.qend}')
        return dict_update  

    # Check is hsp1 inside hsp2 for the subject sequence
    elif hsp1.sstart >= hsp2.sstart and hsp1.send <= hsp2.send:
        # print(f'hsp2 inside hsp1 for subject: {hsp1.sstart} >= {hsp2.sstart} and {hsp1.send} <= {hsp2.send}')
        return dict_update
    
    # Check is hsp2 inside hsp1 for the subject sequence
    elif hsp1.sstart <= hsp2.sstart and hsp1.send >= hsp2.send:
        # print(f'hsp1 inside hsp2 for subject: {hsp1.sstart} <= {hsp2.sstart} and {hsp1.send} >= {hsp2.send}')
        return dict_update 

    # reject HSPs that are in different orientation:
    # Query:    ---- A ---- B -----   A = HSP1 
    # Sbjct:    ---- B ---- A -----   B = HSP2

    if (hsp1.qend - hsp2.qend) * (hsp1.send - hsp2.send) < 0:
        # print(f'HSPs are in different orientation 1: ({hsp1.qend} - {hsp2.qend}) * ({hsp1.send} - {hsp2.send}) ===> {(hsp1.qend - hsp2.qend) * (hsp1.send - hsp2.send)} < 0')
        return dict_update
    elif (hsp1.qstart - hsp2.qstart) * (hsp1.sstart - hsp2.sstart) < 0:
        # print(f'HSPs are in different orientation 2: ({hsp1.qstart} - {hsp2.qstart}) * ({hsp1.sstart} - {hsp2.send}) ===> {(hsp1.qstart - hsp2.qstart) * (hsp1.sstart - hsp2.sstart)} < 0')
        return dict_update
    
    overlap_q = (hsp2.qstart - hsp1.qend) * (hsp2.qend - hsp1.qstart)
    overlap_s = (hsp2.sstart - hsp1.send) * (hsp2.send - hsp1.sstart)
    
    # Accept non-overlapping HSPs in correct orientation
    if overlap_q > 0 and overlap_s > 0 :
        # print(f'No overlap in query and subject: {overlap_q} > 0 and {overlap_s} > 0')
        dict_update['stat'] = 1
        return dict_update
    
    # Test if the query is overlaping
    elif overlap_q < 0:
        # print(f'Overlap in query: {overlap_q} > 0')
        dict_update = remove_overlap_query(hsp1=hsp1, 
                                           hsp2=hsp2)
        hsp2.update(dict_update)
        overlap_s = (hsp2.sstart - hsp1.send) * (hsp2.send - hsp1.sstart)
    
    # Test if the subject is overlaping after possible update of an overlaping query
    if overlap_s < 0:
        # print(f'Overlap in subject: {overlap_s} > 0')
        dict_update = remove_overlap_subject(hsp1=hsp1, 
                                             hsp2=hsp2)
        hsp2.update(dict_update)
    
    # Filter out HSPs that are too short
    if hsp2.lgHSP < HSPMIN:
        # print(f'HSP too short: {hsp2.lgHSP} < {HSPMIN}')
        dict_update['stat'] = 0
        return dict_update
    
    # Set status to 1 for consistent HSPs
    dict_update['stat'] = 1
    return dict_update


##########################################################################

def prepare_df_hsps(list_hsps, blast_names, blast_dtypes, HSPMIN=100):
    """Prepare a dataframe containing HSPs and update values.

    Args:
        list_hsps (list of strings): List of the lines containing the HSPs for the same pair of proteins
        blast_names (list of strings): List of the header in the blast outfmt 6 output
        blast_dtypes (dict): Type of the different columns in the blast outfmt 6 output
        HSPMIN (int, optional): [description]. Defaults to 100.

    Returns:
        pandas.DataFrame: Reduce DataFrame with only the significant hsps
    """
    df_hsps = pd.DataFrame(list_hsps, columns=blast_names)
    df_hsps = df_hsps.astype(blast_dtypes)
    
    # Calculate a more precise percentage of identity for each independent HSP
    df_hsps['id'] = np.floor(df_hsps.pident.values * df_hsps.length.values / 100 + 0.5)
    df_hsps['pos'] = df_hsps['id']
    
    # Get the sens of the HSP on the subject
    df_hsps['sens'] = get_orientation(sstart=df_hsps.sstart.values, send=df_hsps.send.values)
    
    # Update the values of the position in case the sens is changed
    df_hsps['sstart'] = update_values(old_values=df_hsps.sstart.values, sens=df_hsps.sens.values)
    df_hsps['send'] = update_values(old_values=df_hsps.send.values, sens=df_hsps.sens.values)
    
    # Calculate the length of the alignment on the subject and query
    df_hsps['lg_aln_subject'] = length_aln_on_sequence(start=df_hsps.sstart.values, end=df_hsps.send.values)
    df_hsps['lg_aln_query'] = length_aln_on_sequence(start=df_hsps.qstart.values, end=df_hsps.qend.values)
    
    # Calculate length of the HSP
    df_hsps['lgHSP'] = np.max([df_hsps.lg_aln_subject.values, df_hsps.lg_aln_query.values], axis=0)
    
    # Order HSPS by decreasing score
    df_hsps = df_hsps.sort_values('bitscore', ascending=False).reset_index(drop=True)
    
    # Put Stats value to know which HSPS to concerve
    df_hsps['stat'] = [1] + [-1]*(df_hsps.shape[0] - 1)
    
    for hsp1_index in range(df_hsps.shape[0] - 1):
        for hsp2_index in range(hsp1_index + 1, df_hsps.shape[0]):
            if df_hsps.at[hsp1_index, 'stat'] and df_hsps.at[hsp2_index, 'stat'] :                
                values2update = checkHSPS(hsp1=df_hsps.iloc[hsp1_index].copy(), 
                                          hsp2=df_hsps.iloc[hsp2_index].copy(), 
                                          HSPMIN=HSPMIN)

                for value in values2update:
                    df_hsps.at[hsp2_index, value] = values2update[value]
    
    return df_hsps[df_hsps.stat == 1]


##########################################################################
##########################################################################

# Seeds
seed_table = (
        pd.read_table(snakemake.input.seed_file)
        .set_index('protein_id')
)

# Protein lenght dict
# protein_dict = {}

# with open(snakemake.input.protein_file, 'rt') as r_file :
#     r_file.readline()
#     for line in r_file:
#         r_line = line.rstrip()
#         split_line = r_line.split('\t')
        
#         protein_dict[split_line[0]] = int(split_line[-1])

# protein_dict.update(seed_table.length
#                               .to_dict())

# Blast_out and preparation
blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']

# Get the types of te columns for multiple HSPs dataframe
blast_dtypes = {'qseqid':'string',
                'sseqid':'string',
                'pident':'float64',
                'length':'int32',
                'mismatch':'int32',
                'gapopen':'int32',
                'qstart':'int32',
                'qend':'int32',
                'sstart':'int32',
                'send':'int32',
                'evalue':'float64',
                'bitscore':'float64'}

# Calculating the max and min value of interest
max_eval = seed_table.evalue.max()
# min_pident = seed_table.pident.min()
# min_coverage = seed_table.coverage.min()

# Open all the output
list_open_output = [open(tsv_output, 'wt') for tsv_output in snakemake.output]
num_output = len(list_open_output)

# Read the blast hsp by hsp
for sub_blast in iterrator_on_blast_hsp(blast_out=snakemake.input.blast_out) :
    # Variable to get the lines 
    line = ""

    # Get the number of hsps
    num_HSPs = len(sub_blast)

    if num_HSPs == 1:
        evalue_blast = float(sub_blast[0][blast_names.index('evalue')])
        line = "\t".join(sub_blast[0]) + "\n"
    else:
        df_hsps = prepare_df_hsps(list_hsps = sub_blast,
                                blast_dtypes = blast_dtypes, 
                                blast_names = blast_names,
                                HSPMIN = snakemake.params.minimum_length)

        evalue_blast = df_hsps.evalue.max()

        for index in df_hsps.index:
            line += "\t".join(sub_blast[index]) + "\n"                

    if evalue_blast <= max_eval :
    #if evalue_blast <= max_eval and pident_blast >= min_pident :
        # coverage_blast = length / protein_dict[qseqid]

        # Calculating the coverage on the query
        # if coverage_blast >= min_coverage :
        for index in range(num_output) :
            file_out_path = snakemake.output[index]

            constrains = file_out_path.split('--')[-1].split('.out')[0]

            constrains_split = constrains.split('_')

            # As i have the control about the end of the file name it is better to go this way in case "_" in the gene name
            evalue = float(constrains_split[-5])
            # coverage = float(constrains_split[-3])
            # pident = float(constrains_split[-1])

            # if evalue_blast <= evalue and pident_blast >= pident and coverage_blast >= coverage:
            if evalue_blast <= evalue:
                list_open_output[index].write(line)

for file_open in list_open_output :
    file_open.close()
