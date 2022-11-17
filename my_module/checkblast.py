# for multi_padlock_design
# Xiaoyan, 2017

import os
import numpy as np
import pandas as pd

def readblastout(file, seed_seq_len, gene_id, seq_id):
    """ Read and check the blast results and return boolean varialbe specific """
    specific = True
    with open(file, 'r') as fh:
        for line in fh:
            if specific:
                # scores = []
                columns = line.split(',')
                if columns[1] == gene_id and columns[2] == 100 and columns[3] == seed_seq_len:
                    #skip the hit to itself
                    continue

                if seed_seq_len*.5 < int(columns[3]) < seed_seq_len and float(columns[2]) > 80:
                    # more than 50% coverage, 80% homology, ignore ligation site here
                    specific = False
    return specific


def getcandidates(df, target_gene_IDs, dir_blast_out, seed_seq_len):
    """ judge blast resutls and add is_specific column to df """
    is_specific = {}
    for gene in target_gene_IDs:
        blast_output_files = os.listdir(os.path.join(dir_blast_out, gene))

        for f in blast_output_files:
            seq_id = int(f.split('_')[1])
            file_blast = os.path.join(dir_blast_out, gene, f)
            is_specific[seq_id] = readblastout(file_blast, seed_seq_len, gene, seq_id)
           
    df_specific=pd.DataFrame(data={'sequence_id': is_specific.keys(),'is_specific': is_specific.values()})
    df = pd.merge(df, df_specific, on = 'sequence_id', how = 'left')

    return df

def add_id_clust(df, target_gene_IDs, seed_seq_len, step):
    clust_id_list = []
    for gene in target_gene_IDs:
        df_tmp = df[df['target_gene_ID'] == gene]
        tmp_list = []
        num_clust = 0
        for i in df_tmp.sequence_id:
            if num_clust == 0:
                num_clust += 1
                tmp_list.append(num_clust)
                previous_id = i
                max_id = i + seed_seq_len//step
            else:
                id_diff = i - previous_id
                if i > max_id:
                    num_clust += 1
                    max_id = i + seed_seq_len//step
                    tmp_list.append(num_clust)
                    previous_id = i
                else:
                    tmp_list.append(num_clust)
                    previous_id = i
        clust_id_list += tmp_list
    df.loc[:, 'id_clust'] = clust_id_list
    return df