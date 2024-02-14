# for multi_padlock_design
# Xiaoyan, 2017

import os
import numpy as np
import pandas as pd

def readblastout(file, seed_seq_len, gene_id, seq_id, spec_coverage, spec_homology):
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

                if seed_seq_len*spec_coverage < int(columns[3]) < seed_seq_len and float(columns[2]) > spec_homology*100:
                    # more than 50% coverage, 80% homology, ignore ligation site here
                    specific = False
    return specific


def getcandidates(df, target_gene_IDs, dir_blast_out, seed_seq_len, spec_coverage, spec_homology):
    """ judge blast resutls and add is_specific column to df """
    if not isinstance(seed_seq_len, list): #convert seed_seq_len to list for loop if it is a integer.
        seed_seq_len = [seed_seq_len]
    is_specific = {}
    for gene in target_gene_IDs:
        for l in seed_seq_len:
            dname = gene + "_" + str(l)
            blast_output_files = os.listdir(os.path.join(dir_blast_out, dname))

            for f in blast_output_files:
                seq_id = int(f.split('_')[1])
                file_blast = os.path.join(dir_blast_out, dname, f)
                is_specific[seq_id] = readblastout(file_blast, l, gene, seq_id, spec_coverage, spec_homology)
               
    df_specific=pd.DataFrame(data={'sequence_id': is_specific.keys(),'is_specific': is_specific.values()})
    df = pd.merge(df, df_specific, on = 'sequence_id', how = 'left')

    return df


def add_id_clust(df, seed_seq_len):
    # Function to assign cluster IDs within each gene
    def assign_cluster_ids(group):
        # Sort the DataFrame by start position
        group_sorted = group.sort_values(by='start_pos')

        # Initialize cluster dictionary
        cluster_dict = {}
        cluster_id = 0
        cluster_ids = []
        cluster_start = None

        # Iterate over each row in the DataFrame
        for index, row in group_sorted.iterrows():
            # If cluster_start is not set or the start position is beyond 50 bp from the current cluster start
            if cluster_start is None or row['start_pos'] - cluster_start > seed_seq_len:
                # Start a new cluster
                cluster_id += 1
                cluster_start = row['start_pos']
                cluster_ids.append(cluster_id)
            else:
                cluster_ids.append(cluster_id)

        group_sorted['id_clust'] = cluster_ids
        return group_sorted

    df = df.groupby('target_gene_ID').apply(assign_cluster_ids)
    df = df.reset_index(drop=True)

    return df