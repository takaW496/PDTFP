#get experimental design module
from Bio.Seq import Seq
import pandas as pd
import os
import re
import shutil
from my_module import parblast

Processes = []
NextProcess = 0
input_files = 0


def get_exp_design(target_genome, target_gene_IDs, round_table, BridgeProbes, dye_seq, dir_blast_in, dir_blast_out, high_throughput):

	df_bridge = check_bridge_seq(BridgeProbes, target_genome, dir_blast_in, dir_blast_out)
	df_round = pd.read_table(round_table, sep = '\t')
	df_dye = pd.read_table(dye_seq, sep = '\t')
	df_design = assign_gene(target_gene_IDs, df_bridge, df_round, df_dye, high_throughput)

	return df_design


def check_bridge_seq(BridgeProbes, target_genome, dir_blast_in, dir_blast_out):
	df_bridge = pd.read_table(BridgeProbes, sep = '\t')
	parblast.blast_bridge_seq(df_bridge, target_genome, dir_blast_in, dir_blast_out)
	df_bridge = rank_bridge_seq(df_bridge, target_genome, dir_blast_out)

	return df_bridge

def rank_bridge_seq(df_bridge, target_genome, dir_blast_out):
	""" rank barcode and anchor sequences based on blast results against target genome """
	blast_dir = os.path.join(dir_blast_out, 'Bridge_Anchor', os.path.basename(target_genome).split('.')[0])
	blast_files = os.listdir(blast_dir)
	Bp_ID_list, bit_score_list = [], []
	for f in blast_files:
		Bp_ID = f.replace('_blast.txt', '')
		bit_score = [x.split(',')[11] for x in open(os.path.join(blast_dir, f)).readlines()]
		if not bit_score:
			print('blast_result is empty for this sequences: ', Bp_ID)
			print('assign 0 to bit_score.')
			bit_score = [0]
		else:
			bit_score = [s.strip() for s in bit_score]
		Bp_ID_list.append(Bp_ID)
		bit_score_list.append(max(bit_score))
	df_score = pd.DataFrame(data = {'BridgeProbe_ID':Bp_ID_list,'highest_bit_score':list(map(float, bit_score_list))}, columns = ['BridgeProbe_ID', 'highest_bit_score'])
	df_bridge = pd.merge(df_bridge, df_score, on = 'BridgeProbe_ID')
	df_bridge.loc[:,'Rank'] = df_bridge['highest_bit_score'].rank()
	df_bridge = df_bridge.sort_values(by = 'BridgeProbe_ID')
	df_bridge = df_bridge.sort_values(by = 'Rank')

	return df_bridge

def assign_gene(target_gene_IDs, df_bridge, df_round, df_dye, high_throughput):
	if high_throughput:
		print ('high_throughput mode is on. use only the same set of bridge sequences for different rounds.')
		df_Bp_ID = df_bridge.iloc[0:high_throughput]
		df_Bp_ID.loc[:, 'dye'] = df_round.loc[df_round['round'] == 1, 'dye'].tolist()
		num_repeat = int(len(target_gene_IDs)/high_throughput)
		repeated_dfs = [df_Bp_ID.assign(round=i+1) for i in range(num_repeat)]
		df_Bp_ID = pd.concat(repeated_dfs)
		df_Bp_ID.reset_index(drop=True, inplace=True)
		df_design = pd.merge(df_Bp_ID, df_round, on = ['round', 'dye'] )
        
	elif not high_throughput:
		print ('high_throughput mode is off. use different bridge sequences set for different rounds.')
		df_Bp_ID = df_bridge.iloc[0:len(target_gene_IDs)]
		df_Bp_ID = df_Bp_ID.assign(target_gene_ID = target_gene_IDs)
		df_design = pd.merge(df_Bp_ID, df_round, on = 'target_gene_ID', how = 'left')

	df_design = pd.merge(df_design, df_dye, on = 'dye', how = 'left')

	return df_design
