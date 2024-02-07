#finalize probe module
from Bio.Seq import Seq
import pandas as pd
import os
import re
import subprocess
import shutil
from my_module import parblast

Processes = []
NextProcess = 0
input_files = 0

def add_probe_seq(df, target_gene_IDs, target_genome, round_table, BridgeProbes, dye_seq, dir_blast_in, dir_blast_out, high_throughput):
	"""  finalize probe sequence for each candidate seed sequences """
	df_bridge = check_bridge_seq(BridgeProbes, target_genome, dir_blast_in, dir_blast_out)
	df = add_PLP_column(df, target_gene_IDs, df_bridge, round_table, high_throughput)
	df = add_primer_column(df)
	df = add_bridge_probe_column(df, target_gene_IDs, round_table, dye_seq)

	return df


def add_bridge_probe_column(df, target_gene_IDs, round_table, dye_seq):
	df.loc[:,'dye_seq'] = get_dye_seq(df.loc[:,'Seq_Barcode'])
	df_round = pd.read_table(round_table, sep = '\t')
	df_dye = pd.read_table(dye_seq, sep = '\t')
	df_round = pd.merge(df_round, df_dye, on = 'dye', how = 'left')
	df_round.loc[:,'Seq_rev_comp'] = [str(Seq(seq).reverse_complement()) for seq in df_round['detection_oligo']]
	df_round.loc[df_round['dye'] == 'AF750','Seq_rev_comp'] = 'C' + df_round.loc[:, 'Seq_rev_comp'] #irregular design
	df = pd.merge(df,df_round, left_on = 'target_gene_ID', right_on = 'gene_ID')
	df.loc[:,'dye_seq'] = df.loc[:,'dye_seq'] + df.loc[:,'Seq_rev_comp']
	df = df.drop(['gene_ID', 'Seq_rev_comp'], axis = 1)
	df = df.sort_values(by=['round', 'dye'])

	return df

def get_dye_seq(Seq_Barcode):
	dye_seq = [seq[1:] for seq in Seq_Barcode]
	endNN = [seq[-2:] for seq in dye_seq]
	endNN = [str(Seq(seq).complement()) for seq in endNN]
	dye_seq = [seq[:-2] for seq in dye_seq]
	dye_seq = [x + y for (x, y) in zip(dye_seq, endNN)]

	return dye_seq

def add_PLP_column(df, target_gene_IDs, df_bridge, round_table, high_throughput):
	PLP_seq = [str(Seq(seq).reverse_complement()) for seq in df['1st_sequence']]
	PLP_seq = ['/5Phos/ACATTA' + seq for seq in PLP_seq]
	df_bridge.loc[:,'concat_seq'] = df_bridge.loc[:,'Seq_Barcode'] + df_bridge.loc[:,'Seq_Anchor'] + 'AAGATA'
	# high throughput mode
	if high_throughput:
		print ('high_throughput mode is on. use only the same set of bridge sequences for different rounds.')
		df_Bp_ID = df_bridge.iloc[0:high_throughput]
		df_round = pd.read_table(round_table, sep = '\t')
		df_Bp_ID.loc[:, 'dye'] = df_round.loc[df_round['round'] == 1, 'dye'].tolist()
		num_repeat = int(len(target_gene_IDs)/high_throughput)
		repeated_dfs = [df_Bp_ID.assign(round=i+1) for i in range(num_repeat)]
		df_Bp_ID = pd.concat(repeated_dfs)
		df_Bp_ID.reset_index(drop=True, inplace=True)
		df_Bp_ID = pd.merge(df_Bp_ID, df_round, on = ['round', 'dye'] )
		df = pd.merge(df, df_Bp_ID, left_on = 'target_gene_ID', right_on = 'gene_ID')
		df.loc[:,'PLP_seq'] = [x + y for (x, y) in zip(PLP_seq, df['concat_seq'])]
		df = df.drop(['concat_seq','gene_ID', 'dye', 'round'], axis = 1)
        
	elif not high_throughput:
		print ('high_throughput mode is off. use different bridge sequences set for different rounds.')
		df_Bp_ID = df_bridge.iloc[0:len(target_gene_IDs)]
		df_Bp_ID = df_Bp_ID.assign(target_gene_ID = target_gene_IDs)
		df = pd.merge(df, df_Bp_ID, on = 'target_gene_ID', )
		df.loc[:,'PLP_seq'] = [x + y for (x, y) in zip(PLP_seq, df['concat_seq'])]
		df = df.drop('concat_seq', axis = 1)

	return df

def add_primer_column(df):
	primer_seq = [str(Seq(seq).reverse_complement()) for seq in df['2nd_sequence']]
	primer_seq = [seq + 'TAATGTTATCTT' for seq in primer_seq]
	df.loc[:,'primer_seq'] = primer_seq

	return df

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
