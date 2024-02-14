#finalize probe module
from Bio.Seq import Seq
import pandas as pd
import os
import re

Processes = []
NextProcess = 0
input_files = 0

def add_probe_seq(df, df_design):
	"""  finalize probe sequence for each candidate seed sequences """
	df = add_PLP_column(df, df_design)
	df = add_primer_column(df)
	df = add_bridge_probe_column(df)

	return df


def add_bridge_probe_column(df):
	df.loc[:,'dye_seq'] = get_dye_seq(df.loc[:,'Seq_Barcode'])
	df.loc[:,'Seq_rev_comp'] = [str(Seq(seq).reverse_complement()) for seq in df['detection_oligo']]
	df.loc[df['dye'] == 'AF750','Seq_rev_comp'] = 'C' + df.loc[:, 'Seq_rev_comp'] #irregular design
	df.loc[:,'dye_seq'] = df.loc[:,'dye_seq'] + df.loc[:,'Seq_rev_comp']
	df = df.drop('Seq_rev_comp', axis = 1)
	df = df.sort_values(by=['round', 'dye'])

	return df

def get_dye_seq(Seq_Barcode):
	dye_seq = [seq[1:] for seq in Seq_Barcode]
	endNN = [seq[-2:] for seq in dye_seq]
	endNN = [str(Seq(seq).complement()) for seq in endNN]
	dye_seq = [seq[:-2] for seq in dye_seq]
	dye_seq = [x + y for (x, y) in zip(dye_seq, endNN)]

	return dye_seq

def add_PLP_column(df, df_design):
	PLP_seq = [str(Seq(seq).reverse_complement()) for seq in df['1st_sequence']]
	PLP_seq = ['/5Phos/ACATTA' + seq for seq in PLP_seq]
	df = pd.merge(df, df_design, on = 'target_gene_ID')
	df.loc[:,'concat_seq'] = df.loc[:,'Seq_Barcode'] + df.loc[:,'Seq_Anchor'] + 'AAGATA'
	df.loc[:,'PLP_seq'] = [x + y for (x, y) in zip(PLP_seq, df['concat_seq'])]
	df = df.drop('concat_seq', axis = 1)

	return df

def add_primer_column(df):
	primer_seq = [str(Seq(seq).reverse_complement()) for seq in df['2nd_sequence']]
	primer_seq = [seq + 'TAATGTTATCTT' for seq in primer_seq]
	df.loc[:,'primer_seq'] = primer_seq

	return df

def assess_output(df, round_table):
	specific_seq = df[df['is_specific'] == True]
	result = specific_seq.groupby('target_gene_ID')['id_clust'].nunique()
	ids_mt_4 = result[result >= 4].index.tolist()
	df_round = pd.read_table(round_table, sep = '\t')
	missing_id = list(set(df_round['target_gene_ID']) - set(ids_mt_4))
	if not missing_id:
		print("All genes have enough number of candidates")
		return missing_id
	elif missing_id:
	    print("These genes do not have enough number of specific candidates that belong to different id_clust.")
	    print(missing_id)

	return missing_id
