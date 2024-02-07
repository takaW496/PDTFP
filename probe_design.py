#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy
import pandas
import argparse
import os
import re
import subprocess
import multiprocessing
import sys
from Bio.Seq import Seq
#own module
from my_module import screenseq
from my_module import parblast
from my_module import checkblast
from my_module import finalize_probe
from my_module import load_db


# In[2]:


#parameter settings
dir_blast_in = 'blast_in'
dir_blast_out = 'blast_out'
target_genome = 'blastDB/Lotusjaponicus_test.fa' # modified
round_table = 'round_table_Lj_MO.txt' # modified
BridgeProbes = 'BridgeProbes.txt'
dye_seq_file = 'dye_seq.txt'
seed_seq_len = 50
gcmin = 40
gcmax = 60
tmmin = 48
tmmax = 52
outfile = 'Lj_candidates_test2.txt'
threads = 1
args = dict()
args['target_genome'] = target_genome
args['round_table'] = round_table
args['BridgeProbes'] = BridgeProbes
args['dye_seq'] = dye_seq_file
args['seed_seq_len'] = seed_seq_len
args['gcmin'] = gcmin
args['gcmax'] = gcmax
args['tmmin'] = tmmin
args['tmmax'] = tmmax
args['outfile'] = outfile
args['threads'] = threads


# In[4]:

os.makedirs(dir_blast_in, exist_ok=True)
os.makedirs(dir_blast_out, exist_ok=True)

#main workflow
target_gene_IDs = load_db.readgenefile(args['round_table'])
dict_ID2Seq = load_db.load_fastadb(args['target_genome'])
target_seqs = load_db.findseq(target_gene_IDs, dict_ID2Seq)
df = screenseq.search_seed_seq(target_seqs, target_gene_IDs, args['seed_seq_len'], args['gcmin'], args['gcmax'], args['tmmin'], args['tmmax'], dir_blast_in)
parblast.blastdb(target_genome) #prepare blast db
parblast.continueblast(target_gene_IDs, target_genome, dir_blast_in, dir_blast_out)
df = checkblast.getcandidates(df, target_gene_IDs, dir_blast_out, args['seed_seq_len'])
df = checkblast.add_id_clust(df, target_gene_IDs, args['seed_seq_len'], 5)
df = finalize_probe.add_probe_seq(df, target_gene_IDs, args['target_genome'], args['round_table'], args['BridgeProbes'], args['dye_seq'], dir_blast_in, dir_blast_out)
df.to_csv(args['outfile'], sep = '\t', index = False)
print(f"Finished all processes. please check {args['outfile']}")


