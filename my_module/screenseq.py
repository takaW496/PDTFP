import os
import numpy as np
import pandas as pd
import math
import shutil
    
def chopseq(seq, window, step):
    """ Moving window to chop target sequence """
    seqChopped = []
    start = 0
    end = start + window
    while len(seq) >= window:
        seqChopped_tmp = []
        seqChopped_tmp.append(seq[0:window])
        seqChopped_tmp.append(start)
        seqChopped_tmp.append(end)
        seqChopped.append(seqChopped_tmp)
        start += step
        end += step
        seq = seq[step:]
    return seqChopped


def calculategc(seq):
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    gccontent = (gc / len(seq)) * 100
    return(gccontent)

def calculatetm(seq):
    """ Calculate Tm of a target candidate, nearest neighbor model """
    NNlist = chopseq(seq, 2, 1)
    NNlist = [i[0] for i in NNlist]
    NNtable = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    NNendtable = ['A', 'C', 'G', 'T']
    NNcount = np.zeros(16)
    NNend = np.zeros(4)
    for c, NN in enumerate(NNtable):
        NNcount[c] = NNlist.count(NN)
    for c, NN in enumerate(NNendtable):
        NNend[c] = seq[0].count(NN)
    # numbers below from Sugimoto et al. NAR (1996)
    NNEnthalpy = np.array([-8.0, -9.4, -6.6, -5.6, -8.2, -10.9, -11.8, -6.6, -8.8, -10.5, -10.9, -9.4, -6.6, -8.8, -8.2, -8.0])
    NNEntropy = np.array([-21.9, -25.5, -16.4, -15.2, -21.0, -28.4, -29.0, -16.4, -23.5, -26.4, -28.4, -25.5, -18.4, -23.5, -21.0, -21.9])
    NNendEnthalpy = np.array([.6, .6, .6, .6])
    NNendEntropy = np.array([-9.0, -9.0, -9.0, -9.0])

    sumEnthalpy = np.sum(np.multiply(NNcount, NNEnthalpy)) + np.sum(np.multiply(NNend, NNendEnthalpy))
    sumEntropy = np.sum(np.multiply(NNcount, NNEntropy)) + np.sum(np.multiply(NNend, NNendEntropy))
    Tm = (sumEnthalpy * 1000)/(sumEntropy + (1.9872 * math.log(1e-7))) - 273.15       # oligo concentration: 1e-7 M
    #sumSalt = 0.075 + (3.795 * 0.01**0.5)  # monovalent: 0.075 M, bivalent: 0.01 M 
    sumSalt = 0.33   #2xSSC buffer contains 0.3 M NaCl and 0.03 mM sodium citrate -> 330 mM monovalent and no bivalent
    Tm += 16.6 * math.log10(sumSalt)     # salt correction
    #Tm -= 0.72 * 20       # formamide correction
    Tm -= 0.72 * 30     #hybridization buffer contains 30% formamide
    return Tm

def search_balanced_seq_pair(seq):
    seq_1, seq_2 = seq[:len(seq)//2-1], seq[len(seq)//2+1:]
    interval = seq.replace(seq_1, '').replace(seq_2, '')
    tm_1, tm_2 = calculatetm(seq_1), calculatetm(seq_2)
    tm_diff = 1000  #start condition
    tmp_tm_diff = abs(tm_1 - tm_2)
    seq_1_tmp, seq_2_tmp, tm_1_tmp, tm_2_tmp = seq_1, seq_2, tm_1, tm_2 #initialize
    while tm_diff > tmp_tm_diff:
        tm_diff, seq_1, seq_2, tm_1, tm_2 = tmp_tm_diff, seq_1_tmp, seq_2_tmp, tm_1_tmp, tm_2_tmp #update variable with current values
        if tm_1 > tm_2:
            seq_2_tmp = interval[-1] + seq_2
            interval = interval[:-1]
            interval = seq_1[-1] + interval
            seq_1_tmp = seq_1[:-1]
        elif tm_2 > tm_1:
            seq_1_tmp = seq_1 + interval[0]
            interval = interval[1:]
            interval = interval + seq_2[0]
            seq_2_tmp = seq_2[1:]
        elif tm_1 == tm_2:
            break
        #judge diff_tm
        tm_1_tmp, tm_2_tmp = calculatetm(seq_1_tmp), calculatetm(seq_2_tmp)
        tmp_tm_diff = abs(tm_1_tmp - tm_2_tmp)

    return(seq_1, seq_2, tm_1, tm_2)

def runscreen(argin):
    """ Screen gc and write sequences that fulfill gc requirement to file, to use as blast input """
    target_gene_ID, target_seq, seed_seq_len, gcmin, gcmax, tmmin, tmmax, out_dir = argin
    print('target_gene_ID:', target_gene_ID)

    dname = os.path.join(out_dir, target_gene_ID)
    if os.path.exists(dname):
        shutil.rmtree(dname)
    os.makedirs(dname, exist_ok=True)

    siteChopped= {}
    id_list, seq_list, seq_start_list, seq_end_list, GC_list, seq_1_list, seq_2_list, tm_1_list, tm_2_list = \
    [], [], [], [], [], [], [], [], []

    listSeqChopped = chopseq(target_seq, seed_seq_len, 5)
    #evaluate seed seqs based on GC content
    for j, seqChopped in enumerate(listSeqChopped):
        gc_tmp = calculategc(seqChopped[0])
        seqChopped.append(gc_tmp)
        if gcmin < gc_tmp < gcmax:  # GC thresholding
            siteChopped[j] = seqChopped #seq, start, end, gc
            
    c = 0
    balanced_tm_pair = []
    #evaluate seed seqs based on Tm
    for key, value in siteChopped.items():
        seq_1, seq_2, tm_1, tm_2 = search_balanced_seq_pair(value[0])
        balanced_tm_pair.append([tm_1, tm_2])
        #add to the list if balanced seq pair meets the criteria
        if tmmin < tm_1 < tmmax and tmmin < tm_2 < tmmax:
            id_list += [key]
            seq_list += [value[0]]
            seq_start_list += [value[1]]
            seq_end_list += [value[2]]
            GC_list += [value[3]]
            seq_1_list += [seq_1]
            seq_2_list += [seq_2]
            tm_1_list += [tm_1]
            tm_2_list += [tm_2]

            # write files that can be used as input in blastn
            with open(os.path.join(dname, 'target_' + str(key) + '.fasta'), 'w') as fblast:
                fblast.write(">target_%d\n%s\n" % (key, value[0]))
    highest_tm = max([max(i) for i in balanced_tm_pair], default=0) #return 0 if empty
    lowest_tm = min([min(i) for i in balanced_tm_pair], default=0) #return 0 if empty
    if highest_tm == 0 and lowest_tm == 0:
        print('balanced_tm_pair is empty.')
    elif highest_tm < tmmin:
        print('The maximum Tm found in all balanced sequence pairs:', highest_tm)
        print('There is no balanced pairs with Tm above tmmin. Try to increase the seed_seq_len parameter.')
    elif lowest_tm > tmmax:
        print('The minimum Tm found in all balanced sequence pairs:', lowest_tm)
        print('There is no balanced pairs with Tm below tmmax. Try to decrease the seed_seq_len parameter.')

    df = pd.DataFrame(data = {'sequence_id':id_list,'target_gene_ID':target_gene_ID,'seed_sequence':seq_list,'start_pos':seq_start_list,'end_pos':seq_end_list,
        'GC (%)':GC_list,'1st_sequence': seq_1_list, '2nd_sequence': seq_2_list, '1st_Tm': tm_1_list, '2nd_Tm': tm_2_list}, 
    columns = [ 'sequence_id', 'target_gene_ID', 'seed_sequence', 'start_pos', 'end_pos', 'GC (%)', '1st_sequence', '2nd_sequence', '1st_Tm', '2nd_Tm'])
    print(f'number of candidate sequence for {target_gene_ID} is', (df['target_gene_ID'] == target_gene_ID).sum())
    return df

def search_seed_seq (target_seqs, target_gene_IDs, seed_seq_len, gcmin, gcmax, tmmin, tmmax, out_dir):
    """ Parallel process of GC thresholding and subsequent Tm thresholding for input sequences """
    # pack up inputs for multiprocessing
    inputs = []
    for c, i in enumerate(target_gene_IDs):
        inputs.append((i, target_seqs[c], seed_seq_len, gcmin, gcmax, tmmin, tmmax, out_dir))

    siteCandidates = pd.DataFrame()

    for i in inputs:
        argout = runscreen(i)
        if len(argout) < 1:
            print(f'No sequences satisfied criteria for this gene {i[0]}. Please adjust parameters to find seed sequences.')
            continue
        if len(siteCandidates) < 0: 
            siteCandidates = argout
        else:
            siteCandidates = pd.concat([siteCandidates, argout])

    return siteCandidates