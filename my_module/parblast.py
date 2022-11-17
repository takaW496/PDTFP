# for multi_padlock_design
# Xiaoyan, 2017

import os
import re
import subprocess
import shutil

Processes = []
NextProcess = 0
input_files = 0


def newblast(file_input, target_genome, file_output):
    """ Start a new blastn subprocess if there is work to do """
    global Processes
    global NextProcess

    if NextProcess < len(input_files):
        if not os.path.isfile(file_output):     # skip already existing files
            txtcmd = ' '.join(('blastn', '-query', '"' + file_input + '"',
                               '-db', '"' + target_genome + '"',
                               '-outfmt 10',
                               '-out ', '"' + file_output + '"',
                               '-word_size 7 -strand plus'))

            blastnprocess = subprocess.Popen(txtcmd, shell=True)
            Processes.append(blastnprocess)
        NextProcess += 1


def runningblast(file_input, target_genome, file_output):
    """ Check any running processes and start new ones if there are spare slots """
    global Processes
    global NextProcess

    for p in range(len(Processes)-1, -1, -1):     # check the processes in reverse order
        if Processes[p].poll() is not None:     # if the process hasn't finished will return None
            del Processes[p]
    while (len(Processes) < 4) and (NextProcess < len(input_files)):  # more to do and some spare slots
        newblast(file_input, target_genome, file_output)


def blast_bridge_seq(df_Bp, target_genome, dir_blast_in, dir_blast_out):
    global Processes
    global NextProcess
    global input_files


    dir_blast_in_bridge = os.path.join(dir_blast_in, 'Bridge_Anchor', os.path.basename(target_genome).split('.')[0])
    dir_blast_out_bridge = os.path.join(dir_blast_out, 'Bridge_Anchor', os.path.basename(target_genome).split('.')[0])

    os.makedirs(dir_blast_in_bridge, exist_ok=True)
    os.makedirs(dir_blast_out_bridge, exist_ok=True)

    def check_previous_result(dir_blast_out_bridge):
        files = os.listdir(dir_blast_out_bridge)
        files = [f for f in files if not f.startswith(".")]
        if not files:
            return True
        else:
            return False

    if not check_previous_result(dir_blast_out_bridge):
        print('previous blast run for this genome found. use previous results to rank bridge seq. skip check bridge seq step.')
        return
    else:
        print('bridge sequences are not evaluated against this genome. start blast step and rank bridge sequences.')

    for row in df_Bp.itertuples():
        with open(os.path.join(dir_blast_in_bridge, row.BridgeProbe_ID + '.fasta'), 'w') as fblast:
            fblast.write(">%s_PLP_bridge\n%s\n>%s_PLP_anchor\n%s\n" % (row.BridgeProbe_ID, row.Seq_Barcode, row.BridgeProbe_ID, row.Seq_Anchor))

    print ("Starting blast..")
    input_files = os.listdir(dir_blast_in_bridge)
    
    for f in input_files:
        NextProcess = 0
        file_input = os.path.join(dir_blast_in_bridge, f)
        out_file_name = re.sub('.fasta', '_blast.txt', f)
        file_output = os.path.join(dir_blast_out_bridge, out_file_name)

        runningblast(file_input, target_genome, file_output)  # start the max processes running
        while len(Processes) > 0:   # still going on
            runningblast(file_input, target_genome, file_output)
    print ("Blast finished!")

    return

def blastdb(target_genome):
    """ Format fasta sequences to BLAST database """
    if not os.path.isfile(target_genome + '.nin'):
        try:
            # make database file from fna
            txtcmd = ' '.join(('makeblastdb -in', target_genome, '-dbtype nucl'))
            os.system(txtcmd)

        except:
            print(" Could not format BLAST database")

def continueblast(target_gene_IDs, target_genome, dir_blast_in, dir_blast_out):
    global Processes
    global NextProcess
    global input_files

    print ("Starting blast using candidate sequences as queries...")

    for gene in target_gene_IDs:
        print(f'start blast for {gene}')
        input_files = os.listdir(os.path.join(dir_blast_in, gene))
        if os.path.exists(os.path.join(dir_blast_out, gene)):
            shutil.rmtree(os.path.join(dir_blast_out, gene))
        os.makedirs(os.path.join(dir_blast_out, gene), exist_ok=True)
        for f in input_files:
            NextProcess = 0
            file_input = os.path.join(dir_blast_in, gene, f)
            out_file_name = re.sub('.fasta', '_blast.txt', f)
            file_output = os.path.join(dir_blast_out, gene, out_file_name)

            runningblast(file_input, target_genome, file_output)  # start the max processes running
            while len(Processes) > 0:   # still going on
                runningblast(file_input, target_genome, file_output)
    print ("Blast finished!")