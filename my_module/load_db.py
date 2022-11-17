import re

def readgenefile(round_table):
    """ Get the gene list """
    try:
        with open(round_table, 'r') as f:
            genes = [line.split('\t')[0] for line in f]
            genes = [s for s in genes if not s.startswith('gene_ID')]
            genes.sort()
    except IOError:
        print("Gene list file not found. Try again.")
        genes = []
    return genes

def load_fastadb(target_genome):
    """ Read multiple fasta files of a database and make dictionary {gene_id : sequence} """
    # read original fasta files
    d = {}
    header = str()
    seq = str()

    with open(target_genome, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line[0] == '>':      # a new sequence
                d[header] = seq 
                header = line.replace('>', '')
                seq = str()
            else:
                seq += line
                
    d[header] = seq     # the last sequence
    del d['']      # delete first empty element due to appending only
    return d

def findseq(genes, d):
    """ Find target sequences """
    target_seqs = []

    for gene in genes:
        pattern = re.compile(gene)
        pattern_hit = [d[k] for k in d if pattern.search(k)]
        if len(pattern_hit) == 1:
            target_seqs.append(pattern_hit[0])
        elif len(pattern_hit) > 1:
            print(f'multiple entries found with {gene} in the target genome. use the shortest sequence.')
            target_seqs.append(min(pattern_hit)[0])
        else:
            print(f'No hit found with this gene ID: {gene}')

    return target_seqs