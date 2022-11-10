# probe design tool for PHYTOMap

- PHYTOMap paper Nobori et al., 2022 https://www.biorxiv.org/content/10.1101/2022.07.28.501915v1.full
- Python scripts modified from https://github.com/Moldia/multi_padlock_design/tree/master/lib (HYBISS paper, Gyllborg, D. et al., 2020)


## Dependency

- BLAST
- Python
- Bio.seq


## WORKFLOW

- read target gene ids from round_table.txt (`load_db.readgenefile`)
- load a cds fasta (`load_db.load_fastadb`)
- read target gene sequences from a cds fasta (`load_db.findseq`)
- search seed sequences for each gene (`screenseq.threshold_gc_tm`). Seed sequences will be selected based on the GC content, the melting temperatures and seed sequence length. seed sequnces will be stored in `blast_in`.
- run blast search against the target genome using seed sequences as queries (`parblast.continueblast`).  blast results will be stored in `blast_out`.
- assess blast results and judge whether seed sequences are specific enough (`checkblast.getcandidates`). add `is_specific` column to the dataframe
- group seed sequences located closely (`checkblast.add_id_clust`). add `id_clust` column to the dataframe.
- design padlock probes and primers based on seed sequences, then corresponding bridge probes will be assigned to each gene based on the detection dye specified in round table (`finalize_probe.add_probe_seq`).
- manually check suggested probe sequences and determine which sequences you will use for the experiment. choose one seed sequence from one cluster to avoid large overlaps between probes. If there are not enough sequences, change parameters and run program again.s

## input parameters

- target_genome: path to the target genome
- round_table: path to the txt file in which your experimental scheme was written. plese see the input file format section.
- BridgeProbes: path to the txt file in which bridge sequences and anchor sequences are listed. These sequences are derived from the previous publication (HYBISS paper, Gyllborg, D. et al., 2020).
- dye_seq: path to the txt file in which dye and detection oligos are listed. These sequences are derived from the previous publication (HYBISS paper, Gyllborg, D. et al., 2020).
- seed_seq_len: integer, specify seed sequence length. 40-60 nt are usually used. defaul value is 50.
- gcmin: integer, minimum GC content (%) used in `screenseq.threshold_gc_tm`. default value is 40.
- gcmax: integer, maximum GC content (%) used in `screenseq.threshold_gc_tm`. defualt value is 60.
- tmmin: integer, minimum melting temperature (C) used in `screenseq.threshold_gc_tm`. default value is 48.
- tmmax: integer, maximum melting temperature (C) used in `screenseq.threshold_gc_tm`. default value is 52.
- outfile: path and name of the output txt file.
- threads: the number of threads you will use. this parameter is currently not used.

## input file format
### round_table.txt

- tab delimited text file with the header in the first row.
- target gene ID in the first column, detection round in the second column, detection dye in the third column. Forth coulmn onward will be ignored.
- gene ID should be exactly same as a reference cds fasta's header.

### reference_cds.fasta

- normal cds fasta file including all target gene sequences.
- when the gene ID (from round table) pattern hits multiple entries in cds fasta, the shortest sequences will be used as a representative.

### BridgeProbes.txt

- tab delimited text file with the header in the first row.
- BridgeProbe_ID in the first column, barcode sequences in the second column, anchor sequences in the third column.
- these sequences are derived from HYBISS paper and you do not need to modify this file unless you want to use or add other sequences.

### dye_seq.txt

- tab delimited text file with the header in the first row.
- dye name in the first column, detection nucleotide sequences in the second column
- these sequences are derived from HYBISS paper and you do not need to modify this file unless you want to use or add other dye and sequences.

## details
### seed sequence search

- seed sequences are successively generated from the start codon of the gene sequence with the window size of the `seed_seq_len` and the step size of 5 bp. 
- gc content should be gcmin < gc < gcmax.
- seed sequences will be divided into 2 parts with 2 bp intervals inbetween. For example, if `seed_seq_len` is 50, 1st 24 bp, interval 2 bp, 2nd 24 bp. When `seed_seq_len` is odd number, the second half is 1 bp-longer than the 1st half.
- calculate melting temperature of 1st and 2nd halves of the sequences and search the best balanced pair of the sequences by changing the position of the interval. if tm1 > tm2, 
- both sequence' melting temperature should be tmmin < tm < tmmax.
- seed sequences will be stored in `blast_in`.

### melting temperature calculation

- the melting temperature of the sequence will be calculated with the nearest neighbor model implemented in the original script (used in HYBISS paper).
- salt correction part and formamide correction part are modified to adjust to the hybridization buffer used in PHYTOMap protocol.
- Melting temperatures of nucleotide sequences from PHYTOMap paper are around 50C (explained as 60C in the paper) with this calculation method. therefore, default value is set to 50C.

### judging specificity of the sequence

- `checkblast.readblastout` will judge the specificity of the sequence
- sequences will be judged as non-specific if there are any blast hits with more than 50% coverage and 80% homology.

### grouping seed sequences

- `checkblast.add_id_clust` will group the seed sequnce depending on the distance in the gene sequence.
- seed sequences with any overlaps will be grouped in the same cluster.
- this classification depends on the `step` size of the seed sequence search and `seed_seq_len`.
- choose one seed sequence from one cluster to avoid large overlaps between probes.

### finalizing probes

- rank bridge sequences by running blast against target genome using bridge sequences as queries. Sequences will be ranked according to the highest bit score. Sequences with lower bit scores will be used first.
- the first half of the seed sequence will be used for padlock probes. barcode sequences, anchor sequences and misc sequences will be add. 
- the second half of the seed sequence will be used for primer sequences. one part of the padlock probes will be add.
- the bridge probe will be designed depending on the detection dye.


## Note

- 5' end of padlock probes should be phosphorylated. Be careful when you order oligos to the company.

