# PDTFP - Probe Design Tool For PHYTOMap

This is the tool to design probes for PHYTOMap. Currently testing whether probes designed by this tool actually works fine for PHYTOMap.

- PHYTOMap paper Nobori et al., 2022 https://www.biorxiv.org/content/10.1101/2022.07.28.501915v1.full
- Python scripts modified from https://github.com/Moldia/multi_padlock_design (HYBISS paper, Gyllborg, D. et al., 2020)


## Quick start

- prepare two input files: `round_table.txt` and a cds fasta file of your target species. Please see example files.
- configure parameteres and paths to the required files in the second block of the `probes_design.py`.
- run `probe_desing.py`, and check the output file `candidate.tsv` and choose your favourite sequences.


## Dependency

- Python
- BLAST
- Bio.seq


## Workflow

- read target gene ids from the 1st column of the `round_table.txt` (`load_db.readgenefile`)
- load a cds fasta of the target genome (`load_db.load_fastadb`)
- extract target gene sequences from a cds fasta of the target genome (`load_db.findseq`)
- search seed sequences for each gene (`screenseq.search_seed_seq`). Seed sequences (of defined length by `seed_seq_len`) will be selected based on the GC content and the melting temperatures. seed sequnces will be stored in `blast_in` for the blast search in the nexxt step. Return a dataframe listing candidate seed sequences.
- run blast search against the target genome using seed sequences in `blast_in` as queries (`parblast.continueblast`).  blast results will be stored in `blast_out`.
- assess blast results in `blast_out` and judge whether seed sequences are specific enough (`checkblast.getcandidates`). add `is_specific` column to the dataframe
- group seed sequences located closely (`checkblast.add_id_clust`). add `id_clust` column to the dataframe.
- generate padlock probes and primers based on seed sequences, then corresponding bridge probes will be assigned to each gene based on the detection dye specified in the round table (`finalize_probe.add_probe_seq`).
- manually check suggested probe sequences and determine which sequences you will use for the experiment. choose one seed sequence from one cluster to avoid large overlaps between probes. If there are not enough sequences for certain genes, change round tables and parameters, and run program again.

## Input parameters

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

## Input file format
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
### probe design principal

- PLP
/5Phos/ACATTA + rev_compl{target_sequence1} + {unique_bridge_seq_from_HybISS} + {anchor} + AAGATA

- primer
rev_compl{target_sequence2} + TAATGTTATCTT

- bridge probe
{unique_bridge_seq} + (C when probe is AF750) + rev_compl{fluorescent_probe1-4}
unique bridge seq is a bit modified compared to the one in PLP. one nucleotide at the 5' end is trimmed off and two nucleotides at the 3' end are changed to complemented nuleotides.


### seed sequence search

- seed sequences are successively generated from the start codon of the gene sequence with the window size of the `seed_seq_len` and the step size of 5 bp. 
- if there are previous results in blast_in folder, those edirectories will be deleted.
- gc content should be gcmin < gc < gcmax.
- seed sequences will be divided into 2 parts with 2 bp intervals inbetween. For example, if `seed_seq_len` is 50, 1st 24 bp, interval 2 bp, 2nd 24 bp. When `seed_seq_len` is odd number, the second half is 1 bp-longer than the 1st half.
- calculate melting temperature of 1st and 2nd halves of the sequences and search the best balanced pair by changing the length of the sequence pair so that the difference between two melting temperatures is minimized.
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
