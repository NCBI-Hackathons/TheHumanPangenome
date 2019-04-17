#!/bin/bash

# Error correction of ONT reads using CCS reads
# Output will be in {H1/H2/untagged}/corrected_ONT_{H1/H2/untagged}.fa and {H1/H2/untagged}/corrected_ONT_{H1/H2/untagged}_split.fa

# The script assumes that the CCS and ONT reads are split by haplotype.
# To split the reads by haplotype see: https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/MHC/phasing-notes.md
CCS_H1=/data/newsplit/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.newsplit.28498559.H1.fastq.gz
CCS_H2=/data/newsplit/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.newsplit.28498559.H2.fastq.gz
CCS_UNTAGGED=/data/newsplit/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.newsplit.28498559.untagged.fastq.gz
ONT_H1=/data/newsplit/HG002_Promethion.MHConly.28498559.H1.fastq.gz
ONT_H2=/data/newsplit/HG002_Promethion.MHConly.28498559.H2.fastq.gz
ONT_UNTAGGED=/data/newsplit/HG002_Promethion.MHConly.28498559.untagged.fastq.gz

# https://github.com/maickrau/GraphAligner
# The extract scripts are NOT included in the bioconda release. You have to compile them from source. The aligner can be either compiled or from bioconda
GRAPHALIGNER_PATH=GraphAligner
EXTRACT_PATH=./ExtractCorrectedReads
EXTRACT_SPLIT_PATH=./ExtractPathSequence

# https://github.com/GATB/bcalm
# Or bioconda: https://anaconda.org/bioconda/bcalm/files
BCALM_PATH=bcalm

# https://github.com/GATB/bcalm/blob/master/scripts/convertToGFA.py
BCALM_CONVERT_PATH=./convertToGFA.py

#############################

# correct the first haplotype

mkdir H1
cd H1
# include CCS reads from the first haplotype and the untagged CCS
ls -1 $CCS_H1 $CCS_UNTAGGED > files
# build DBG
/usr/bin/time -v $BCALM_PATH -in files -kmer-size 63 -abundance-min 3 -nb-cores 40 1> bcalm_stdout.txt 2> bcalm_stderr.txt
# remove temporary files
rm files files.unitigs.fa.* files.h5
# convert the DBG to a GFA
$BCALM_CONVERT_PATH files.unitigs.fa MHC-H1-k63-a3.gfa 63
# remove temporary file
rm files.unitigs.fa
# align the ONT reads from the first haplotype
/usr/bin/time -v $GRAPHALIGNER_PATH -g MHC-H1-k63-a3.gfa -f $ONT_H1 -a alns.gam -b 50 -t 40 1> aligner_stdout.txt 2> aligner_stderr.txt
# replace the aligned read sequence with the graph sequence, keep unaligned sequence in lowercase
$EXTRACT_PATH MHC-H1-k63-a3.gfa alns.gam $ONT_H1 > corrected_ONT_H1.fa 2> /dev/null
# extract the paths of the alignments, unaligned sequence is discarded and splits the read
$EXTRACT_SPLIT_PATH MHC-H1-k63-a3.gfa alns.gam > corrected_ONT_H1_split.fa 2> /dev/null

# do the same for haplotype 2

cd ..
mkdir H2
cd H2
ls -1 $CCS_H2 $CCS_UNTAGGED > files
/usr/bin/time -v $BCALM_PATH -in files -kmer-size 63 -abundance-min 3 -nb-cores 40 1> bcalm_stdout.txt 2> bcalm_stderr.txt
rm files files.unitigs.fa.* files.h5
$BCALM_CONVERT_PATH files.unitigs.fa MHC-H2-k63-a3.gfa 63
rm files.unitigs.fa
/usr/bin/time -v $GRAPHALIGNER_PATH -g MHC-H2-k63-a3.gfa -f $ONT_H2 -a alns.gam -b 50 -t 40 1> aligner_stdout.txt 2> aligner_stderr.txt
$EXTRACT_PATH MHC-H2-k63-a3.gfa alns.gam $ONT_H2 > corrected_ONT_H2.fa 2> /dev/null
$EXTRACT_SPLIT_PATH MHC-H2-k63-a3.gfa alns.gam > corrected_ONT_H2_split.fa 2> /dev/null

# for untagged reads, build the graph using both haplotypes and untagged

cd ..
mkdir untagged
cd untagged
ls -1 $CCS_H1 $CCS_H2 $CCS_UNTAGGED > files
/usr/bin/time -v $BCALM_PATH -in files -kmer-size 63 -abundance-min 3 -nb-cores 40 1> bcalm_stdout.txt 2> bcalm_stderr.txt
rm files files.unitigs.fa.* files.h5
$BCALM_CONVERT_PATH files.unitigs.fa MHC-untagged-k63-a3.gfa 63
rm files.unitigs.fa
/usr/bin/time -v $GRAPHALIGNER_PATH -g MHC-untagged-k63-a3.gfa -f $ONT_UNTAGGED -a alns.gam -b 50 -t 40 1> aligner_stdout.txt 2> aligner_stderr.txt
$EXTRACT_PATH MHC-untagged-k63-a3.gfa alns.gam $ONT_UNTAGGED > corrected_ONT_untagged.fa 2> /dev/null
$EXTRACT_SPLIT_PATH MHC-untagged-k63-a3.gfa alns.gam > corrected_ONT_untagged_split.fa 2> /dev/null
cd ..
