#!/bin/bash

# download the reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
gzip -d GCA_000001405.28_GRCh38.p13_genomic.fna.gz

# download the reference annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gzip -d gencode.v29.annotation.gtf.gz

# download the 1000G variants for chromosome 21, corrected to GRCh28 coordinates
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr21_GRCh38_sites.20170504.vcf.gz
gzip -d ALL.chr21_GRCh38_sites.20170504.vcf.gz

# subset reference genome to chromosome 21
samtools faidx -o GRCh38_chr21.fa GCA_000001405.28_GRCh38.p13_genomic.fna "CM000683.2"

# rename sequence to 'chr21'
sed -i 's/^>CM000683.2/>chr21/' GRCh38_chr21.fa

# subset reference annotation to chromosome 21
awk '{if ($3 == "chr21"){print $0}}' gencode.v29.annotation.gtf > gencode.v29.chr21.gtf

# rename chrom from 21 to chr21 in vcf
sed -i 's/^21/chr21/' ALL.chr21_GRCh38_sites.20170504.vcf

# subset reads to those mapping to chr21
# need to align SRA307005 to chr21 outside of the pipeline
samtools view -h SRA307005.bam | awk 'if ($3 == "chr21" || substr($0,1,1) == "@"){print $0}' | samtools view -b > SRA307005.chr21.bam

# convert read alignments to fastq
samtools fastq

# construct graph from variants and reference sequence
vg construct

# construct splice graph
vg rna

# map reads
vg map

# calculate per base-pair coverage
vg pack

# sum over genes
python3 gene_quant.py

