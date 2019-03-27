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
grep "^chr21" gencode.v29.annotation.gtf > gencode.v29.chr21.annotation.gtf

# rename chrom from 21 to chr21 in vcf
sed -i 's/^21/chr21/' ALL.chr21_GRCh38_sites.20170504.vcf

# subset reads to those mapping to chr21
# need to align SRR307005 to chr21 outside of the pipeline
samtools view -h SRR307005.bam | awk '{if ($3 == "chr21" || substr($0,1,1) == "@"){print $0}}' | samtools view -b > SRR307005.chr21.bam

# convert read alignments to fastq
samtools fastq -1 SRR307005.chr21.1.fq -2 SRR307005.chr21.2.fq -0 /dev/null -s /dev/null -n -F 0x900 SRR307005.chr21.bam

# construct graph from variants and reference sequence, using the docker
docker run --mount type=bind,src=`pwd`,dst=/wd quay.io/vgteam/vg:ci-249-439f22abc8f915513b8c9365cd146481c97ac842 /vg/bin/vg construct -t 32 -r /wd/GRCh38_chr21.fa -v /wd/ALL.chr21_GRCh38_sites.20170504.vcf > GRCh38_chr21.vg

# construct splice graph
docker run --mount type=bind,src=`pwd`,dst=/wd quay.io/vgteam/vg:ci-249-439f22abc8f915513b8c9365cd146481c97ac842 /vg/bin/vg rna -n /wd/gencode.v29.chr21.annotation.gtf -e -t 32 -p /wd/GRCh38_chr21.vg > GRCh38_chr21.splice.vg

# index the splice graph
docker run --mount type=bind,src=`pwd`,dst=/wd quay.io/vgteam/vg:ci-249-439f22abc8f915513b8c9365cd146481c97ac842 /vg/bin/vg index -p -x /wd/GRCh38_chr21.splice.xg /wd/GRCh38_chr21.splice.vg
docker run --mount type=bind,src=`pwd`,dst=/wd quay.io/vgteam/vg:ci-249-439f22abc8f915513b8c9365cd146481c97ac842 /vg/bin/vg index -p -g /wd/GRCh38_chr21.splice.gcsa /wd/GRCh38_chr21.splice.vg

# map reads
vg map

# calculate per base-pair coverage
vg pack

# sum over genes
python3 gene_quant.py

