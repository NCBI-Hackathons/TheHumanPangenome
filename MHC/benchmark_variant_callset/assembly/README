# Assembly and Dipcall output

## Code
The dipcall code was retrieved from https://github.com/lh3/dipcall/releases version v0.1 with GitHub commit 7746f33

## Operation

Here is a series of command that were performed. 

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
dipcall.kit/samtools faidx hg38.fa
/home/dnanexus/dipcall.kit/run-dipcall hg002_denovo_grch38 /home/dnanexus/hg38.fa /home/dnanexus/hap1_fasta.gz /home/dnanexus/hap2_fasta.gz -x /home/dnanexus/dipcall.kit/hs38.PAR.bed
make -j2 -f hg002_denovo_grch38.mak


#The command `make -j2 -f hg002_denovo_grch38.mak` performs a series of operation under the hood.

/home/dnanexus/dipcall.kit/minimap2 -a -xasm5 --cs -r2k -t8 /home/dnanexus/hg38.fa /home/dnanexus/hap1_fasta.gz 2> hg002_denovo_grch38.hap1.sam.gz.log | gzip > hg002_denovo_grch38.hap1.sam.gz
/home/dnanexus/dipcall.kit/minimap2 -a -xasm5 --cs -r2k -t8 /home/dnanexus/hg38.fa /home/dnanexus/hap2_fasta.gz 2> hg002_denovo_grch38.hap2.sam.gz.log | gzip > hg002_denovo_grch38.hap2.sam.gz
gzip -dc hg002_denovo_grch38.hap1.paf.gz | sort -k6,6 -k8,8n | /home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/paftools.js call - 2> hg002_denovo_grch38.hap1.var.gz.vst | gzip > hg002_denovo_grch38.hap1.var.gz
gzip -dc hg002_denovo_grch38.hap2.paf.gz | sort -k6,6 -k8,8n | /home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/paftools.js call - 2> hg002_denovo_grch38.hap2.var.gz.vst | gzip > hg002_denovo_grch38.hap2.var.gz
/home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/dipcall-aux.js samflt hg002_denovo_grch38.hap1.sam.gz | /home/dnanexus/dipcall.kit/samtools sort -m4G --threads 4 -o hg002_denovo_grch38.hap1.bam -
gzip -dc hg002_denovo_grch38.hap1.var.gz | grep ^R | cut -f2- > hg002_denovo_grch38.hap1.bed
gzip -dc hg002_denovo_grch38.hap2.var.gz | grep ^R | cut -f2- > hg002_denovo_grch38.hap2.bed
/home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/dipcall-aux.js samflt hg002_denovo_grch38.hap2.sam.gz | /home/dnanexus/dipcall.kit/samtools sort -m4G --threads 4 -o hg002_denovo_grch38.hap2.bam -
/home/dnanexus/dipcall.kit/bedtk isec -m hg002_denovo_grch38.hap1.bed hg002_denovo_grch38.hap2.bed > hg002_denovo_grch38.dip.bed
/home/dnanexus/dipcall.kit/htsbox pileup -q5 -evcf /home/dnanexus/hg38.fa hg002_denovo_grch38.hap1.bam hg002_denovo_grch38.hap2.bam | /home/dnanexus/dipcall.kit/htsbox bgzip > hg002_denovo_grch38.pair.vcf.gz
/home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/dipcall-aux.js vcfpair hg002_denovo_grch38.pair.vcf.gz | /home/dnanexus/dipcall.kit/htsbox bgzip > hg002_denovo_grch38.dip.vcf.gz

#For hg19 the command is similar except the ref were taken from

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

and the PAR file is different dipcall.kit/hs37d5.PAR.bed


