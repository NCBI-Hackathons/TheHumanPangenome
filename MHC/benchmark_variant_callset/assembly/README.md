# Assembly and Dipcall output

## Code
The dipcall code was retrieved from https://github.com/lh3/dipcall/releases version v0.1 with GitHub commit 7746f33. We modify the https://github.com/lh3/dipcall/blob/master/run-dipcall#L40 by adding -z400,200 to increase mappability around MHC region

## Operation

Here is a series of commands that were performed. 

wget -O GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz
dipcall.kit/samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
/home/dnanexus/dipcall.kit/run-dipcall grch38_hg002 /home/dnanexus/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa /home/dnanexus/hap1_fasta.gz /home/dnanexus/hap2_fasta.gz -x /home/dnanexus/dipcall.kit/hs38.PAR.bed
make -j2 -f grch38_hg002.mak


#The command `make -j2 -f hg002_denovo_grch38.mak` performs a series of operation under the hood.

/home/dnanexus/dipcall.kit/minimap2 -c --paf-no-hit -xasm5 -z400,200 --cs -r2k -t8 /home/dnanexus/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa /home/dnanexus/hap1_fasta.gz 2> grch38_hg002.hap1.paf.gz.log | gzip > grch38_hg002.hap1.paf.gz
/home/dnanexus/dipcall.kit/minimap2 -c --paf-no-hit -xasm5 -z400,200 --cs -r2k -t8 /home/dnanexus/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa /home/dnanexus/hap2_fasta.gz 2> grch38_hg002.hap2.paf.gz.log | gzip > grch38_hg002.hap2.paf.gz
/home/dnanexus/dipcall.kit/minimap2 -a -xasm5 -z400,200 --cs -r2k -t8 /home/dnanexus/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa /home/dnanexus/hap1_fasta.gz 2> grch38_hg002.hap1.sam.gz.log | gzip > grch38_hg002.hap1.sam.gz
/home/dnanexus/dipcall.kit/minimap2 -a -xasm5 -z400,200 --cs -r2k -t8 /home/dnanexus/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa /home/dnanexus/hap2_fasta.gz 2> grch38_hg002.hap2.sam.gz.log | gzip > grch38_hg002.hap2.sam.gz
gzip -dc grch38_hg002.hap1.paf.gz | sort -k6,6 -k8,8n | /home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/paftools.js call - 2> grch38_hg002.hap1.var.gz.vst | gzip > grch38_hg002.hap1.var.gz
gzip -dc grch38_hg002.hap2.paf.gz | sort -k6,6 -k8,8n | /home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/paftools.js call - 2> grch38_hg002.hap2.var.gz.vst | gzip > grch38_hg002.hap2.var.gz
/home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/dipcall-aux.js samflt grch38_hg002.hap1.sam.gz | /home/dnanexus/dipcall.kit/samtools sort -m4G --threads 4 -o grch38_hg002.hap1.bam -
gzip -dc grch38_hg002.hap1.var.gz | grep ^R | cut -f2- > grch38_hg002.hap1.bed
gzip -dc grch38_hg002.hap2.var.gz | grep ^R | cut -f2- > grch38_hg002.hap2.bed
/home/dnanexus/dipcall.kit/bedtk isec -m grch38_hg002.hap1.bed grch38_hg002.hap2.bed > grch38_hg002.dip.bed
/home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/dipcall-aux.js samflt grch38_hg002.hap2.sam.gz | /home/dnanexus/dipcall.kit/samtools sort -m4G --threads 4 -o grch38_hg002.hap2.bam -
/home/dnanexus/dipcall.kit/htsbox pileup -q5 -evcf /home/dnanexus/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa grch38_hg002.hap1.bam grch38_hg002.hap2.bam | /home/dnanexus/dipcall.kit/htsbox bgzip > grch38_hg002.pair.vcf.gz
/home/dnanexus/dipcall.kit/k8 /home/dnanexus/dipcall.kit/dipcall-aux.js vcfpair grch38_hg002.pair.vcf.gz | /home/dnanexus/dipcall.kit/htsbox bgzip > grch38_hg002.dip.vcf.gz

#For hg19 the command is similar except the ref were taken from

wget -O hs37d5.fa.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

and the PAR file is different dipcall.kit/hs37d5.PAR.bed


