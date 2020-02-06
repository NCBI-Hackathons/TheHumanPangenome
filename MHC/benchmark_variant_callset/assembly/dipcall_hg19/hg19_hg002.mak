ROOT=/home/dnanexus/dipcall.kit
N_THREADS=8
MM2_IDX=/home/dnanexus/hs37d5.fa
MM2_OPT=-xasm5 -z400,200 --cs -r2k -t$(N_THREADS)
REF_FA=/home/dnanexus/hs37d5.fa

all:hg19_hg002.dip.bed hg19_hg002.dip.vcf.gz

hg19_hg002.hap1.paf.gz:/home/dnanexus/hap1_fasta.gz
	$(ROOT)/minimap2 -c --paf-no-hit $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@
hg19_hg002.hap2.paf.gz:/home/dnanexus/hap2_fasta.gz
	$(ROOT)/minimap2 -c --paf-no-hit $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@

hg19_hg002.hap1.sam.gz:/home/dnanexus/hap1_fasta.gz
	$(ROOT)/minimap2 -a $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@
hg19_hg002.hap2.sam.gz:/home/dnanexus/hap2_fasta.gz
	$(ROOT)/minimap2 -a $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@

hg19_hg002.hap1.bam:hg19_hg002.hap1.sam.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js samflt $< | $(ROOT)/samtools sort -m4G --threads 4 -o $@ -
hg19_hg002.hap2.bam:hg19_hg002.hap2.sam.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js samflt $< | $(ROOT)/samtools sort -m4G --threads 4 -o $@ -

hg19_hg002.pair.vcf.gz:hg19_hg002.hap1.bam hg19_hg002.hap2.bam
	$(ROOT)/htsbox pileup -q5 -evcf $(REF_FA) $^ | $(ROOT)/htsbox bgzip > $@
hg19_hg002.dip.vcf.gz:hg19_hg002.pair.vcf.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js vcfpair $< | $(ROOT)/htsbox bgzip > $@

hg19_hg002.hap1.var.gz:hg19_hg002.hap1.paf.gz
	gzip -dc $< | sort -k6,6 -k8,8n | $(ROOT)/k8 $(ROOT)/paftools.js call - 2> $@.vst | gzip > $@
hg19_hg002.hap2.var.gz:hg19_hg002.hap2.paf.gz
	gzip -dc $< | sort -k6,6 -k8,8n | $(ROOT)/k8 $(ROOT)/paftools.js call - 2> $@.vst | gzip > $@
hg19_hg002.hap1.bed:hg19_hg002.hap1.var.gz
	gzip -dc $< | grep ^R | cut -f2- > $@
hg19_hg002.hap2.bed:hg19_hg002.hap2.var.gz
	gzip -dc $< | grep ^R | cut -f2- > $@
hg19_hg002.dip.bed:hg19_hg002.hap1.bed hg19_hg002.hap2.bed
	$(ROOT)/bedtk isec -m $^ > $@
