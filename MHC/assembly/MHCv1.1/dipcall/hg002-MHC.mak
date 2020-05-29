ROOT=/home/dnanexus/dipcall.kit
N_THREADS=8
MM2_IDX=/home/dnanexus/hs37d5.fa
MM2_OPT=-xasm5 -z200000,10000 --cs -r2k -t$(N_THREADS)
REF_FA=/home/dnanexus/hs37d5.fa

all:hg002-MHC.dip.bed hg002-MHC.dip.vcf.gz

hg002-MHC.hap1.paf.gz:/home/dnanexus/hap1_fasta.gz
	$(ROOT)/minimap2 -c --paf-no-hit $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@
hg002-MHC.hap2.paf.gz:/home/dnanexus/hap2_fasta.gz
	$(ROOT)/minimap2 -c --paf-no-hit $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@

hg002-MHC.hap1.sam.gz:/home/dnanexus/hap1_fasta.gz
	$(ROOT)/minimap2 -a $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@
hg002-MHC.hap2.sam.gz:/home/dnanexus/hap2_fasta.gz
	$(ROOT)/minimap2 -a $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@

hg002-MHC.hap1.bam:hg002-MHC.hap1.sam.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js samflt $< | $(ROOT)/samtools sort -m4G --threads 4 -o $@ -
hg002-MHC.hap2.bam:hg002-MHC.hap2.sam.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js samflt $< | $(ROOT)/samtools sort -m4G --threads 4 -o $@ -

hg002-MHC.pair.vcf.gz:hg002-MHC.hap1.bam hg002-MHC.hap2.bam
	$(ROOT)/htsbox pileup -q5 -evcf $(REF_FA) $^ | $(ROOT)/htsbox bgzip > $@
hg002-MHC.dip.vcf.gz:hg002-MHC.pair.vcf.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js vcfpair $< | $(ROOT)/htsbox bgzip > $@

hg002-MHC.hap1.var.gz:hg002-MHC.hap1.paf.gz
	gzip -dc $< | sort -k6,6 -k8,8n | $(ROOT)/k8 $(ROOT)/paftools.js call - 2> $@.vst | gzip > $@
hg002-MHC.hap2.var.gz:hg002-MHC.hap2.paf.gz
	gzip -dc $< | sort -k6,6 -k8,8n | $(ROOT)/k8 $(ROOT)/paftools.js call - 2> $@.vst | gzip > $@
hg002-MHC.hap1.bed:hg002-MHC.hap1.var.gz
	gzip -dc $< | grep ^R | cut -f2- > $@
hg002-MHC.hap2.bed:hg002-MHC.hap2.var.gz
	gzip -dc $< | grep ^R | cut -f2- > $@
hg002-MHC.dip.bed:hg002-MHC.hap1.bed hg002-MHC.hap2.bed
	$(ROOT)/bedtk isec -m $^ > $@
