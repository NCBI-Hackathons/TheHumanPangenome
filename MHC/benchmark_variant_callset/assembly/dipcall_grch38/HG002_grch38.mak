ROOT=/home/dnanexus/dipcall.kit
N_THREADS=8
MM2_IDX=/home/dnanexus/hg38.fa
MM2_OPT=-xasm5 --cs -r2k -t$(N_THREADS)
REF_FA=/home/dnanexus/hg38.fa

all:HG002_grch38.dip.bed HG002_grch38.dip.vcf.gz

HG002_grch38.hap1.paf.gz:/home/dnanexus/hap1_fasta.gz
	$(ROOT)/minimap2 -c --paf-no-hit $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@
HG002_grch38.hap2.paf.gz:/home/dnanexus/hap2_fasta.gz
	$(ROOT)/minimap2 -c --paf-no-hit $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@

HG002_grch38.hap1.sam.gz:/home/dnanexus/hap1_fasta.gz
	$(ROOT)/minimap2 -a $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@
HG002_grch38.hap2.sam.gz:/home/dnanexus/hap2_fasta.gz
	$(ROOT)/minimap2 -a $(MM2_OPT) $(MM2_IDX) $< 2> $@.log | gzip > $@

HG002_grch38.hap1.bam:HG002_grch38.hap1.sam.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js samflt $< | $(ROOT)/samtools sort -m4G --threads 4 -o $@ -
HG002_grch38.hap2.bam:HG002_grch38.hap2.sam.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js samflt $< | $(ROOT)/samtools sort -m4G --threads 4 -o $@ -

HG002_grch38.pair.vcf.gz:HG002_grch38.hap1.bam HG002_grch38.hap2.bam
	$(ROOT)/htsbox pileup -q5 -evcf $(REF_FA) $^ | $(ROOT)/htsbox bgzip > $@
HG002_grch38.dip.vcf.gz:HG002_grch38.pair.vcf.gz
	$(ROOT)/k8 $(ROOT)/dipcall-aux.js vcfpair $< | $(ROOT)/htsbox bgzip > $@

HG002_grch38.hap1.var.gz:HG002_grch38.hap1.paf.gz
	gzip -dc $< | sort -k6,6 -k8,8n | $(ROOT)/k8 $(ROOT)/paftools.js call - 2> $@.vst | gzip > $@
HG002_grch38.hap2.var.gz:HG002_grch38.hap2.paf.gz
	gzip -dc $< | sort -k6,6 -k8,8n | $(ROOT)/k8 $(ROOT)/paftools.js call - 2> $@.vst | gzip > $@
HG002_grch38.hap1.bed:HG002_grch38.hap1.var.gz
	gzip -dc $< | grep ^R | cut -f2- > $@
HG002_grch38.hap2.bed:HG002_grch38.hap2.var.gz
	gzip -dc $< | grep ^R | cut -f2- > $@
HG002_grch38.dip.bed:HG002_grch38.hap1.bed HG002_grch38.hap2.bed
	$(ROOT)/bedtk isec -m $^ > $@
