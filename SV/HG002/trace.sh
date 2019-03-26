## Original VCF used to make the graph (SVs for 15 individuals)
VCF_GRAPH=/data/vg/sv-pop-explicit.vcf.gz

## SV calls from different dataset
VCF_SR=/data/vg/HG00514.vcf.gz ## temporary (for testing on a vg call output)
VCF_SRRD=illumina/vg-call-regional-downsample.vcf.gz
VCF_CCS=XX/vg-call.vcf.gz
VCF_LR=XX/vg-call.vcf.gz

REF=/data/GRCh38/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
REGIONS=/data/SVs/SVSummaryPop_EUR-AFR.extract-regions.bed

## Normalize original VCF
bcftools norm -m -both -f $REF $VCF_GRAPH | bgzip -c > ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz

## Short-reads
bcftools norm -m -both -f $REF $VCF_SR | bgzip -c > ${VCF_SR%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_SR%.vcf.gz}.norm.vcf.gz variantsInInsertions.SR.tsv $REGIONS

## Short-reads regional downsample
bcftools norm -m -both -f $REF $VCF_SRRD | bgzip -c > ${VCF_SRRD%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_SRRD%.vcf.gz}.norm.vcf.gz variantsInInsertions.SRRD.tsv $REGIONS

## CCS
bcftools norm -m -both -f $REF $VCF_CCS | bgzip -c > ${VCF_CCS%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_CCS%.vcf.gz}.norm.vcf.gz variantsInInsertions.CCS.tsv $REGIONS

## Long-reads
bcftools norm -m -both -f $REF $VCF_LR | bgzip -c > ${VCF_LR%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_LR%.vcf.gz}.norm.vcf.gz variantsInInsertions.LR.tsv $REGIONS

## COMPARE variantsInInsertions.*.tsv FILES (Venn diagram, bar plot)
## Both (SR vs SRRD) and (SR vs LR vs CCS)
