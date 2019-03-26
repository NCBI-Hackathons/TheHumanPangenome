## Original VCF used to make the graph (SVs for 15 individuals)
VCF_GRAPH=graph.explicit.vcf.gz

## SV calls from different dataset
VCF_SR=illumina/vg-call.vcf.gz
VCF_SRRD=illumina/vg-call.vcf.gz
VCF_CCS=XX/vg-call.vcf.gz
VCF_LR=XX/vg-call.vcf.gz

REF=xxx.fa

## Normalize original VCF
bcftools norm -m -both -f $REF $VCF_GRAPH > ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz

## Short-reads
bcftools norm -m -both -f $REF $VCF_SR > ${VCF_SR%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_SR%.vcf.gz}.norm.vcf.gz variantsInInsertions.SR.tsv

## Short-reads regional downsample
bcftools norm -m -both -f $REF $VCF_SRRD > ${VCF_SRRD%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_SRRD%.vcf.gz}.norm.vcf.gz variantsInInsertions.SRRD.tsv

## CCS
bcftools norm -m -both -f $REF $VCF_CCS > ${VCF_CCS%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_CCS%.vcf.gz}.norm.vcf.gz variantsInInsertions.CCS.tsv

## Long-reads
bcftools norm -m -both -f $REF $VCF_LR > ${VCF_LR%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_LR%.vcf.gz}.norm.vcf.gz variantsInInsertions.LR.tsv

## COMPARE variantsInInsertions.*.tsv FILES (Venn diagram, bar plot)
## Both (SR vs SRRD) and (SR vs LR vs CCS)
