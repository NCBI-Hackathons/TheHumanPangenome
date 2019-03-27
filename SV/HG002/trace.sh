## Original VCF used to make the graph (SVs for 15 individuals)
VCF_GRAPH=/data/vg/sv-pop-explicit.vcf.gz

## SV calls from different dataset
VCF_SR=/data/vg/HG002-30x-vg.vcf.gz
VCF_SRRD=/data/vg/HG002-sv-regions-vg.vcf.gz
VCF_CCS=XX/vg-call.vcf.gz
VCF_LR=XX/vg-call.vcf.gz

REF=/data/GRCh38/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
REFG=/data/GRCh38/hg38.fa
REGIONS=/data/SVs/SVSummaryPop_EUR-AFR.extract-regions.bed

## Normalize original VCF
bcftools norm -m -both -f $REF $VCF_GRAPH | bgzip -c > ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz

## Short-reads
bcftools norm -m -both -f $REFG $VCF_SR | bgzip -c > ${VCF_SR%.vcf.gz}.norm.vcf.gz
zcat ${VCF_SR%.vcf.gz}.norm.vcf.gz | awk '{if(substr($0, 1, 1)=="#" || (length($4)==1 && length($5)>20)){print$0}}' | bgzip -c > ${VCF_SR%.vcf.gz}.norm.nosnp.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_SR%.vcf.gz}.norm.nosnp.vcf.gz /data/variants_in_insertions/HG002/variantsInInsertions.SR.norm $REGIONS
Rscript filterInternalVariants.R /data/variants_in_insertions/HG002/variantsInInsertions.SR.norm.tsv /data/variants_in_insertions/HG002/variantsInInsertions.SR.norm.filtered.tsv

## Short-reads regional downsample
bcftools norm -m -both -f $REFG $VCF_SRRD | bgzip -c > ${VCF_SRRD%.vcf.gz}.norm.vcf.gz
zcat ${VCF_SRRD%.vcf.gz}.norm.vcf.gz | awk '{if(substr($0, 1, 1)=="#" || (length($4)==1 && length($5)>20)){print$0}}' | bgzip -c > ${VCF_SRRD%.vcf.gz}.norm.ins.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_SRRD%.vcf.gz}.norm.ins.vcf.gz /data/variants_in_insertions/HG002/variantsInInsertions.SRRD $REGIONS

## Align insertions SR vs SRRD
Rscript alignCalledInsertions.R ${VCF_SR%.vcf.gz}.norm.nosnp.vcf.gz ${VCF_SRRD%.vcf.gz}.norm.ins.vcf.gz /data/variants_in_insertions/HG002/align.SR.SRRD.tsv $REGIONS

## CCS
bcftools norm -m -both -f $REF $VCF_CCS | bgzip -c > ${VCF_CCS%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_CCS%.vcf.gz}.norm.vcf.gz /data/variants_in_insertions/HG002/variantsInInsertions.CCS $REGIONS

## Long-reads
bcftools norm -m -both -f $REF $VCF_LR | bgzip -c > ${VCF_LR%.vcf.gz}.norm.vcf.gz
Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF_LR%.vcf.gz}.norm.vcf.gz /data/variants_in_insertions/HG002/variantsInInsertions.LR $REGIONS

## COMPARE variantsInInsertions.*.tsv FILES (Venn diagram, bar plot)
## Both (SR vs SRRD) and (SR vs LR vs CCS)
