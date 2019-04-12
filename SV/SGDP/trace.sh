## Original VCF used to make the graph (SVs for 15 individuals)
VCF_GRAPH=/data/vg/sv-pop-explicit.vcf.gz

## Normalize original VCF
bcftools norm -m -both -f $REF $VCF_GRAPH | bgzip -c > ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz

REF=/data/GRCh38/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
REFG=/data/GRCh38/hg38.fa
REGIONS=/data/SVs/SVSummaryPop_EUR-AFR.extract-regions.bed


## SV calls from different dataset
for SAMP in  HG00128 HG02574
do
    echo $SAMP
    VCF=/data/vg/${SAMP}-30x-vg.vcf.gz
    bcftools norm -m -both -f $REFG $VCF | bgzip -c > temp
    sudo cp temp ${VCF%.vcf.gz}.norm.vcf.gz  # to avoid permission denied problems
    zcat ${VCF%.vcf.gz}.norm.vcf.gz | awk '{if(substr($0, 1, 1)=="#" || (length($4)==1 && length($5)>20)){print$0}}' | bgzip -c > temp
    sudo cp temp ${VCF%.vcf.gz}.norm.nosnp.vcf.gz
    Rscript callVariantsInInsertedSeq.R ${VCF_GRAPH%.vcf.gz}.norm.vcf.gz ${VCF%.vcf.gz}.norm.nosnp.vcf.gz /data/variants_in_insertions/SGDP/variantsInInsertions.${SAMP}.norm $REGIONS
    Rscript filterInternalVariants.R /data/variants_in_insertions/SGDP/variantsInInsertions.${SAMP}.norm.tsv /data/variants_in_insertions/SGDP/variantsInInsertions.${SAMP}.norm.filtered.tsv
done



