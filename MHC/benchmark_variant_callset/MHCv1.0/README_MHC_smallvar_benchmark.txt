#Commands for creating benchmark variants and regions for small variants in the MHC

sed 's/^chr//;s/1\/1/1|1/;s/0\/1/0|1/;s/1\/0/1|0/;s/2\/1/2|1/;s/1\/2/1|2/' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_nohead.VCF | cat /Volumes/giab/analysis/mhc/smallvarbenchmark/header_GRCh37.vcf - > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.VCF 

/Applications/bioinfo/htslib-1.3/bgzip -c /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.VCF > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz 

/Applications/bioinfo/htslib-1.3/tabix /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz 
 
cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.0.vcf.gz

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz.tbi /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.0.vcf.gz.tbi


awk '{FS=" ";OFS="\t"} {if($3-$2>200) print $1,$2-int(($3-$2)/10),$3+int(($3-$2)/10)}' /Volumes/giab/analysis/mhc/smallvarbenchmark/low_F_MHC_GRCh37_N2.txt | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -i stdin -d 10000 | awk '$3-$2>1000' > /Volumes/giab/analysis/mhc/smallvarbenchmark/low_F_MHC_GRCh37_N2_gt200_slop10perc_merge10kb_gt1kb.bed


awk 'length($4)>49 || length($5)>49' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head.VCF | awk '{FS=OFS="\t"} {print $1,$2-1,$2+length($4)}' | sed 's/^chr//' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_SVsgt49bp.bed

/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -a /Volumes/giab/analysis/benchmarking-tools/resources/stratification-bed-files/LowComplexity/AllTandemRepeatsandHomopolymers_slop5.bed.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_SVsgt49bp.bed | /Applications/bioinfo/bedtools2.26.0/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_SVsgt49bp.bed |  awk '{FS=OFS="\t"} {print $1,$2-50,$3+50}' |  /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -i stdin -d 1000 > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_SVsgt49bp_repeatexpanded_slop50_merge1000.bed

/Applications/bioinfo/bedtools2.26.0/bin/mergeBed -i /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head.VCF -d 10 -c 3 -o count |  awk '$4>10' | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -i stdin -d 1000 -c 4 -o sum | awk '$4>20' | awk '{FS=OFS="\t"} {print $1,$2-10,$3+10}'  | sed 's/^chr//' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_merge10gt10merge1kgt20_slop10.bed

/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -a /Volumes/giab/analysis/benchmarking-tools/resources/stratification-bed-files/LowComplexity/SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/SimpleRepeat_imperfecthomopolgt10_slop5_MHC_GRCh37.bed 

/Applications/bioinfo/bedtools2.26.0/bin/multiIntersectBed -i /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_merge10gt10merge1kgt20_slop10.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_SVsgt49bp_repeatexpanded_slop50_merge1000.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/low_F_MHC_GRCh37_N2_gt200_slop10perc_merge10kb_gt1kb.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/HLADRB_GRCh37.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/SimpleRepeat_imperfecthomopolgt10_slop5_MHC_GRCh37.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -i stdin -d 100 > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_allexclude.bed

/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -a /Volumes/giab/analysis/benchmarking-tools/resources/stratification-bed-files/LowComplexity/AllTandemRepeatsandHomopolymers_slop5.bed.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_allexclude.bed | /Applications/bioinfo/bedtools2.26.0/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_allexclude.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -i stdin -d 100 > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_allexclude_repeatexpanded.bed

/Applications/bioinfo/bedtools2.26.0/bin/subtractBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_allexclude_repeatexpanded.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.0.bed

/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -v -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_allexclude_repeatexpanded.bed | wc -l
21499
/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed | wc -l
21508

 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed
4581456
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_allexclude_repeatexpanded.bed
389101
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_merge10gt10merge1kgt20_slop10.bed
63500
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_SVsgt49bp_repeatexpanded_slop50_merge1000.bed 
157610
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/low_F_MHC_GRCh37_N2_gt200_slop10perc_merge10kb_gt1kb.bed 
99504
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/SimpleRepeat_imperfecthomopolgt10_slop5_MHC_GRCh37.bed 
87318

/Volumes/DroboZook/bioinfo/rtg-tools-3.10.1/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz -c /Volumes/giab/analysis/v4-development/dv-pacbio-giab-vcfs/HG002.sequelii.case-study.grch37.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002.sequelii.case-study.grch37vsmhcref -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --ref-overlap 
There were 31 problematic baseline variants skipped during loading (see vcfeval.log for details).
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
    3.000              19320          20074        366       2187     0.9821       0.8983     0.9383
     None              19320          20074        366       2187     0.9821       0.8983     0.9383

/Volumes/DroboZook/bioinfo/rtg-tools-3.10.1/rtg vcfeval -b /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002.sequelii.case-study.grch37vsmhcref/fn.vcf.gz -c /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002.sequelii.case-study.grch37vsmhcref/fp.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002.sequelii.case-study.grch37vsmhcref/fnfp -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf  --squash-ploidy --ref-overlap 
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
    3.000                264            290         76       1923     0.7923       0.1207     0.2095
     None                264            290         76       1923     0.7923       0.1207     0.2095



awk '{FS=OFS="\t"} {$2+=32223;$3+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4_GRCh38.bed

#Convert to GRCh38
awk '{FS=OFS="\t"} {$2+=32223;$3+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed | sed 's/^/chr/' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_GRCh38.bed

zgrep -v ^# /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz | awk '{FS=OFS="\t"} {$2+=32223; print}' | sed 's/^/chr/' | cat /Volumes/giab/analysis/mhc/smallvarbenchmark/header_GRCh38.vcf - | /Applications/bioinfo/htslib-1.3/bgzip -c > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz

/Applications/bioinfo/htslib-1.3/tabix /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz 

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_GRCh38.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh38_v1.0.bed

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh38_v1.0.vcf.gz

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz.tbi /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh38_v1.0.vcf.gz.tbi

#combine with v4 benchmark
/Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/subtractBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromassembly.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4.bed

awk '{FS=OFS="\t"} {$2+=32223;$3+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4.bed | sed 's/^/chr/' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4_GRCh38.bed

/Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/subtractBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_head_nochr.vcf.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromassembly.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_nochr_forv4.vcf

awk '{FS=OFS="\t"} {$2+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_nochr_forv4.vcf | sed 's/^/chr/' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_forv4.vcf

#replace GRCh37 MHC with assembly benchmark
/Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_removepartialrepeats.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4.bed | /Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4.bed | /Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/mergeBed -i stdin > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_removepartialrepeats_replaceMHC.bed

/Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf.vcf -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4.bed | cat /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh37_190609_nochr_forv4.vcf - | sort -k2,2n > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC_nohead.vcf

grep ^# /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf.vcf | cat - /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC_nohead.vcf | awk 'length($4)<200 && length($5)<200' | /Volumes/DroboZook/bioinfo/htslib-1.3.2/bgzip -c > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC.vcf.gz

/Volumes/DroboZook/bioinfo/htslib-1.3.2/tabix /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC.vcf.gz 



#replace GRCh38 MHC with assembly benchmark
/Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_removepartialrepeats.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4_GRCh38.bed | /Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4_GRCh38.bed | /Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/mergeBed -i stdin > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_removepartialrepeats_replaceMHC.bed

/Volumes/DroboZook/bioinfo/bedtools2.25.0/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf.vcf -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4_GRCh38.bed | cat /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_forv4.vcf - | sort -k2,2n > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC_nohead.vcf

grep ^# /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf.vcf | cat - /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC_nohead.vcf | awk 'length($4)<200 && length($5)<200' | /Volumes/DroboZook/bioinfo/htslib-1.3.2/bgzip -c > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC.vcf.gz

/Volumes/DroboZook/bioinfo/htslib-1.3.2/tabix /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC.vcf.gz 

/Volumes/DroboZook/bioinfo/rtg-tools-3.10.1/rtg vcfeval -b /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC.vcf.gz -c /Volumes/giab/analysis/v4-development/dv-pacbio-giab-vcfs/HG002.pacbio.sequelii.dv0.8.grch38.vcf.gz -o /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002.pacbio.sequelii.dv0.8.grch38vsv4replacemhc -t /Volumes/giab/data/reference_genomes/rtg_sdf/GRCh38.sdf -e /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_removepartialrepeats_replaceMHC.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh38.bed 
