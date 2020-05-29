#Commands for creating benchmark variants and regions for small variants in the MHC

#BEDTOOLS=$BEDTOOLS
#HTSLIB=/Volumes/DroboZook/bioinfo/htslib-1.3.2
#RTGTOOLS=/Volumes/DroboZook/bioinfo/rtg-tools-3.10.1
BEDTOOLS=/Applications/bioinfo/bedtools2.26.0
HTSLIB=/Applications/bioinfo/htslib-1.3
RTGTOOLS=/Applications/bioinfo/rtg-tools-3.10.1

#download Chai's dipcall vcf from https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/MHC/benchmark_variant_callset/assembly/dipcall_hg19/repeat_resolved_3.dip.vcf.gz
#Commit https://github.com/NCBI-Hackathons/TheHumanPangenome/commit/a9cd26270a08f8fb7cfced87349b18c29dd6c12a

$HTSLIB/tabix -f /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz

awk '{FS=" ";OFS="\t"} {if($3-$2>200) print $1,$2-int(($3-$2)/10)-2000,$3+int(($3-$2)/10)+2000}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/low_F_MHC_GRCh37_N2_repeat_resolved_3.txt | $BEDTOOLS/bin/mergeBed -i stdin -d 10000 | awk '$3-$2>1000' > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/low_F_MHC_GRCh37_N2_repeat_resolved_3_gt200_slop10perc_merge10kb_gt1kb.bed


gunzip -c /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz | awk 'length($4)>49 || length($5)>49' | awk '{FS=OFS="\t"} {print $1,$2-1,$2+length($4)}' | sed 's/^chr//' > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_SVsgt49bp.bed

$BEDTOOLS/bin/intersectBed -wa -a /Volumes/giab/analysis/benchmarking-tools/resources/stratification-bed-files/LowComplexity/AllTandemRepeatsandHomopolymers_slop5.bed.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_SVsgt49bp.bed | $BEDTOOLS/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_SVsgt49bp.bed |  awk '{FS=OFS="\t"} {print $1,$2-50,$3+50}' |  $BEDTOOLS/bin/mergeBed -i stdin -d 1000 > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_SVsgt49bp_repeatexpanded_slop50_merge1000.bed

$BEDTOOLS/bin/mergeBed -i /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -d 10 -c 3 -o count |  awk '$4>10' | $BEDTOOLS/bin/mergeBed -i stdin -d 1000 -c 4 -o sum | awk '$4>20' | awk '{FS=OFS="\t"} {print $1,$2-10,$3+10}'  | sed 's/^chr//' > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_merge10gt10merge1kgt20_slop10.bed

$BEDTOOLS/bin/intersectBed -wa -a /Volumes/giab/analysis/benchmarking-tools/resources/stratification-bed-files/LowComplexity/SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/SimpleRepeat_imperfecthomopolgt10_slop5_MHC_GRCh37.bed 

$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.complement.bed

$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/mhcvsself -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf  --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --ref-overlap
Evaluation too complex (50002 unresolved paths, 127638 iterations) at reference region 6:29903733-29903877. Variants in this region will not be included in results.
Evaluation too complex (50001 unresolved paths, 82778 iterations) at reference region 6:31310688-31310951. Variants in this region will not be included in results.
Evaluation too complex (50001 unresolved paths, 49887 iterations) at reference region 6:31395645-31395887. Variants in this region will not be included in results.
Evaluation too complex (50002 unresolved paths, 127121 iterations) at reference region 6:31396494-31396633. Variants in this region will not be included in results.
Evaluation too complex (50001 unresolved paths, 79147 iterations) at reference region 6:31396801-31396884. Variants in this region will not be included in results.
Evaluation too complex (50002 unresolved paths, 49794 iterations) at reference region 6:32593702-32594013. Variants in this region will not be included in results.
Evaluation too complex (50001 unresolved paths, 49882 iterations) at reference region 6:32632439-32632599. Variants in this region will not be included in results.
Evaluation too complex (50001 unresolved paths, 101355 iterations) at reference region 6:32670171-32670179. Variants in this region will not be included in results.
There were 19 problematic baseline variants skipped during loading (see vcfeval.log for details).
There were 19 problematic called variants skipped during loading (see vcfeval.log for details).
There were 25851 variants not thresholded in ROC data files due to missing or invalid GQ (FORMAT) values.
Could not maximize F-measure from ROC data, only un-thresholded statistics will be computed. Consider selecting a different scoring attribute with --vcf-score-field
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
     None              25851          25851          0          0     1.0000       1.0000     1.0000

awk '{FS=":\|-";OFS="\t"} { print $5,$6,$7}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/mhcvsself/vcfeval.log | grep 'reference region' | sed 's/.*reference.region.\(.*\)..Varian.*/\1/' > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/mhcvsselftoocomplex.bed

$BEDTOOLS/bin/multiIntersectBed -i /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_merge10gt10merge1kgt20_slop10.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_SVsgt49bp_repeatexpanded_slop50_merge1000.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/low_F_MHC_GRCh37_N2_repeat_resolved_3_gt200_slop10perc_merge10kb_gt1kb.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/HLADRB_GRCh37.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/SimpleRepeat_imperfecthomopolgt10_slop5_MHC_GRCh37.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.complement.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/mhcvsselftoocomplex.bed | $BEDTOOLS/bin/mergeBed -i stdin -d 100 > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_allexclude.bed

$BEDTOOLS/bin/intersectBed -wa -a /Volumes/giab/analysis/benchmarking-tools/resources/stratification-bed-files/LowComplexity/AllTandemRepeatsandHomopolymers_slop5.bed.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_allexclude.bed | $BEDTOOLS/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_allexclude.bed | $BEDTOOLS/bin/mergeBed -i stdin -d 100 > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_allexclude_repeatexpanded.bed

$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_allexclude_repeatexpanded.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.bed

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.vcf.gz

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz.tbi /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.vcf.gz.tbi

$BEDTOOLS/bin/intersectBed -v -a /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_allexclude_repeatexpanded.bed | wc -l
#22445
$BEDTOOLS/bin/intersectBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed | wc -l
#22368

gunzip -c /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz | awk 'length($4)>49 || length($5)>49' | wc -l
#     126
gunzip -c /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz | awk 'length($4)>49' | wc -l
#      63


 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed
#4647177
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_allexclude_repeatexpanded.bed
#323380
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_merge10gt10merge1kgt20_slop10.bed
#18574
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_SVsgt49bp_repeatexpanded_slop50_merge1000.bed 
#68834
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.complement.bed
#119086
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/low_F_MHC_GRCh37_N2_repeat_resolved_3_gt200_slop10perc_merge10kb_gt1kb.bed 
#37157
 awk '{sum+=$3-$2} END {print sum}' /Volumes/giab/analysis/mhc/smallvarbenchmark/SimpleRepeat_imperfecthomopolgt10_slop5_MHC_GRCh37.bed 
#87318


$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/v4-development/dv-pacbio-giab-vcfs/HG002.sequelii.case-study.grch37.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/HG002.sequelii.case-study.grch37vsmhcref -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --ref-overlap 
Evaluation too complex (50001 unresolved paths, 82790 iterations) at reference region 6:31310689-31310951. Variants in this region will not be included in results.
Evaluation too complex (50002 unresolved paths, 49954 iterations) at reference region 6:31395645-31395910. Variants in this region will not be included in results.
Evaluation too complex (50002 unresolved paths, 49816 iterations) at reference region 6:32593702-32594013. Variants in this region will not be included in results.
Evaluation too complex (50001 unresolved paths, 49888 iterations) at reference region 6:32632440-32632599. Variants in this region will not be included in results.
There were 18 problematic baseline variants skipped during loading (see vcfeval.log for details).
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
    3.000              20152          20168        256       2213     0.9875       0.9011     0.9423
     None              20152          20168        256       2213     0.9875       0.9011     0.9423



$RTGTOOLS/rtg vcfeval -b /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/HG002.sequelii.case-study.grch37vsmhcref/fn.vcf.gz -c /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/HG002.sequelii.case-study.grch37vsmhcref/fp.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/HG002.sequelii.case-study.grch37vsmhcref/fnfp -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf  --squash-ploidy --ref-overlap 
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
    3.000                 27             26         24        188     0.5200       0.1256     0.2023
     None                 27             26         24        188     0.5200       0.1256     0.2023

$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/v4-development/v4.1_draft_evaluations/CCS_Clair/hs37d5_sorted.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/CCSClair.grch37vsmhcref -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --ref-overlap 
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
  103.000              21189          21191        346       1176     0.9839       0.9474     0.9653
     None              21215          21218        432       1150     0.9800       0.9486     0.9641



$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/v4-development/v4.1_before_MHC_replacement/HG002_GRCh37_1_22_v4.1_python_highconf.vcf.gz  -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/v4.1beforeMHCvsdipcallmhcrefbothbeds -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed -e /Volumes/giab/analysis/v4-development/v4.1_before_MHC_replacement/HG002_GRCh37_1_22_v4.1_python_highconf.bed --ref-overlap 
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
   33.000              14477          14470          3         22     0.9998       0.9985     0.9991
     None              14477          14470          3         22     0.9998       0.9985     0.9991

$BEDTOOLS/bin/intersectBed -a /Volumes/giab/analysis/v4-development/v4.1_before_MHC_replacement/HG002_GRCh37_1_22_v4.1_python_highconf.vcf.gz -b /Volumes/giab/analysis/v4-development/v4.1_before_MHC_replacement/HG002_GRCh37_1_22_v4.1_python_highconf.bed -header | $BEDTOOLS/bin/intersectBed -a stdin -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed | wc -l
#14999

$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/v4-development/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/v4.1vsdipcallmhcrefbothbeds -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed -e /Volumes/giab/analysis/v4-development/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed --ref-overlap 
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
     None              15588          15161         13          9     0.9991       0.9994     0.9993

$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/v4-development/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/v4.1vsdipcallmhcref4.1bedonly -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --bed-regions /Volumes/giab/analysis/v4-development/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed --ref-overlap 



$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/v4-development/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/v4.1vsdipcallmhcref -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --ref-overlap 

#compare to v1.0
$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/v1.0vsdipcallmhcrefbothbeds -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.bed --ref-overlap 
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
     None              17039          16520         11          8     0.9993       0.9995     0.9994

$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/v1.0vsdipcallmhcrefv1.0bedonly -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.bed --ref-overlap 

$RTGTOOLS/rtg vcfeval -b  /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip.vcf.gz -c /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh37_v1.1.vcf.gz -o /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/v1.0vsdipcallmhcrefdipcallbedonly -t /Volumes/giab/data/reference_genomes/rtg_sdf/1000g_v37_phase2.sdf -e /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/MHC_highconf.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh37.bed --ref-overlap 


awk '{FS=OFS="\t"} {$2+=32223;$3+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4_GRCh38.bed

#Convert to GRCh38
awk '{FS=OFS="\t"} {$2+=32223;$3+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed | sed 's/^/chr/' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_GRCh38.bed

zgrep -v ^# /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_head_nochr.vcf.gz | awk '{FS=OFS="\t"} {$2+=32223; print}' | sed 's/^/chr/' | cat /Volumes/giab/analysis/mhc/smallvarbenchmark/header_GRCh38.vcf - | /Applications/bioinfo/htslib-1.3/bgzip -c > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz

/Applications/bioinfo/htslib-1.3/tabix /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz 

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_GRCh38.bed /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh38_v1.1.bed

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh38_v1.1.vcf.gz

cp /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_GRCh38_head.vcf.gz.tbi /Volumes/giab/analysis/mhc/smallvarbenchmark/HG002_MHCdiploidassembly_GRCh38_v1.1.vcf.gz.tbi

#combine with v4 benchmark
$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromassembly.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4.bed

awk '{FS=OFS="\t"} {$2+=32223;$3+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4.bed | sed 's/^/chr/' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4_GRCh38.bed

$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_head_nochr.vcf.gz -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromassembly.bed > /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_nochr_forv4.vcf

awk '{FS=OFS="\t"} {$2+=32223; print}' /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_nochr_forv4.vcf | sed 's/^/chr/' > /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_forv4.vcf

#replace GRCh37 MHC with assembly benchmark
$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_removepartialrepeats.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4.bed | $BEDTOOLS/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4.bed | $BEDTOOLS/bin/mergeBed -i stdin > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_removepartialrepeats_replaceMHC.bed

$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf.vcf -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4.bed | cat /Volumes/giab/analysis/mhc/smallvarbenchmark/dipcall/repeat_resolved_3.dip_nochr_forv4.vcf - | sort -k2,2n > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC_nohead.vcf

grep ^# /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf.vcf | cat - /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC_nohead.vcf | awk 'length($4)<200 && length($5)<200' | $HTSLIB/bgzip -c > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC.vcf.gz

$HTSLIB/tabix /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh37_6_v4beta_highconf_replaceMHC.vcf.gz 



#replace GRCh38 MHC with assembly benchmark
$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_removepartialrepeats.bed -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4_GRCh38.bed | $BEDTOOLS/bin/multiIntersectBed -i stdin /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_highconf_forv4_GRCh38.bed | $BEDTOOLS/bin/mergeBed -i stdin > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_removepartialrepeats_replaceMHC.bed

$BEDTOOLS/bin/subtractBed -a /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf.vcf -b /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_excludeforintegrationfromv4_GRCh38.bed | cat /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_0525.minimap2GRCh38_190609_forv4.vcf - | sort -k2,2n > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC_nohead.vcf

grep ^# /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf.vcf | cat - /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC_nohead.vcf | awk 'length($4)<200 && length($5)<200' | $HTSLIB/bgzip -c > /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC.vcf.gz

$HTSLIB/tabix /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC.vcf.gz 

$RTGTOOLS/rtg vcfeval -b /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_replaceMHC.vcf.gz -c /Volumes/giab/analysis/v4-development/dv-pacbio-giab-vcfs/HG002.pacbio.sequelii.dv0.8.grch38.vcf.gz -o /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002.pacbio.sequelii.dv0.8.grch38vsv4replacemhc -t /Volumes/giab/data/reference_genomes/rtg_sdf/GRCh38.sdf -e /Volumes/giab/analysis/v4-development/MHC_benchmark_insertion/HG002_GRCh38_6_v4beta_highconf_removepartialrepeats_replaceMHC.bed --bed-regions /Volumes/giab/analysis/mhc/smallvarbenchmark/MHC_GRCh38.bed 
