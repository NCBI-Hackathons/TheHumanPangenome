Notes on haplotype partitioning reads in the MHC
================================================
Note: this is not fully complete, in particular commands to sort/index CRAMs/VCFs are omitted

Coordinate to be used: 6:28477797-33448354 for GRCh37
 
Pulling software
```
conda create -n pangenome-hackathon whatshap minimap2 samtools bcftools
```
 
Getting the input files
```
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh37/NA24385_300G/NA24385.GRCh37.phased_variants.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
```

Merging FASTQs again
```
zcat data/HG002_ONTrel2_16x_RG_HP10xtrioRTG.MHConly.*.fastq.gz | pigz > reads/HG002_ONTrel2_16x_RG_HP10xtrioRTG.MHConly.fastq.gz
zcat data/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.*.fastq.gz |pigz > reads/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.fastq.gz
```

Mapping CCS reads

```
minimap2 -t 4 -R '@RG\tID:1\tSM:HG002' -ax asm20 ref/hs37d5.fa.gz reads/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.fastq.gz | samtools view -T ref/hs37d5.fa.gz -C - > cram/HG002.PacBio.15kbCCS.cram
```

Mapping ONT reads

```
minimap2 -t 4 -R '@RG\tID:1\tSM:HG002' -ax map-ont ref/hs37d5.fa.gz reads/HG002_ONTrel2_16x_RG_HP10xtrioRTG.MHConly.fastq.gz | samtools view -T ref/hs37d5.fa.gz -C - > cram/HG002_ONTrel2_16x_RG_HP10xtrioRTG.cram
```

Select het SNVs from DeepVariant call set. The goal is to use this as a starting point to construct a super confident set of heterozygous SNVs that we can use to partition reads
```
bcftools view -g het -v snps vcf/deepvariant/pacbio-15kb-hapsort-wgs.vcf.gz 6:28477797-33448354 > vcf/deepvariant/mhc.het.snps.vcf
```

Re-genotype these SNV from ONT reads
```
whatshap genotype --ignore-read-groups --chromosome 6 --reference ref/hs37d5.fa -o whatshap/genotype/HG002.ONTrel2_16x.vcf vcf/deepvariant/mhc.het.snps.vcf cram/HG002_ONTrel2_16x_RG_HP10xtrioRTG.sorted.cram > whatshap/genotype/HG002.ONTrel2_16x.log 2>&1
```

Re-genotype these SNV from CCS reads
```
whatshap genotype --ignore-read-groups --chromosome 6 --reference ref/hs37d5.fa -o whatshap/genotype/HG002.PacBio.15kbCCS.vcf vcf/deepvariant/mhc.het.snps.vcf cram/HG002.PacBio.15kbCCS.sorted.cram > whatshap/genotype/HG002.PacBio.15kbCCS.log 2>&1
```

Find the intersection, i.e. only retain those that are het accoring to DV/CCS, WhatHap/ONT, and WhatsHap/CCS
```
bcftools view -g het whatshap/genotype/HG002.ONTrel2_16x.vcf > whatshap/genotype/HG002.ONTrel2_16x.hetonly.vcf
bcftools view -g het whatshap/genotype/HG002.PacBio.15kbCCS.vcf > whatshap/genotype/HG002.PacBio.15kbCCS.hetonly.vcf
bcftools isec -p whatshap/genotype-intersect whatshap/genotype/HG002.ONTrel2_16x.hetonly.vcf.gz whatshap/genotype/HG002.PacBio.15kbCCS.hetonly.vcf.gz
cp ~/whatshap/genotype-intersect/0003.vcf ~/whatshap/confident-hets.vcf
```

Extract 10x phased blocks produced by LongRanger 
```
cd 10x
bcftools view NA24385.GRCh37.phased_variants.vcf.gz chr6:28477797-33448354 | awk 'BEGIN {OFS="\t"} $1 == "chr6" {$1="6"} $1=="#CHROM" {$10="HG002"} {print}'|bgzip > NA24385.GRCh37.phased_variants.mhc.vcf.gz
```

Now run phasing from multiple combinations of data sources
```
whatshap phase --ignore-read-groups --chromosome 6 --reference ref/hs37d5.fa -o whatshap/phase/HG002.PacBio.15kbCCS.vcf ~/whatshap/confident-hets.vcf cram/HG002.PacBio.15kbCCS.sorted.cram > whatshap/phase/HG002.PacBio.15kbCCS.log 2>&1
whatshap phase --ignore-read-groups --chromosome 6 --reference ref/hs37d5.fa -o whatshap/phase/HG002.ONTrel2_16x.vcf ~/whatshap/confident-hets.vcf cram/HG002_ONTrel2_16x_RG_HP10xtrioRTG.sorted.cram > whatshap/phase/HG002.ONTrel2_16x.log 2>&1
whatshap phase --ignore-read-groups --chromosome 6 --reference ref/hs37d5.fa -o whatshap/phase/HG002.PacBio+10x.vcf ~/whatshap/confident-hets.vcf cram/HG002.PacBio.15kbCCS.sorted.cram 10x/NA24385.GRCh37.phased_variants.mhc.vcf.gz > whatshap/phase/HG002.PacBio+10x.log 2>&1
```

Tag reads and split FASTQs
```
whatshap haplotag --ignore-read-groups --output-haplotag-list cram/HG002.PacBio.15kbCCS.sorted.tags.tsv --reference ref/hs37d5.fa -o cram/HG002.PacBio.15kbCCS.sorted.tagged.cram whatshap/phase/HG002.PacBio+10x.vcf.gz cram/HG002.PacBio.15kbCCS.sorted.cram
whatshap haplotag --ignore-read-groups --output-haplotag-list cram/HG002_ONTrel2_16x_RG_HP10xtrioRTG.sorted.tags.tsv --reference ref/hs37d5.fa -o cram/HG002_ONTrel2_16x_RG_HP10xtrioRTG.sorted.tagged.cram whatshap/phase/HG002.PacBio+10x.vcf.gz cram/HG002_ONTrel2_16x_RG_HP10xtrioRTG.sorted.cram
./split-fastq.sh reads/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.fastq.gz cram/HG002.PacBio.15kbCCS.sorted.tags.tsv
```

Assigning Jason's haplotigs based on parental data
```
cd jason_assemblies
cp /data/assemblies_from_Tobias_partitioned/*.fa .
for f in p_ctg_cns-* ; do sed -E "s/^>(.*)/>${f}.\1/g" $f ; done > partitioned-assemblies.fa 
minimap2 -t 4 -R '@RG\tID:1\tSM:HG002' -ax asm20 ref/hs37d5.fa.gz jason_assemblies/partitioned-assemblies.fa | samtools view -T ref/hs37d5.fa.gz -C - > cram/HG002.partitioned-assemblies.cram
whatshap unphase AJ.6.vcf.gz |bgzip > AJ.6.unphased.vcf.gz
bcftools isec -p rtg-isec whatshap/confident-hets.vcf.gz rtg/AJ.6.unphased.vcf.gz
whatshap phase --ped ped/AJ.ped rtg-isec/0003.vcf |bgzip > rtg-isec/trio-phased.vcf.gz
whatshap haplotag --reference ref/hs37d5.fa --output-haplotag-list cram/HG002.partitioned-assemblies.sorted.taglist.tsv -o cram/HG002.partitioned-assemblies.sorted.tagged.cram  rtg-isec/trio-phased.vcf.gz cram/HG002.partitioned-assemblies.sorted.cram
```

Now switch to brand new PromethIon data
=======================================
```
minimap2 -t 4 -R '@RG\tID:1\tSM:HG002' -ax map-ont ref/hs37d5.fa.gz reads/HG002_Promethion.MHConly.fastq.gz | samtools view -T ref/hs37d5.fa.gz -C - > cram/HG002_Promethion.cram
 
whatshap phase --chromosome 6 --reference ref/hs37d5.fa -o whatshap/phase/HG002.Promethion+10x.vcf ~/whatshap/confident-hets.vcf cram/HG002_Promethion.sorted.cram 10x/NA24385.GRCh37.phased_variants.mhc.vcf.gz > whatshap/phase/HG002.Promethion+10x.log 2>&1

whatshap haplotag --ignore-read-groups --output-haplotag-list cram/HG002.PacBio.15kbCCS.sorted.tags-newsplit.tsv --reference ref/hs37d5.fa -o cram/HG002.PacBio.15kbCCS.sorted.tagged-newsplit.cram whatshap/phase/HG002.Promethion+10x.vcf.gz cram/HG002.PacBio.15kbCCS.sorted.cram

whatshap haplotag --ignore-read-groups --output-haplotag-list cram/HG002.Promethion.sorted.tags-newsplit.tsv --reference ref/hs37d5.fa -o cram/HG002.Promethion.sorted.tagged-newsplit.cram whatshap/phase/HG002.Promethion+10x.vcf.gz cram/HG002_Promethion.sorted.cram

./split-fastq.sh reads/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.newsplit.fastq.gz cram/HG002.PacBio.15kbCCS.sorted.tags-newsplit.tsv
```

`split-fastq.sh` script:

```
#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.fastq.gz> <taglist.tsv>"
    exit 1
fi

fastq=$1
listfile=$2


for block in $(tail -n +2 $listfile |cut -f 3|grep -v 'none' | sort |uniq) ; do
	echo "Block $block"
	whatshap split --output-h1 ${fastq%.fastq.gz}.$block.H1.fastq.gz --output-h2 ${fastq%.fastq.gz}.$block.H2.fastq.gz $fastq <(awk "\$3==$block" $listfile)
done
```
