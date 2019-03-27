The data was downloaded from the following links:

PacBio CCS 15Kb Sequel I - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/alignment/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
PacBio CCS 10Kb Sequel I - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
PacBio CCS 11Kb Sequel II - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
ONT "ultralong" - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/combined_2018-08-10/HG002_ONTrel2_16x_RG_HP10xtrioRTG.cram
ONT Promethion - https://s3-us-west-2.amazonaws.com/human-pangenomics/agbt2019/fastq/GM24385.fastq
10X Genomics - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/10XGenomics/NA24385_phased_possorted_bam.bam


ONT Promethion data was aligned and indexed using the following commands:

$: minimap2 -t 16 -aL -z 600,200 -x map-ont -R '@RG\tID:ucsc1\tPU:1\tPL:nanopore\t:PM:promethion\tLB:1\tSM:HG002' hs37d5.fa GM24385.fastq | samtools sort -m 1G -@16 -O bam --reference hs37d5.fa > GM24385.bam
$: samtools index GM24385.bam

For PacBio CCS 10Kb, 11Kb, 15Kb, ONT "ultralong", and 10X Genomics the reads covering the MHC were selected using the following commands (with bam file name updated for each data) for haplotype 1 and haplotype 2:

$: samtools view -bh HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.bam 6:28477797-33448354 | bamtools filter -in stdin -tag HP:1 | samtools bam2fq - | bgzip -c > HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.HP1.fastq.gz
$: samtools view -bh HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.bam 6:28477797-33448354 | bamtools filter -in stdin -tag HP:2 | samtools bam2fq - | bgzip -c > HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.HP2.fastq.gz

For the ONT Promethion data reads covering the MHC was generated using following command:

$: samtools view -bh GM24385.bam 6:28477797-33448354 | samtools bam2fq - | bgzip -c > HG002_Promethion.MHConly.fastq.gz

