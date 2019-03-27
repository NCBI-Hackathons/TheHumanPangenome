#!/bin/bash

zcat /Projects/DS/grch38+jg1/JG1.0.0beta.fa.gz | bgzip -c > JG1.fa.gz
samtools faidx JG1.fa.gz
samtools faidx JG1.fa.gz chr21 > JG1.chr21.fa.gz

samtools faidx GRCh38.fa.gz CM000683.2 | bgzip -c > GRCh38.chr21.fa.gz
zcat *.chr21.fa.gz | bgzip -c > chr21.fa.gz

minimap2 -x asm5 -X -c -t 8 chr21.fa.gz chr21.fa.gz > chr21.paf 2> err
seqwish-6e4fe705 -p chr21.paf -s chr21.fa.gz -b chr21.graph -g chr21.gfa -t 8

cat <(zcat GRCh38.chr21.fa.gz) /home/ys/data/PanTro_Clint_PTRv2/PanTro_Clint_PTRv2_Chr21.fasta /home/ys/data/PonAbe3_Susie_PABv2/PonAbe3_Susie_PABv2_Chr21.fasta /home/ys/data/CHM1_BAC_clones/AC244111.3.fasta /home/ys/data/CHM1_BAC_clones/AC244144.2.fasta /home/ys/data/CHM1_BAC_clones/AC244518.2.fasta /home/ys/data/CHM1_BAC_clones/AC245051.3.fasta /home/ys/data/CHM1_BAC_clones/AC245314.2.fasta /home/ys/data/CHM1_BAC_clones/AC246819.2.fasta /home/ys/data/CHM1_BAC_clones/AC255431.1.fasta /home/ys/data/CHM1_BAC_clones/AC256301.1.fasta /home/ys/data/CHM1_BAC_clones/AC277730.1.fasta /home/ys/data/CHM1_BAC_clones/AC277802.1.fasta /home/ys/data/CHM1_BAC_clones/AC277887.1.fasta | bgzip -c > panApe.chr21.fa.gz

minimap2 -x asm5 -X -c -t 8 panApe.chr21.fa.gz panApe.chr21.fa.gz > panApe.chr21.paf 2> panApe.err
./seqwish-6e4fe705 -s panApe.chr21.fa.gz -p panApe.chr21.paf -b panApe.chr21.graph -g panApe.chr21.gfa -t 8 > panApe.out 2> panApe.err

minimap2 -x asm10 -X -c -t 8 panApe.chr21.fa.gz panApe.chr21.fa.gz > panApe.chr21.relaxed.paf 2> panApe.relaxed.err
./seqwish-6e4fe705 -s panApe.chr21.fa.gz -p panApe.chr21.relaxed.paf -b panApe.chr21.relaxed.graph -g panApe.chr21.relaxed.gfa -t 8 >> panApe.out 2>> panApe.err

####On vg docker
/vg/bin/vg view -F -v /home/panApe.chr21.gfa > /home/panApe.chr21.vg
/vg/bin/vg index -x /home/panApe.chr21.xg /home/panApe.chr21.vg
