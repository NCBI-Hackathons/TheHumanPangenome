zcat /home/tobi/reads/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.newsplit.28498559.H1.fastq.gz | awk 'NR % 4 == 1 {print ">"substr($1,2,100)}; NR % 4 == 2 {print}' > tmp.fa
zcat /home/tobi/reads/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.MHConly.newsplit.28498559.untagged.fastq.gz | awk 'NR % 4 == 1 {print ">"substr($1,2,100)}; NR % 4 == 2 {print}' >> tmp.fa
cat ../other_unphased.fa >> tmp.fa
python ../convert_reads.py > reads.fa
echo reads.fa > input.fofn
