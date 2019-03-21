# Pangenomics data definition

This is a description of the human read sets that will be used to build the "graphs." The RNAseq data will be a subset of the the same data loaded to GCP for the RNAseq hackathon. 

## GRCh38 Data


| Assembly Accession/Version | Path | Comment |
|----------|------------|----------|
| GCA_000001405.28 (GRCh38.p13) | [FTP link](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13) | The patches related data is in the follwing sub-directory: [FTP link](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_structure/) | It's not ready-made for analysis pipelines, but I will be encouraging teams to use the patch scaffolds in the graph building effort.|
|GCA_000001405.14 (GRCh38) | [FTP link](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines) | This data was probably in the RNAseq hackathon: can we confirm what we've got already? We want to be sure we've got the versions of the analysis sets that include the alternate loci, as that is considered content relevant for the graph building.|

## GOLD data in SRA
(Gold: https://www.genome.wustl.edu/items/reference-genome-improvement/)

|BioSample|Gold_Sample|Origin|TrioMember|SRA_Accession|Assembly_Accession|Comment|
|----|----|----|----|----|----|----|
|SAMN04229549|HG00731|Puerto Rican|father|SRX2802379|NA|3kb Nextera illumina library 2 of HG00731|
|SAMN04229549|HG00731|Puerto Rican|father|SRX2802378|NA|3kb Nextera illumina library 1 of HG00731|
|SAMN04229549|HG00731|Puerto Rican|father|SRX2802347|NA|PCR-free paired end library with overlapping reads of HG00731|
|SAMN04229549|HG00731|Puerto Rican|father|SRX2798634|NA|TruSeq PCR-free Illumina library of HG00731|
|SAMN04229550|HG00732|Puerto Rican|mother|SRX2802376|NA|3kb Nextera illumina library 2 of HG00732|
|SAMN04229550|HG00732|Puerto Rican|mother|SRX2802372|NA|3kb Nextera illumina library 1 of HG00732|
|SAMN04229550|HG00732|Puerto Rican|mother|SRX2802057|NA|PCR-free paired end library with overlapping reads of HG00732|
|SAMN04229550|HG00732|Puerto Rican|mother|SRX2798624|NA|TruSeq PCR-free Illumina library of HG00732|
|SAMN04229548|HG00733|Puerto Rican|child|NA|GCA_002208065.1|de novo assembly (draft)|
|SAMN04229548|HG00733|Puerto Rican|child|SRX2804194|NA|3kb Nextera illumina library 2 of HG00733|
|SAMN04229548|HG00733|Puerto Rican|child|SRX2804182|NA|3kb Nextera illumina library 1 of HG00733|
|SAMN04229548|HG00733|Puerto Rican|child|SRX2802055|NA|PCR-free paired end library with overlapping reads of HG00733|
|SAMN04229548|HG00733|Puerto Rican|child|SRX2798623|NA|TruSeq PCR-free Illumina library of HG00733|
|SAMN04229548|HG00733|Puerto Rican|child|SRX1458652|NA|Other Sequencing of human Puerto Rican (PacBio RS)|
|SAMN04229548|HG00733|Puerto Rican|child|SRX1458080|NA|Other Sequencing of human Puerto Rican (PacBio RS)|
|SAMN04229548|HG00733|Puerto Rican|child|SRX1457808|NA|Other Sequencing of human Puerto Rican (PacBio RS)|
|SAMN03838746|NA19240|Yoruban|child|NA|GCA_001524155.4|de novo assembly (draft)|
|SAMN03838746|NA19240|Yoruban|child|SRX1098167|NA|TruSeq PCR-free Illumina library of NA19240|
|SAMN03838746|NA19240|Yoruban|child|SRX1098166|NA|PCR-free paired end library with overlapping reads of NA19240|
|SAMN03838746|NA19240|Yoruban|child|SRX1098165|NA|3kb Nextera illumina library 2|
|SAMN03838746|NA19240|Yoruban|child|SRX1098164|NA|3kb Nextera illumina library 1 of NA19240|
|SAMN03838746|NA19240|Yoruban|child|SRX1098163|NA|SWIFT illumina library with PCR and Blue Pippin of NA19240|
|SAMN03838746|NA19240|Yoruban|child|SRX1098162|NA|PCR-free SWIFT illumina library of NA19240|
|SAMN03838746|NA19240|Yoruban|child|SRX1096798|NA|Other Sequencing of human Yoruban (PacBio RS)|
|SAMN03838746|NA19240|Yoruban|child|SRX1094388|NA|Other Sequencing of human Yoruban (PacBio RS)|
|SAMN03838746|NA19240|Yoruban|child|SRX1094374|NA|Other Sequencing of human Yoruban (PacBio RS)|
|SAMN03838746|NA19240|Yoruban|child|SRX1094289|NA|Other Sequencing of human Yoruban (PacBio RS)|
|SAMN03838746|NA19240|Yoruban|child|SRX1093654|NA|Other Sequencing of human Yoruban (PacBio RS)|
|SAMN03838746|NA19240|Yoruban|child|SRX1093555|NA|Other Sequencing of human Yoruban (PacBio RS)|
|SAMN03838746|NA19240|Yoruban|child|SRX1093000|NA|Other Sequencing of human Yoruban (PacBio RS)|
|SAMN04229552|HG00514|Han Chinese|child|NA|GCA_002180035.3|de novo assembly (draft)|
|SAMN04229552|HG00514|Han Chinese|child|SRX2802354|NA|3kb Nextera illumina library 2 of HG00514|
|SAMN04229552|HG00514|Han Chinese|child|SRX2802353|NA|3kb Nextera illumina library 1 of HG00514|
|SAMN04229552|HG00514|Han Chinese|child|SRX2802056|NA|PCR-free paired end library with overlapping reads of HG00514|
|SAMN04229552|HG00514|Han Chinese|child|SRX2798895|NA|TruSeq PCR-free Illumina library of HG00514|
|SAMN04229552|HG00514|Han Chinese|child|SRX1619615|NA|Other Sequencing of human Han Chinese(Southern) (PacBio RS)|
|SAMN04229552|HG00514|Han Chinese|child|SRX1619015|NA|Other Sequencing of human Han Chinese(Southern) (PacBio RS)|
|SAMN04229552|HG00514|Han Chinese|child|SRX1617978|NA|Other Sequencing of human Han Chinese(Southern) (PacBio RS)|
|SAMN05181962|NA12878|Utah/CEPH|mother|NA|GCA_002077035.3|de novo assembly (draft)|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX2802351|NA|3kb Nextera illumina library 2 of NA12878|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX2802349|NA|3kb Nextera illumina library 1 ofNA12878|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX2802054|NA|PCR-free paired end library with overlapping reads of NA12878|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX2798622|NA|TruSeq PCR-free Illumina library of NA12878|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX1837675|NA|Other Sequencing of human Ceph/Utah/Mormon (PacBio RS)|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX1837653|NA|Other Sequencing of human Ceph/Utah/Mormon (PacBio RS)|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX1837275|NA|Other Sequencing of human Ceph/Utah/Mormon (PacBio RS)|
|SAMN05181962|NA12878|Utah/CEPH|mother|SRX1837266|NA|Other Sequencing of human Ceph/Utah/Mormon (PacBio RS)|
|SAMN05603729|HG01352|Columbian|child|NA|GCA_002209525.2|de novo assembly (draft)|
|SAMN05603729|HG01352|Columbian|child|SRX2830480|NA|PCR-free illumina library for HG01352|
|SAMN05603729|HG01352|Columbian|child|SRX2095531|NA|Other Sequencing of human COLOMBIAN IN MEDELLIN, COLOMBIA (PacBio RS)|
|SAMN05603729|HG01352|Columbian|child|SRX2095500|NA|Other Sequencing of human COLOMBIAN IN MEDELLIN, COLOMBIA (PacBio RS)|
|SAMN05603729|HG01352|Columbian|child|SRX2095346|NA|Other Sequencing of human COLOMBIAN IN MEDELLIN, COLOMBIA (PacBio RS)|
|SAMN06885952|NA19434|Luhya|child|NA|GCA_002872155.1|de novo assembly (draft)|
|SAMN06885952|NA19434|Luhya|child|SRX2830521|NA|PCR-free illumina library for NA19434|
|SAMN06885952|NA19434|Luhya|child|SRX4118367|NA|Whole Genome Sequencing of human Webuye, KENYA (PacBio RSII)|
|SAMN06885952|NA19434|Luhya|child|SRX3102059|NA|Whole Genome Sequencing of human Webuye, KENYA (PacBio RSII)|
|SAMN05603847|HG02059|Kinh-Vietnamese|child|NA|GCA_003070785.1|de novo assembly (draft)|
|SAMN05603847|HG02059|Kinh-Vietnamese|child|SRX2830492|NA|PCR-free illumina library for HG02059|
|SAMN05603847|HG02059|Kinh-Vietnamese|child|SRX2537696|NA|Other Sequencing of human VIETNAMESE FROM KINH IN HO CHI MINH CITY, VIETNAM (PacBio RSII)|
|SAMN05603847|HG02059|Kinh-Vietnamese|child|SRX2537695|NA|Other Sequencing of human VIETNAMESE FROM KINH IN HO CHI MINH CITY, VIETNAM (PacBio RSII)|
|SAMN05603847|HG02059|Kinh-Vietnamese|child|SRX2537694|NA|Other Sequencing of human VIETNAMESE FROM KINH IN HO CHI MINH CITY, VIETNAM (PacBio RSII)|
|SAMN08723473|HG03486|Mende|child|NA|GCA_003086635.1|de novo assembly (draft)|
|SAMN08723473|HG03486|Mende|child|SRX4867869|NA|Other Sequencing of human (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|NA|GCA_003574075.1|de novo assembly (draft)|
|SAMN05603745|HG02818|Gambian|child|SRX2830485|NA|PCR-free illumina library for HG02818|
|SAMN05603745|HG02818|Gambian|child|SRX4904906|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3203304|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3203303|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3203302|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3115514|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3115513|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3115512|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3115511|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN05603745|HG02818|Gambian|child|SRX3115510|NA|Other Sequencing of human GAMBIAN IN WESTERN DIVISION OF GAMBIA (PacBio RSII)|
|SAMN10026989|HG03807|Bengali|child|NA|GCA_003601015.1|de novo assembly (draft)|
|SAMN10026989|HG03807|Bengali|child|SRX4905637|NA|Whole Genome Sequencing of human (PacBio RSII)|
|SAMN10026989|HG03807|Bengali|child|SRX4808110|NA|Whole Genome Sequencing of human (PacBio RSII)|
|SAMN10026989|HG03807|Bengali|child|SRX4807417|NA|Whole Genome Sequencing of human (PacBio RSII)|
|SAMN10026989|HG03807|Bengali|child|SRX4807281|NA|Whole Genome Sequencing of human (PacBio RSII)|
|SAMN10026989|HG03807|Bengali|child|SRX4807258|NA|Whole Genome Sequencing of human (PacBio RSII)|
|SAMN10026989|HG03807|Bengali|child|SRX4806918|NA|Whole Genome Sequencing of human (PacBio RSII)|
|SAMN10026989|HG03807|Bengali|child|SRX4804704|NA|Whole Genome Sequencing of human (PacBio RSII)|
|SAMN09690649|HG04217|Telegu|child|295 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=481794|
|SAMN09651199|HG02106|Peruvian|child|45 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=480858|
|SAMN09643900|HG00268|Finnish|child|30 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=480712|
|SAMN10390358|NA19836|African American|child|NA|NA|No SRA data|
|SAMN10432089|HG03125|Esan|child|14 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=505678|
|SAMN02744161|CHM1|NA|NA|NA|GCA_001297185.2|de novo assembly (Pilon corrected); there are additional uncorrected assemblies from the SRA experiments noted below|
|SAMN02744161|CHM1|NA|NA|140 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=246220|
|SAMN02744161|CHM1|NA|NA|SRX533609|NA|Pacific Biosciences Model Organism Genome Sequencing-Human 54X (P5/C3)|
|SAMN03255769|CHM13|NA|NA|NA|GCA_000983455.2|de novo assembly (Pilon corrected); there are additional uncorrected assemblies from the SRA experiments noted below|
|SAMN03255769|CHM13|NA|NA|128 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=269593|

## Genome in a bottle data in SRA

|----|----|----|----|----|----|----|
|SAMN03283347|NA24385/HG002|Ashkenazi|child|NA|GCA_001542345.1|de novo assembly|
|SAMN03283347|NA24385/HG002|Ashkenazi|child|188 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?LinkName=biosample_sra&from_uid=3283347|
|SAMN03283345|NA24149/HG003|Ashkenazi|father|NA|GCA_001549605.1|de novo assembly|
|SAMN03283345|NA24149/HG003|Ashkenazi|father|193 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?LinkName=biosample_sra&from_uid=3283345|
|SAMN03283346|NA24143/HG004|Ashkenazi|mother|NA|GCA_001549595.1|de novo assembly|
|SAMN03283346|NA24143/HG004|Ashkenazi|mother|186 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?LinkName=biosample_sra&from_uid=3283346|
|SAMN03283350 |NA24631/HG005|Chinese|child|SRX4739017|NA|Chinese trio PacBio Sequel data|
|SAMN03283348 |NA24694/HG006|Chinese|father|SRX4739121|NA|Chinese trio PacBio Sequel data|
|SAMN03283349 |NA24695/HG007|Chinese|mother|SRX4739122|NA|Chinese trio PacBio Sequel data|

## Other Data in SRA:

|SAMEA104349931|NA12878|Utah/CEPH|mother|NA|GCA_900232925.2|de novo assembly; not haplotype resolved. ONT sequencing.|
|SAMEA104349931|NA12878|Utah/CEPH|mother|53 experiments|NA|https://www.ncbi.nlm.nih.gov/sra?LinkName=biosample_sra&from_uid=7787279|

## Haplotype Resolved BioNano Maps
Links to the public BioNano map files can be found in this file on [the FTP site](http://ftp.ncbi.nlm.nih.gov/pub/supplementary_data/bionanomaps.csv)
(That file contains the SUPPF_accns).

|BioSample|Gold_Sample|Origin|BioProject|Suppfile_Accn|File Types|
|----|----|----|----|----|----|
|SAMN10026989|HG03807|Bengali|PRJNA490190|SUPPF_0000002915; SUPPF_0000002926|.bnx, .cmap|
|SAMN08723473|HG03486|Mende|PRJNA438669|SUPPF_0000002914; SUPPF_0000002924|.bnx, .cmap|
|SAMN06885952|NA19434|Luhya|PRJNA385272|SUPPF_0000002909; SUPPF_0000002919|.bnx, .cmap|
|SAMN05603847|HG02059|Kinh-Vietnamese|PRJNA339726|SUPPF_0000002913; SUPPF_0000002923|.bnx, .cmap|
|SAMN05603745|HG02818|Gambian|PRJNA339722|SUPPF_0000002915; SUPPF_0000002925|.bnx, .cmap|
|SAMN05603729|HG01352|Columbian|PRJNA339719|SUPPF_0000002912; SUPPF_0000002922|.bnx, .cmap|
|SAMN04229552|HG00514|Han Chinese|PRJNA300843|SUPPF_0000002910; SUPPF_0000002920|.bnx, .cmap|
|SAMN04229548|HG00733|Puerto Rican|PRJNA300840|SUPPF_0000002911; SUPPF_0000002921|.bnx, .cmap|
|SAMN03838746|NA19240|Yoruban|PRJNA288807|SUPPF_0000002908; SUPPF_0000002918|.bnx, .cmap|
|SAMN03255769|CHM13|NA|PRJNA269593|SUPPF_0000002917|.bnx, .cmap|
