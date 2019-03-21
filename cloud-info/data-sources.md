# Pangenomics data definition

This is a description of the human read sets that will be used to build the "graphs." The RNAseq data will be a subset of the the same data loaded to GCP for the RNAseq hackathon. 

## GRCh38 Data


| Assembly Accession/Version | Path | Comment |
|----------|------------|----------|
| GCA_000001405.28 (GRCh38.p13) | [FTP link](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13) | The patches related data is in the follwing sub-directory: [FTP link](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_structure/) | It's not ready-made for analysis pipelines, but I will be encouraging teams to use the patch scaffolds in the graph building effort.|
|GCA_000001405.14 (GRCh38) | [FTP link](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines) | This data was probably in the RNAseq hackathon: can we confirm what we've got already? We want to be sure we've got the versions of the analysis sets that include the alternate loci, as that is considered content relevant for the graph building.|

## GOLD data in SRA
(Gold: https://www.genome.wustl.edu/items/reference-genome-improvement/)
