## Motivation

* Association of enhancers and promoter regions with specifically expressed alleles

* Unambiguous elucidation of both cis- and trans- eQTLs


## Aims

![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA_project_day3_aims.jpg)



## Workflow

### 1. Contruct spliced variant graphs

#### 1a. Single whole genome graph

#### 1b. Multiple disjoint gene (splice) graphs

### 2. Map short reads to graph(s)

### 3. Surject mapped reads to haplotype-specific transcripts

### 4. Run RSEM or similar

![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA_project_day_3_pipeline.jpg)


## Implementation

### Download the relevant pipeline scripts
```
wget https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/RNA/wdl_pipeline/vg_rna.wdl -O vg_rna.wdl
wget https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/RNA/wdl_pipeline/vg_rna.inputs.json -O vg_rna.inputs.json
```

### Download cromwell (requires java 8 and Docker to run)
#### Quickstart intro on using Cromwell: (https://cromwell.readthedocs.io/en/develop/tutorials/FiveMinuteIntro/)
```
wget https://github.com/broadinstitute/cromwell/releases/download/36.1/cromwell-36.1.jar .
```
### Download example input files
```
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz .
bgzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

wget ...
```
### Run the WDL workflow
```
sudo java -jar cromwell-36.1.jar run vg_rna.wdl -i vg_rna.inputs.json

```

## Results

![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA_project_day3_SMIM11A.jpg)
![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA_project_day3_SMIM11A_var.jpg)
![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA_project_day3_mapstats.jpg)
***Note that the numbers shown above for the simulated data (first two rows) are for alignments that have not been filtered according to mapping quality**  

## Future goals

* Complete allele-specific expression pipeline

* Comprehensive evaluation of splice graph mapping using vg




