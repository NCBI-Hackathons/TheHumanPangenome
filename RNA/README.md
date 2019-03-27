# Motivation

## Association of Enhancers and Promoter Regions with Specifically Expressed Alleles

## Unambiguous elucidation of both cis- and trans- eQTLs

# Summary

## 1. Create proof of concept pipeline for estimating allele-specific expression

## 2. Wrap it all in WDL

![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA_project_day3_aims.jpg)

# Implementation

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


# Workflow

## 1. Contruct spliced variant graphs from HG002, Chromosome 21

### 1a. Single whole genome graph

### 1b. Multiple disjoint gene graph

## 2. Map Short Reads 

## 3. Run RSEM or similar

## 4. Wrap it all in WDL!

![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA_project_day3_pipeline.jpg)

# Results

![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA%20project-2.2.jpg?raw=true "Title")

![Alt text](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/RNA/RNA%20project%20-%20day%203-2.jpg?raw=true "Title")





