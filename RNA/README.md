# Motivation

## Association of Enhancers and Promoter Regions with Specifically Expressed Alleles

## Unambiguous elucidation of both cis- and trans- eQTLs

# Summary

## 1. Create proof of concept pipeline for estimating allele-specific expression

## 2. Verify Phasing in Genic Regions

## 3. Wrap it all in WDL

# Workflow

![Alt text](?raw=true "Title")

## 1. Contruct spliced variant graphs from HG002, Chromosome 21

### 1a. Single whole genome graph

### 1b. Multiple disjoint gene graph

## 2. Map Short Reads 

## 3. Run RSEM or similar (maybe HISAT2-gene)

## 4. Map Isoseq

## 5. Eval Phasing with Isoseq

## 6. Compare 'RSEM' with Isoseq

## 7. Define "Haplotype-Guided Novel Transcripts" (stretch)

### 7a. Possibly integrate: https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/tree/master/Splice

## 8. Wrap it all in WDL!

# Results

# Presentations

### Day 2 presentation
https://docs.google.com/presentation/d/1wWK4d39ZiJW09IytHxArcYOA3xTiF3rgxKcM5b62BmA/edit?usp=sharing


