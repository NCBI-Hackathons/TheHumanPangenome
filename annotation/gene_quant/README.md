# Quantify RNAseq reads per gene

A common pipeline that researchers perform for differential expression analysis is mapping RNAseq reads to the graph and quantifying their counts per gene.
A useful graphical pangenome would also be able to perform this functionality while also enabling extra features like allele-specific transcription profiling.

## Pipeline

1. Construct a reference graph from chromosome 21 of GRCh38 and variants from the 1000 Genomes Project
2. Convert the graph to a splice-aware graph using `vg rna` and index it (running)
3. Map RNAseq reads to the graph (`vg map`)
4. Count coverage with `vg pack`
5. Quantify coverage per gene
