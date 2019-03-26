# Team: Building haplotypes and graphs from full de novo assemblies
Our team focus on tow key question: 1) investigate what is the best way to build genome graphs from de novo assemblies, 2) define a graph coordinate system that will be retro-compatible with GRCh38.

Slides: https://docs.google.com/presentation/d/1mocR6xLbRdgga79rXrkVModWd7YZiVCvTMXyQHhlBqU/edit?usp=sharing

## AIM 1 - Build graphs from de *novo assemblies*
To explore the best strategy to build genome graphs from *de novo* assembled genomes, we designed two experiments. The first one is a simple proof of concept graph build using the Japanese genome assembly jg1 and GRCh38, and the second graph is a great-Apes graph using GRCh38, the Chimp (Clint), and Orangutan (Susie). We are comparing the graphs built using three different strategies: 

1) [Sequish](https://github.com/ekg/seqwish) [in progress]
We created a graph based on sequences from chromosome 21 from GRCh38 (CM000683.2), Clint the Chimpanzee (CM009259.2), Susie the Sumatran orangutan (CM009283.2), and CHM1 (AC244111.3, AC244144.2, AC244518.2, AC245051.3, AC245314.2, AC246819.2, AC255431.1, AC256301.1, AC277730.1, AC277802.1, AC277887.1).  We used first used minimap2 (v2.16-r922) with the parameter preset asm5 to do an all-vs-all alignment of the sequences. Next, we used seqwish (6e4fe705) to induce a graph in GFAv1 format.

We use vg view and index to convert the GFAv1 graph into vg format and build the different required indexes (Fig. 1).

(insert pipeline here)

2) [SibeliaZ](https://github.com/medvedevgroup/SibeliaZ) [in progress]
We are creating a graph using chromosome 1 of the Japanese *de novo* assembly jg1 and GRCh38. We use SibeliaZ a fast 

3) [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). [TODO]

## AIM 2 - 

