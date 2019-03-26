# Annotation in reference genome graphs

## Motivation

Current methods of annotating genomes rely on genomic intervals as a core formalism.
There are some difficulties in generalizing this formalism to reference graphs.
A genomic interval corresponds to a path in the graph.
However, if we restrict the annotation to one path in the graph, the alternate alleles included in the graph are not included in the annotation.
We argue that connected subgraphs are a more appropriate formalism for genome graphs.

Using a new core formalism for annotation necessarily means that there is not currently infrastructure built out to support it.
We need exchangeable representations of the data, software support, and analysis tools to make the formalism useful for practitioners.
This is what we are proposing to develop during this hackathon.

## Project

We are proposing to develop a proof-of-concept system for utilizing continuous-valued annotations with genome graphs. Our goals are to...

* Agree on a basic data representation for annotating subgraphs with continuous valued data
* Specify a human-readable file format for expressing the data representation
* Produce a demonstration set of annotation data for a genome graph
* Augment VG to produce and work with these data
* Hook the file format into an existing graph visualization tool

## Updates

### Day 2

#### Annotation Import - Easy Case

![Annotation Import, The Easy Case](fig/annotation_easy_case.svg)

#### Annotation Import - Harder Cases

![Annotation Import, The Harder Cases](fig/annotation_harder_cases.svg)

#### Gene-level RNAseq quantification pipeline

See the subproject-specific [README](gene_quant/README.md).

#### Visualization using MoMIG

![MoMIG Genome Graph Visualization software screenshot](fig/momig_screenshot.png)
