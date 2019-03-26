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

## Slides

###Day 2

https://docs.google.com/presentation/d/1UqPvWb7xUCGaKIv6qgX95_uLZxRAZOzXvDO4IkpGn6U/edit?usp=sharing
