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

#### The gGFF format

We have defined a generalization of the GFF3 format that replaces genomic intervals with a subraph. It is a text-based, tab-separated file. Every line contains each of the following fields. If a field is to be ignored, it can be replaced with a "." (without quotes). The fields are

* subgraph: a comma separated list of intervals of sequences on nodes, along with orientation in the format `ID[start:end](+/-/?)`.
* source: the name of the program or database that generated the annotation
* type: the type of feature
* score: a floating point value
* phase: 0, 1, or 2 indicating the first base of the feature that is a codon, measuring from the source node in the subgraph
* attributes: a semi-colon separated list of tag-value pairs, with tags separated from the values by an "="

## Team Members

* Travis Wrightsman (tw493@cornell.edu)
* Jordan Eizenga (jeizenga@ucsc.edu)
* Rajeeva Musunuri (rmusunuri@nygenome.org)
* Toshiyuki Yokoyama (toshiyuki.t.yokoyama@gmail.com)

## Future Directions

* local alignment of nodes could help resolve annotations that span outside of snarls
  * both distance or best alignment methods could be used to traverse the graph and liftover annotation from the reference path to other paths
