# Pangenomics Hackathon Team Focus Questions and User Stories
## Graph Building (from the GRC/NCBI, at this point, please add/modify)
  * User story: As the GRC, we can curate the pan-genome, so that it provides the most up-to-date and accurate data representation 
    * What is the mechanism for adding (or removing) an entire assembly from the graph? 
      [] Can this be done incrementally, or must everything start from scratch?
    * Same question, but for a partial assembly?
    * If a portion of the graph is determined to be wrong, how will that be communicated (e.g. so that tools do not use that portion of the graph)
## QC
  * As a graph builder, I have quality metrics for the pan-genome representation 
    * e.g. path errors, sequence errors
    * support for bases, SVs, haplotypes (paths in general). 
  * As a graph builder, I can tell whether my changes to the graph have improved or woresened the pan-genome representation
## Tools
  * As a researcher, I can convert/remap between the linear coordinates from two or more different assemblies that are consituents of the graph 
    * Consider the graph as an interpreter between linear assemblies...
    * Or from GRCh37 into the graph coordinates (if I'm a clinician)
  * As a researcher, I can generate a linear representation of an assembly from the pan-genome (e.g. GRCh38 path, or some other path) 
    * note: could just use assembly fasta, but would be cool to build something that can reconstruct those sequences from the graph itself
  * As the GRC, I can provide mark-up for regions of the pan-genome, so we can tell users about regions of interest
  * As a user, I can find matches to my sequence of interest in the pan-genome/graph 
    * best match, other matches
    * meta data that tells me in which populations/samples those paths exist
    * coordinates within the pan-genome or a specific consituent assembly

## Visualization
  * As a clinical researcher, I can tell which parts of the graph are relevant to the consituent assemblies, so I can tell whether its likely that my sample will share those characteristics
  * As a population geneticist, I can visualize a sample/population/haplotype of interest within the graph
  * As a user/GRC, I can visualize regions of uncertainty/under curation in the graph  
    * uncertainty refers in this case to potential sequencing errors, chimerism or other contig assembly construction errors
  * As a researcher with a genomic region of interest, I can visualize the wide variety of data generated on the GRCh38 linear reference in the context the pan-genome (genes, epi markers, RNAseq, etc) 
    * per consituent assembly?
    * on the graph?
  * As the GRC, we can tell which assemblies contributed to each part of the graph, so we know how their curation is likely to impact the graph/pan-genome
