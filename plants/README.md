


#What about plants?

Use case​ : What about Plants?
Whereas many of the technologies being used or developed at the PanGenome Hackathon are
based on utilities primarily applied for the Human Genome; these technologies may be used for
other research organisms for study, such as other mammals, reptiles, insects, plants, and
microbes (viruses too?). Human PanGenome efforts appear designed for understanding of
genome architecture by way of structural variation; however, for the instance of plants these
technologies may also be designed for genetic change. Relationships of proximal and distal
phenotype-relevant events are key in understanding complex genome structure.


The combination of improved computational thresholds and sequencing technologies allowing
for different approaches to sequencing, and the use of hybrid-assembly methods to effectively
generate more full-coverage genomes has allowed the graph-based approaches to catalog
genome architecture to unprecedented levels. It is up to the graph-based research community
to guide structure in tool development and utility.


In recent years, several genomes of plants have been sequenced, and having such
technologies available, these have research and end-use applications for agriculture. Crop
plants form the foundation for the world's natural food and textile resources, and plant breeding
efforts are often focused on improving several quality traits. A graph-based sequence-centric
view of genomes set the stage for key decisions that can be made to improve the crop
infrastructure.


When considering the plant genome, there are quite an array of genome types with regard to
species identity, genome size, chromosome number, and ploidy level.
PanGenome studies have commenced on many of the model plant species for their attributes of
small genomes and short-timish research-time constraints, such as Arabidopsis thaliana
(flowering plants)(Clark et al., Science 2007), Medicago truncatula (legumes)(Zhou et al, 2017;
Miller et al, 2017), and Brachypodium distachyon (grasses)(Gordon et al., 2017). Likewise
studies have also initiated with pangenome studies on the larger cousin plants of crop economic
importance such as crucifers, soybean, and wheat (Montenegro et al, 2017), respectively.
These previous studies used highly developed sequence analyses, but were not necessarily
graph-based in approach. Several pangenome-related papers appear to be in the pipelines for
other important plant species; whether they use graph-based methods remain to be seen. The
exercise of testing graph-based sequence views will help formulate use-case scenarios.


Immediate interests would be to validate hypothetical evolutionary tree diagrams assigned to
species, and perhaps address instances where species are proposed to be ancient tetraploids,
or to compare genome changes in observed polyploid genomes. RNA-Seq methods may also
be matched against graph-based maps to quantify expression of the genomes; for instance, can
nutritional or medicinal trait changes be tracked to genome structural variation using
graph-based methods targeted on key metabolic pathway-associated genes. The tracking of
highly-repetitive transposons-initiated events may also explain some of the alterations observed
in different genome species and their evolutionary consequence resulting in gene duplication,
rearrangements, and the like. Using graph-based methods to map out highly variable regions
may also provide strategies toward implementing targeted engineering of species, or assist in
classic breeding strategies where known attributes are known to structurally exist. Similarly,
many wild ancestor lines are sought to bring in new gene function to serve as sources for
disease resistance, quality traits, and nutrition.


Visualization of genomes by graph-based methods appear to have a valuable role in the future
of agricultural improvements.

#References:

Clark et al. Science (2007) Arabidopsis...


Gordon et al. Nature Communications (2017) 8:2174.
Extensive gene content variation in the Brachypodium distachyon pan-genome correlates with population structure.
DOI: 10.1038/s41467-017-02292-8


Miller et al. BMC Genomics (2017) 18:541. Hybrid assembly with long and short reads
improves discovery of gene family expansions.
DOI 10.1186/s12864-017-3927-8


Montenegro et al. The Plant Journal (2017) 90:1007-1013.
The pangenome of hexaploid bread wheat.
DOI 10.1111/tpj.13515


Zhou et al. BMC Genomics (2017) 18:261
Exploring structural variation and gene family architecture with De Novo assemblies of 15 Medicago genomes.
DOI 10.1186/s12864-017-3654-1


#USE CASES

```
-Disease Resistance
-Evolution
-Role of Transposons in Genome Structure
-Inter-specifc Crosses / Hybridization
-Imputation
-Quality traits
-Food Quality Nutrition
-Metabolic Pathway Gene Tracking
-Modeling for Engineering (e.g. CRISPR-Cas9)
-Following Trait vs. Genotype
-Hyper-variable Genome Regions
```

#Graphs for Plants

Although a diploid, there is gene duplication to the extend that maize may be a hybridization of two genomes, an ancient allotetraploid. Looking at different plants that have different pathlines in evolution may help uncover unique biological mechanisms.


A sample graph alignment of sequenced maize genome was developed for Maize chromosome 10 P. Bradbury).


Application of genome graph maps to improve understandings of plant polyploid genomes. Data preparation was performed for the following. Due to some time limitations, this work will be moved forward post-hackathon. In plant breeding it is useful to have an understanding of the contributing genomes. Using the tools developed at the PanGenome Hackathon sample data from the wheat genome was used as a test example. The following should provide some insights for the utility of graph-based mappings:


The initial approach was to use sample reference maps for chromosome 1.

```
-Chromosome 1A Triticum aestivum cv. Chinese Spring
-Chromosome 1B Triticum aestivum cv. Chinese Spring
-Chromosome 1D Triticum aestivum cv. Chinese Spring
-Chromosome 1A Triticum turgidum cv. Svevo
-Chromosome 1B Triticum turgidum cv. Svevo
-Chromosome 1D Aegilops tauschii
```

Other possibility for the grass species would be to use the model plant Brachypodium spp.; over 40 species have been sequenced and assembled (Gordon et al., 2017). It would be of interest how this species as a model organism can make connections to a distant relative within an evolutionary tree; other of interest would be to compare a distant species such as Avena sativa (a distant grass), or a more closely associated species within the Triticeae, such as Hordeum vulgare, or other Triticeae species.


Applications:


Understanding evolution by comparing ancestor lines to cultivated lines. Focus on key regions may validate hypotheses developed through genome surveys. As an example is the storage proteins of wheat where...


Studies continue.
