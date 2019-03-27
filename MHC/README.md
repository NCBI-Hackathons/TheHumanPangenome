Resolve MHC regions from HG002
===============================


General Goal During the hackathon
=================================
Focus on developing diploid assemblies of MHC in HG002 – these could be only the third and fourth authoritative 
haplotype assemblies produced for MHC!

Build a graph from these haplotypes and probably haplotypes from the population as well

Potentially use 10x genomics linked-reads to refine this graph by “calling variants on the graph” in HG002 (e.g., to fix homopolymers)

Potentially project graph onto primary reference and produce vcfs for GRCh37 and 38. 


Data preparation
================

* Fetch reads from PacBio HG002 CCS data mapped to MHC regions
    * recruit read mapping to GRCh37 MHC region (minimap2?), 
    * recruit read mapping to assembled contig of the MHC
    * recruit pre-phased reads from Whatshap (using pre-mapped BAM files)

* Fetch reads from ONT mapped to MHC region
    * using pre-alignment data from GIAB for phasing SVs

* Fetch reads from 10x genomics to MHC region
    * using pre-alignment data for phasing SNPs to connect phasing blocks 

Assembly approach
=================
* Using reference:
    * bin reads with their concordance with the reference
    * de novo assembling discordant reads
    * map assembled contig to ONT reads to scafflods, and/or phase specific scaffold by large SVs

* De Novo:
    * Approach from Shilpa's recent publication (see Tool Integration too)
    * FALCON-Unzip style first + other information to extend haplotype
        * directly layout contigs using string graph for both haplotype
    * Binning reads by haplotype using phasing data 
        * Just used phased reads / phased blocks from HG002 -> de novo assembly

* Using parental reads/kmers to bin reads

* Using 10x phasing information to chain disconnect phasing blocks

Tool Integration
================
* Reads to Graph mapping (Adam, Tobias, and Shilpa)
    * What kind graph, can use this to resolve haplotypes

* Quality evaluation and (Structure) variant calling:
    * assembly-assembly alignment to VCF?
    * compare to current GIAB call set?

Visualization
=============
* Read mapping approach

* Assemblytic
            
* MHC comparison between different individual:
    - Simulate reads from GRCh38 alts of MHC -> assembly
    - What other dataset allows us to reconstruct MHC fast? Maybe get some assembled genome from the NCBI_Data_Upload


After the hackathon
===================

* Explore ways these analyses could be used to benchmark methods used to characterize MHC

* Can these haplotypes be represented in standard VCF with respect to the primary GRCh37/38 reference in GIAB benchmark sets?

* Or do we need new representations and benchmarking tools?

* vg can project haplotypes into a vcf with respect to the primary reference, but we’ll have to see if this is compatible with current benchmarking tools for small variants and structural variants

* Probably hold off on KIR, since it has more copy number variation




Presentation
==================
https://docs.google.com/presentation/d/1cueyPgZLsPW1Fyimu5dt3OgsTj3KM0xL3h6fKJvWcXI/edit?usp=sharing




