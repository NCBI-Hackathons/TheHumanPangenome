# Building haplotypes and graphs from full *de novo* assemblies

The goal of our team is to establish pipelines to build graph genomes from *de novo* assembly level data, as well as define a new coordinate system that will be retro-compatible with GRCh38 to enhance interoperability.

![panApe](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/ds-vis1.png?raw=true "The path of GRCh38 through a graph")

**The path of GRCh38 through a graph**


## Projects
### 1) define the best strategy to build genome graphs from *de novo* assemblies.
To explore the best strategy to build genome graphs from *de novo* assembled genomes, we designed two experiments. The first one is a simple proof of concept graph build using the Japanese genome assembly jg1 and GRCh38 (chromosome 21), and the second graph is a great-Apes graph using GRCh38, the Chimp (Clint), and Orangutan (Susie) (Also using chromosome 21). We are comparing the graphs built using three different strategies: 

![flowchart](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/graph2xg_all.png?raw=true)

1) [Seqwish](https://github.com/ekg/seqwish)
We created a graph based on sequences from chromosome 21 from GRCh38 (CM000683.2), Clint the Chimpanzee (CM009259.2), Susie the Sumatran orangutan (CM009283.2), and CHM1 (AC244111.3, AC244144.2, AC244518.2, AC245051.3, AC245314.2, AC246819.2, AC255431.1, AC256301.1, AC277730.1, AC277802.1, AC277887.1).  We used first used minimap2 (v2.16-r922) with the parameter preset asm5 to do an all-vs-all alignment of the sequences. Next, we used seqwish (6e4fe705) to induce a graph in GFAv1 format.  See [script](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/scripts/seqwish.sh) We use vg view and index to convert the GFAv1 graph into vg format and build the different required indexes (See flowchart above).

```
# run the seqwish pipeline
./scripts/seqwish.sh

```

2) [SibeliaZ](https://github.com/medvedevgroup/SibeliaZ)
We created a graph from chromosomes 1 of the Japanese *de novo* assembly jg1 and GRCh38.
To do so we first ran SibeliaZ to produce a whole-genome alignment of these two sequences in MAF format.
Then we converted MAF to GFA1 and imported it into vg.
To make the conversion, we wrote a script producing a GFA1 graph from a MAF alignment (scripts/maf_to_gfa1.py).

An updated version of the script along with exact commands on how to generate the GFA file will be included in the next minor release of SibeliaZ.

```
# construct the graph
sibeliaz-lcb --graph chr1_hg38_jg1_deBruijngraph.dbg --fasta hg38_chr1.fasta -k 25 -b 200 -o chr1_gddd -m 500 -t 4 --abundance 16

# convert MAF to GFAv1
scripts/maf_to_gfa1.py <in.maf> <out.gfa1>

# convert GFAv1 to vg
vg view <in.gfa1> <out.vg>
```

3) [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus)
Cactus allows the definition of an evolutionary tree and uses the phylogenetic relationship between the samples in the graph construction process. We build a small great apes pangenome (panApe) for chromosome 21 using the tree defined bellow (figure). 

![cactus_tree](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/cactus_panApe_tree.png?raw=true "Apes Tree")

```
# Generate the HAL graph using cactus
cactus jobStore ./panApe21_cactus.txt panApe_chr21.hal --buildHal

# Convert HAL to vg using hal2vg
hal2vg panApe_chr21.hal > panApe_chr21.vg
```

#### Graph visualisation 
Following the conversion we loaded the graph onto Bandage for visualisation. Another alternative is Momi-G (https://github.com/MoMI-G/MoMI-G).

| Example 1 | Example 2 |
| - | - |
| ![bandage1](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/ds-vis1.png?raw=true "panApe21_1") | ![bandage2](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/ds-vis2.png?raw=true "panApe21_2") |

### 2) define a graph coordinate system that will be retro-compatible with GRCh38.

Representing haplotype information in reference genome is beneficial in increasing mappability and reducing bias. The major concerns for representing haplotype in the existing reference genome are the alteration of coordinate, redundant representation, and ambiguous sequence inference. Our proposed notation tackle these issues with the following design philosophies:

1) The haplotype contigs are coordinated and defined as an add-on outside the reference genome coordinate. This allow the haplotype contigs set to be updated separately and the inclusion of haplotype sequence does not alter genomic coordinate. This design also allow the user to use patch in reference genome or recreating custom sequence using their haplotype of interest.

2) Each haplotype and nested haplotype are defined as a unique segment based on the reference genome or the closest haplotype; therefore, the number of bases that need to be stored for each haplotype sequence is minimalized.

3) Each haplotype can be uniquely represent using Graphical Fragment Assembly-liked notation which can track back into node storing specific sequence for each haplotype.

![Proposed graph coordinate system](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/figure1.png)

Proposed graph coordinate system to represent haplotype. A) GFA file that includes nodes and edges for reference genome and haplotypes. B) Accompany path file include path for reference genome, haplotig 1, and haplotig 2. These paths have corresponding sequences stored separately. C) Visualization of A using path label from B. An individual with haplotype 1 should be compared against sequence 1-3-5-6-4. An individual with type 2 should be compared against sequence 1-3-7-6-4. 

![Updating graph](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/figure2.png)

Expanding the model Our proposed model allows nodes and edges representation of the GFA to change while the sequence corresponding to each haplotive can remain unchanged 

![coordinates](https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/master/DS/graph/Genome_graph_compare.png?raw=true "Apes Tree")

## Useful links

- Choosing which vg index to use: https://github.com/vgteam/vg/wiki/Index-Construction

- Example vg graphs and indexes: https://github.com/vgteam/sequenceTubeMap/tree/master/exampleData

- Slides: https://docs.google.com/presentation/d/1mocR6xLbRdgga79rXrkVModWd7YZiVCvTMXyQHhlBqU/edit?usp=sharing

## TODOs:

- Write code to browse the graph and compute summary statistics.

- Finish Cactus runs.

- Map Japanese short reads and quantify the improvements

- Compute Fst, Pi, Theta, ...

## ISSUES: 

- Seqwish doesn't output correct CIGAR strings.


## Contributors
|Name         |Email        |
| ----------- | ----------- |
| Benedict Paten | bpaten@ucsc.edu |
| Arkarachai Fungtammasan | arkarachai.af@gmail.com |
| Medhat Mahmoud | helmy.medhat@gmail.com |
| Jana Ebler |jebler@mpi-inf.mpg.de |
| Yassine Souilmi | yassine.souilmi@adelaide.edu.au |
| Valerie Schneider | schneiva@ncbi.nlm.nih.gov |
| Zev Kronenberg | zkronenberg@pacificbiosciences.com |
| Allison Regier | aregier@wustl.edu |
| Ilia Minkin | ivminkin@gmail.com |
