
## Specific aim:  A Faster, Better Short-Read Mapper with Hit Chaining.

Our work modifies [this function](https://github.com/vgteam/vg/blob/master/src/subcommand/gaffe_main.cpp) in `vg`.

#### A BETTER Short-Read Mapper
This allows us to deal with cross-overs and indels.  Anything that we can't deal with on the first pass.

The current behavior is to ignore indels and cross-overs and only provide gapless alignments (the default allows for four mismatched bases).

#### A FASTER Short-Read Mapper
When creating an alignment from numerous graph pathways, the speed of alignment is slowed exponentially with each node present in "snarl" regions.  These "snarl" regions constitute a bubble on the graph where multiple nodes can be chosen.  Each of these nodes can exponentially increase the number of paths within the "snarl", and thus the alignment time.  We attempt to improve upon this by excluding nodes not associated with a given haplotype.

Our algorithm speeds up the clustering after the search hits are gathered.

## Running the Code

You will need to run `vg` to see the updated function at work.

We'll be using `valgrind` which can be installed with: `sudo apt install valgrind`

Next, clone `vg` with: `git clone https://github.com/vgteam/vg.git`

At the top of [src/subcommand/gaffe_main.cpp](https://github.com/vgteam/vg/blob/master/src/subcommand/gaffe_main.cpp) please uncomment this line to use with valgrind: https://github.com/vgteam/vg/blob/master/src/subcommand/gaffe_main.cpp#L25

Then follow [the installation instructions at the vg repo to compile the repo from source](https://github.com/vgteam/vg) (warning: this can take 30 minutes+ to compile the first time).

Data files for testing are [located here](http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing).  They can be downloaded with:
```
#!/usr/bin/env bash

mkdir data
mkdir data/reads
curl -o data/snp1kg-CHR21_filter.dist http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.dist
curl -o data/snp1kg-CHR21_filter.gbwt http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.gbwt
curl -o data/snp1kg-CHR21_filter.gcsa http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.gcsa
curl -o data/snp1kg-CHR21_filter.gcsa.lcp http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.gcsa.lcp
curl -o data/snp1kg-CHR21_filter.min http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.min
curl -o data/snp1kg-CHR21_filter.snarls http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.snarls
curl -o data/snp1kg-CHR21_filter.vg http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.vg
curl -o data/snp1kg-CHR21_filter.xg http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.xg
curl -o data/reads/sim.fq.gz http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/sim.fq.gz
curl -o data/reads/sim.gam http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/sim.gam
curl -o data/reads/toil-vg-sim.txt http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/toil-vg-sim.txt
curl -o data/reads/true.pos http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/true.pos
```

Now that `vg` is installed, you can run the following using the sample data above:
```
#!/usr/bin/env bash

valgrind --tool=callgrind --instr-atstart=no \
    ./bin/vg gaffe \
    -x data/snp1kg-CHR21_filter.xg \
    -m data/snp1kg-CHR21_filter.min \
    -d data/snp1kg-CHR21_filter.dist \
    -s data/snp1kg-CHR21_filter.snarls \
    -H data/snp1kg-CHR21_filter.gbwt \
    -G data/reads/sim.gam \
    > data/mapped.gam
```

Then to annotate, run:

```
valgrind --tool=callgrind --instr-atstart=no \
     ./bin/vg annotate -p -x data/snp1kg-CHR21_filter.xg -a data/mapped.gam > data/annotated.gam
```

Finally, to compare the annotated reads to the truth and set the mapped_correctly field:

```
valgrind --tool=callgrind --instr-atstart=no \
    ./bin/vg gamcompare  -r 100 data/annotated.gam data/reads/sim.gam > data/compared.gam
```

## Slides
![Image00](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/00.png)
![Image01](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/01.png)
![Image02](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/02.png)
![Image03](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/03.png)
![Image04](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/04.png)

## Stretch Goals:
- Fast Clustering (Stretch Goal #1)
- Haplotype-based Hit Joining (Stretch Goal #2)

## Fallback Goal:
Dump hit coverage by node.
