
## Specific aim:  A Faster, Better Short-Read Mapper with Hit Chaining.

Our work modifies [this function](https://github.com/vgteam/vg/blob/master/src/subcommand/gaffe_main.cpp) in `vg`.

Improve a prototype minimizer-based mapper, by adding a faster clustering function to cluster minimizer hits, and hit extension logic for handling clusters that have no good full-length gapless alignment.

The clustering algorithm has been improved by reducing the amount of data copying in the implementation. Additionally, we have devised an improved algorithm for comparing sets of clusters.

## Running the Code

You will need to run `vg` to see the updated function at work.

We'll be using `valgrind` which can be installed with: `sudo apt install valgrind kcachegrind`

At the time of writing, the current vg `master` branch did not contain the updated gaffe function and was cloned with: `git clone --recursive https://github.com/vgteam/vg.git`

A repo repo was cloned and compiled with the new gaffe function using:
```
git clone --recursive https://github.com/vgteam/vg.git && cd vg
git pull https://github.com/xchang1/vg.git seed_clustering
```

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

Optionally, one can run:
```
./bin/vg annotate -p -x data/snp1kg-CHR21_filter.xg -a data/mapped.gam > data/annotated.gam
```

Finally, to compare the annotated reads to the truth and set the mapped_correctly field:
```
./bin/vg gamcompare -r 100 data/annotated.gam data/reads/sim.gam > data/compared.gam
```

## Final Slides

### Clusterer Speedup
 * The clustere is not slower.
 * The clusterer is slightly faster
 
### Maize Graph
 * A maize graph has been built.
 * It contains a snarl with ~107,000 nodes, which is too big for the distance index to handle.
 * `odgi` draws this beautiful picture: ![ODGI maize graph hairball](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/odgi.png)
 

## Day 1 Slides
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
