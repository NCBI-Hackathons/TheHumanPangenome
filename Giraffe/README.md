
## Specific aim:  A Faster, Better Short-Read Mapper with Hit Chaining.

This allows us to deal with cross-overs and indels.  Anything that we can't deal with on the first pass.

The current behavior is to ignore indels and cross-overs and only provide gapless alignments (the default allows for four mismatched bases).

Our work modifies [this function](https://github.com/vgteam/vg/blob/master/src/subcommand/gaffe_main.cpp) in `vg`.

## Running the Code

You will need to run `vg` to see the updated function at work.

Follow [the installation instructions at the vg repo](https://github.com/vgteam/vg) (or install a precompiled binary since it can take 30 minutes+ to compile vg).

Data files for testing are [located here](http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing).  They can be downloaded with:
```
curl -o snp1kg-CHR21_filter.dist http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.dist
curl -o snp1kg-CHR21_filter.gbwt http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.gbwt
curl -o snp1kg-CHR21_filter.gcsa http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.gcsa
curl -o snp1kg-CHR21_filter.gcsa.lcp http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.gcsa.lcp
curl -o snp1kg-CHR21_filter.min http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.min
curl -o snp1kg-CHR21_filter.snarls http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.snarls
curl -o snp1kg-CHR21_filter.vg http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.vg
curl -o snp1kg-CHR21_filter.xg http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/snp1kg-CHR21_filter.xg
curl -o reads/sim.fq.gz http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/sim.fq.gz
curl -o reads/sim.gam http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/sim.gam
curl -o reads/toil-vg-sim.txt http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/toil-vg-sim.txt
curl -o reads/true.pos http://public.gi.ucsc.edu/~anovak/graphs/gaffe/basic-testing/reads/true.pos
```

Now that `vg` is installed, you can run the following using the sample data above:
```
vg gaffe \
    -x snp1kg-CHR21_filter.xg \
    -m snp1kg-CHR21_filter.min \
    -d snp1kg-CHR21_filter.dist \
    -s snp1kg-CHR21_filter.snarls \
    -H snp1kg-CHR21_filter.gbwt \
    -G reads/sim.gam \
    > mapped.gam
```

## Slides
![Image00](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/00.png)
![Image01](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/01.png)
![Image02](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/02.png)
![Image03](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/03.png)
![Image04](https://raw.githubusercontent.com/NCBI-Hackathons/TheHumanPangenome/master/Giraffe/images/04.png)

## Stretch Goals:
#### Fast Clustering (Stretch Goal #1)
When creating an alignment from numerous graph pathways, the speed of alignment is slowed exponentially with each node present in "snarl" regions.  These "snarl" regions constitute a bubble on the graph where multiple nodes can be chosen.  Each of these nodes can exponentially increase the number of paths within the "snarl", and thus the alignment time.  We attempt to improve upon this by excluding nodes not associated with a given haplotype.

Our algorithm speeds up the clustering after the search hits are gathered.

#### Haplotype-based Hit Joining (Stretch Goal #2)

## Fallback Goal:
Dump hit coverage by node.
