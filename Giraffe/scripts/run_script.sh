#!/usr/bin/env bash

sudo apt install valgrind kcachegrind

git clone --recursive https://github.com/vgteam/vg.git && cd vg && make get-deps
sed -i 's/\/\/#define USE_CALLGRIND/#define USE_CALLGRIND/g' src/subcommand/gaffe_main.cpp

. ./source_me.sh && make

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

valgrind --tool=callgrind --instr-atstart=no \
    ./bin/vg gaffe \
    -x data/snp1kg-CHR21_filter.xg \
    -m data/snp1kg-CHR21_filter.min \
    -d data/snp1kg-CHR21_filter.dist \
    -s data/snp1kg-CHR21_filter.snarls \
    -H data/snp1kg-CHR21_filter.gbwt \
    -G data/reads/sim.gam \
    > data/mapped.gam
