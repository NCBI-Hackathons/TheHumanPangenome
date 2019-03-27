# vg_wdl

Mike Lin @mlin and Charlie Markello @cmarkello hacked on the [vg_wdl](https://github.com/vgteam/vg_wdl) repository, which (while incipient) is accumulating vg-related docker+[WDL](http://openwdl.org/) workflows and automatically testing them on small examples. 

Here's how you can build a vg genome graph, complete with GBWT index of 1000 Genomes haplotypes on small variants, for the *ABO* locus which influences blood type. System requirements: git, python3-pip, Java 11 JRE, docker.

```
pip3 install miniwdl
git clone git@github.com:vgteam/vg_wdl.git
miniwdl cromwell --no-quant-check vg_wdl/workflows/vg_construct_and_index.wdl \
    graph_name=ABOlocus contigs=ABOlocus use_haplotypes=true \
    ref_fasta_gz=vg_wdl/tests/ABOlocus/ABOlocus.fa.gz \
    contigs_vcf_gz=tests/ABOlocus/ABOlocus_small.vcf.gz
```

Beyond SNVs/indels, there are several interesting ABO structural variants standing in the human population, which give rise to unusual blood histocompatibility phenotypes. We selected a known 3.8Kbp deletion (relevant publication: [doi:10.1111/vox.12613](https://onlinelibrary.wiley.com/doi/full/10.1111/vox.12613)), which can be modeled in a vg graph:

```
miniwdl cromwell --no-quant-check vg_wdl/workflows/vg_construct_and_index.wdl \
    graph_name=ABOlocus_SV contigs=ABOlocus \
    ref_fasta_gz=vg_wdl/tests/ABOlocus/ABOlocus.fa.gz \
    contigs_vcf_gz=tests/ABOlocus/ABOlocus_SV.vcf.gz
```

A few 1000 Genomes samples exhibit this deletion. The [vg_wdl continuous integration tests](https://github.com/vgteam/vg_wdl/blob/master/tests/ABOlocus/vg_ABOlocus_test_SV.wdl) map the reads from one such sample to this graph, and confirm that vg aligns some reads along the deletion edge.
