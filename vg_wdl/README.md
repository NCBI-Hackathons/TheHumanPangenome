# vg_wdl

Mike Lin and Charlie Markello hacked on the [vg_wdl](https://github.com/vgteam/vg_wdl) repository, which (while incipient) is accumulating vg-related docker+[WDL](http://openwdl.org/) workflows and automatically testing them on small examples. 

Here's how you can build a vg genome graph, complete with GBWT index of 1000 Genomes haplotypes on small variants, for the *ABO* locus which influences blood type. This should run in a few minutes on your laptop, given system requirements: git, python3-pip, Java 11 JRE, docker set up for [non-root](https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user) invocation.

```
pip3 install miniwdl
git clone git@github.com:vgteam/vg_wdl.git
miniwdl cromwell vg_wdl/workflows/vg_construct_and_index.wdl \
    graph_name=ABOlocus contigs=ABOlocus use_haplotypes=true \
    ref_fasta_gz=vg_wdl/tests/ABOlocus/ABOlocus.fa.gz \
    contigs_vcf_gz=vg_wdl/tests/ABOlocus/ABOlocus_small.vcf.gz
```

The [vg_construct_and_index.wdl](https://github.com/vgteam/vg_wdl/blob/master/workflows/vg_construct_and_index.wdl) workflow puts [prose instructions from the vg wiki](https://github.com/vgteam/vg/wiki/Index-Construction) into portably executable form. The next step will be to feed the graph+indices produced from this into [Charlie's mapping+calling workflow](https://github.com/vgteam/vg_wdl/blob/8eea6a9dd078e8110cb2e12cea1748fdbfd6b3e0/workflows/vg_pipeline.workingexample.wdl).

Beyond SNVs/indels, there are several interesting ABO structural variants standing in human populations, which give rise to unusual blood histocompatibility phenotypes. We selected a known 3.8Kbp deletion (relevant publication: [doi:10.1111/vox.12613](https://onlinelibrary.wiley.com/doi/full/10.1111/vox.12613)), which can be modeled in a vg graph:

```
miniwdl cromwell vg_wdl/workflows/vg_construct_and_index.wdl \
    graph_name=ABOlocus_SV contigs=ABOlocus \
    ref_fasta_gz=vg_wdl/tests/ABOlocus/ABOlocus.fa.gz \
    contigs_vcf_gz=vg_wdl/tests/ABOlocus/ABOlocus_SV.vcf.gz
```

A few 1000 Genomes samples carry this deletion. The [vg_wdl tests](https://github.com/vgteam/vg_wdl/blob/master/tests/ABOlocus/vg_ABOlocus_test_SV.wdl) map the reads from one such carrier (HG01308) to this graph, and confirm that vg aligns some reads along the deletion edge.
