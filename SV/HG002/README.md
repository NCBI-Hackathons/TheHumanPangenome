# HG002 comparison

We want to compare variant calls using three types of sequencing data:

1. Illumina short reads.
1. CCS long reads.
1. PacBio long reads.

## Mapping

Using vg map. For long reads also [GraphAligner](https://github.com/maickrau/GraphAligner).

## Comparing variant calls

### Re-injecting path information in the graph

If we had the SVs, e.g. insertions, annotated as paths in the graph we could directly get variant calls relative to this path.

### Comparing the inserted sequence

If we don't have the SVs annotated as path, vg outputs the inserted sequence in the ALT.
To get internal SNPs and compare datasets we will have to compare the inserted sequences.

Code in `callVariantsInInsertedSeq.R`.

#### Cluster insertions

Before comparing the sequences we must decide which insertions to compare. 
To be considered the same insertions, two insertions should be:

1. Close together: e.g. distance between insertion location <30 bp.
1. More or less the same size

We could also catch wrongly assigned pairs of insertions when the sequences don't align well in the next step.

This kind of operations were already implemented to compare SVs in [sveval](https://github.com/jmonlong/sveval).

#### Smith-Waterman

The inserted sequence called and in the original VCF are aligned to detect SNVs and indels.


## All reads vs regionally downsampled reads

See [report](hg002-allVsRegionalDownsampling.md)
