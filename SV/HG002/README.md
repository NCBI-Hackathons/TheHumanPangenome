# HG002 comparison

We want to compare variant calls using three types of sequencing data:

1. Illumina short reads.
1. CCS long reads.
1. PacBio long reads.

For the short reads we also want to compare the results using:

1. All the reads.
1. Reads downsampled around regions of interest.

## Mapping

Using vg map. For long reads also [GraphAligner](https://github.com/maickrau/GraphAligner).

## Identifying internal variants in insertions

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

The inserted sequence called and the inserted sequence in the original VCF are aligned to detect SNVs and indels.

For example from an alignment like below we identify an insertion of `T`, a SNV (`C>T`) and a deletion of `GT`:

```
Original insertion: AGCTAGGCTA-TTACGACGATCGAGTGCTATCGACTAGCATCG--
                                     x
Called insertion:   --CTAGGCTATTTACGATGATCGA--GCTATCGACTAGCATCGAT
```

## Comparison results

## All reads vs regionally downsampled reads

Variants were called using the full set of reads from the 30x and using the regionally downsampled reads.
We want to know if the variants identified are similar so that we could use regional downsampling to speed up the analysis on more samples.

TL;DR: 

- We found important differences, including many variants unique to the downsampled analysis.
- Some might differences might be due to some insertions being called only when all reads are used (overlap increased when focusing on variants in insertions called in both).
- Aligning directly the inserted sequences called in both runs show many imperfect alignment so the divergence might be due to the reads/calling rather than the calling of the internal SNVs/indels.

More in the [report](hg002-allVsRegionalDownsampling.md).
