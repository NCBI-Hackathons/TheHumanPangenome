# Analysis of the short-reads from HG002

## Reads from GiaB

LINKS

## Reducing the amount of reads to map

Later, we want to map to SVs across a larger set of samples so we want to try to reduce the amount of reads to map.
We'll first test the effect of doing this on HG002 reads.
In addition to mapping all the reads we can try to strategies to reduce the amount of input reads.

Here we start from the BAM file with reads mapped to the GRCh37 procided by GiaB.
On the other samples we'll have BAM files with reads mapped to GRCh38??

### Keeping reads that could be relevant

For example keep: 

1. Reads in the regions of interest.
1. Unmapped reads.

Advantage: Smaller set of reads selected.

Drawback: Maybe missing reads with low mapping quality that would map to the SV.

Something like:

```
samtools view -h -L svs.bed BAM > kept.sam 
samtools view -f 4 BAM >> kept.sam 
samtools sort -O BAM kept.sam > kept.bam
```

### Filtering reads that we know are not relevant

Filter out reads that map uniquely to regions far from regions of interest.

Advantage: Safer because we keep all reads that could potentially benefit from the graph.

Drawback: More reads, larger file/cost/computation time.

E.g. filter out by region, quality 60 and flags *paired* and *in proper pair*:

```
samtools view -h -b -U filtered.bam -L filterout.bed -q 60 -f 3 BAM ## Maybe need to redirect output > /dev/null
```
