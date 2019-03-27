library(methods)
library(GenomicRanges)
library(VariantAnnotation)
library(sveval)
library(dplyr)
library(parallel)

args = commandArgs(TRUE)
call.vcf = args[1] 
call2.vcf = args[2] 
output.file = args[3]
nb.cores = 10

## Optional parameter: BED file with regions of interest
regions.bed = NULL
if(length(args)==4){
  regions.bed = args[4]
}

## Read VCF with called variants
message('Read calls VCF 1')
call = readSVvcf(call.vcf, keep.ins.seq=TRUE)
## Maybe some quality filtering here ??????????????
## Keep insertions
call = call[which(call$type=='INS')]

## Read VCF with called variants
message('Read calls VCF 2')
call2 = readSVvcf(call2.vcf, keep.ins.seq=TRUE)
## Maybe some quality filtering here ??????????????
## Keep insertions
call2 = call2[which(call2$type=='INS')]

## Optional: Keep only variants around the subset insersions (+- 100)
if(!is.null(regions.bed)){
  bed = read.table(regions.bed, as.is=TRUE, sep='\t')
  colnames(bed)[1:3] = c('chr','start','end')
  bed = makeGRangesFromDataFrame(bed)
  call2 = subsetByOverlaps(call2, bed, maxgap=100)
  call = subsetByOverlaps(call, bed, maxgap=100)
}

## Remove insertions that are too large to align (by the pairwiseAlignment function)
## Max 45kbp (only a couple of insertions)
message(sum(call2$size >= 45000), ' large insertions filtered because too large to align.')
call2 = call2[which(call2$size < 45000)]

## Make sure ALT sequences are in a consistent format (oh the joys of R types...)
if(class(call$ALT)=='DNAStringSet'){
  call$ALT = split(call$ALT, 1:length(call$ALT))
}
if(class(call2$ALT)=='DNAStringSet'){
  call2$ALT = split(call2$ALT, 1:length(call2$ALT))
}

## Match original insertions with called insertions by location and size
message('Overlap')
max.ins.gap = 30 # Maximum distance to cluster insertions
rsim.ins.size = .8 # Minimum reciprocal size similarity to match insertions
ol = findOverlaps(call2, call, maxgap=max.ins.gap)
ol.df = as.data.frame(ol) %>%
  mutate(call2.size=call2$size[queryHits],
         call.size=call$size[subjectHits]) %>%
  filter(call2.size / call.size > rsim.ins.size, call.size / call2.size > rsim.ins.size)

## Align sequence with Smith-Waterman
message("Align")
call2.seq = call2$ALT[ol.df$queryHits]
calls.seq = call$ALT[ol.df$subjectHits]
al.df = mclapply(1:nrow(ol.df), function(ii){
  ## message(ii)
  ## Alignment
  call2.seq = call2.seq[[ii]][[1]]
  calls.seq = calls.seq[[ii]][[1]]
  pas = pairwiseAlignment(call2.seq,
                          calls.seq,
                          type='global')
  ## How well were the sequence aligned (potentially for post-filtering)
  return(tibble(align.prop.match = nmatch(pas) / max(ol.df$call2.size[ii], ol.df$call.size[ii]),
                id = names(call2)[ol.df$queryHits[ii]]))
}, mc.cores=nb.cores)

al.df = do.call(rbind, al.df)

## Output variants
write.table(al.df, file=output.file, sep='\t', quote=FALSE, row.names=FALSE)

