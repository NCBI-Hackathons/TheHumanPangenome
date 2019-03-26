library(methods)
library(GenomicRanges)
library(VariantAnnotation)
library(sveval)
library(dplyr)
library(parallel)

## Rscript callVariantsInInsertedSeq.R ORIG_VCF CALL_VCF OUTPUT_TSV

## Make sure both VCFs are split and normalized
## E.g. "bcftools norm -m -both -f REF.fa"

args = commandArgs(TRUE)
svpop.vcf = args[1] # e.g. 'SVPOP-explicit.vcf.gz'
call.vcf = args[2] # 'vg-call-HG002-illumina.vcf.gz'
output.file = args[3]
nb.cores = 3

## Optional parameter: BED file with regions of interest
regions.bed = NULL
if(length(args)==4){
  regions.bed = args[4]
}

## Read original VCF
svpop = readVcf(svpop.vcf)
svpop = rowRanges(svpop)
## Keep insertions defined as REF=1bp + ALT>1bp
svpop$size = unlist(Biostrings::nchar(svpop$ALT))
ref.size = Biostrings::nchar(svpop$REF)
svpop = svpop[which(ref.size==1 & svpop$size>1)]

## Read VCF with called variants
call = readSVvcf(call.vcf, keep.ins.seq=TRUE)
## Maybe some quality filtering here ??????????????
## Keep insertions
call = call[which(call$type=='INS')]

## Optional: Keep only variants around the subset insersions (+- 1kbp)
if(!is.null(regions.bed)){
  bed = read.table(regions.bed, as.is=TRUE, sep='\t')
  colnames(bed)[1:3] = c('chr','start','end')
  bed = makeGRangesFromDataFrame(bed)
  svpop = subsetByOverlaps(svpop, bed, maxgap=1e3)
  call = subsetByOverlaps(call, bed, maxgap=1e3)
}

## Match original insertions with called insertions by location and size
max.ins.gap = 30 # Maximum distance to cluster insertions
rsim.ins.size = .8 # Minimum reciprocal size similarity to match insertions
ol = findOverlaps(svpop, call, maxgap=max.ins.gap)
ol.df = as.data.frame(ol) %>%
  mutate(svpop.size=svpop$size[queryHits],
         call.size=call$size[subjectHits]) %>%
  filter(svpop.size / call.size > rsim.ins.size, call.size / svpop.size > rsim.ins.size)

## Align sequence with Smith-Waterman
svpop.seq = svpop$ALT[ol.df$queryHits]
calls.seq = call$ALT[ol.df$subjectHits]
vars.df = mclapply(1:nrow(ol.df), function(ii){
  ## Alignment
  svpop.seq = svpop.seq[[ii]][[1]]
  calls.seq = calls.seq[[ii]][[1]]
  pas = pairwiseAlignment(svpop.seq,
                          calls.seq,
                          type='global')
  ## List variants
  vars = tibble()
  ## SNV as mismatches
  al.mm = mismatchTable(pas)
  if(nrow(al.mm)>0){
    vars = rbind(vars, tibble(pos=al.mm$PatternStart, ref=al.mm$PatternSubstring, alt=al.mm$SubjectSubstring, type='SNV'))
  }
  ## Gaps for indels
  al.svpop = as.character(pattern(pas))
  al.call = as.character(subject(pas))
  al.df = tibble(svpop.pos=rep(NA, nchar(al.svpop)), call.pos=rep(NA, nchar(al.svpop)))
  gap.idx = gregexpr('-', al.svpop, fixed=TRUE)[[1]]
  if(gap.idx[1] != -1){
    al.df$svpop.pos[gap.idx] = '-'
  }
  nogap.idx = which(is.na(al.df$svpop.pos))
  al.df$svpop.pos[nogap.idx] = 1:length(nogap.idx) + pas@pattern@range@start - 1
  gap.idx = gregexpr('-', al.call, fixed=TRUE)[[1]]
  if(gap.idx[1] != -1){
    al.df$call.pos[gap.idx] = '-'
  }
  nogap.idx = which(is.na(al.df$call.pos))
  al.df$call.pos[nogap.idx] = 1:length(nogap.idx) + pas@subject@range@start - 1
  ## Insertions
  ins = rle(al.df$svpop.pos=='-')
  ins.df = tibble(x=cumsum(c(1,ins$lengths[-length(ins$lengths)])), width=ins$lengths, gap=ins$values)
  if(ins.df$gap[1]){
    ins.df = ins.df[-1]
  }
  if(ins.df$gap[nrow(ins.df)]){
    ins.df = ins.df[-nrow(ins.df)]
  }
  ins.df = subset(ins.df, gap)
  if(nrow(ins.df)){
    ins.pos = as.numeric(al.df$svpop.pos[ins.df$x-1])
    refs = sapply(ins.pos, function(pos) as.character(svpop.seq[pos]))
    alts = sapply(1:nrow(ins.df), function(ii){
      s = as.numeric(al.df$call.pos[ins.df$x[ii]-1])
      e = as.numeric(al.df$call.pos[ins.df$x[ii]-1]) + ins.df$width[ii]
      as.character(calls.seq[s:e])
    })
    vars = rbind(vars, tibble(pos=ins.pos, ref=refs, alt=alts, type='INS'))
  }
  ## Deletions
  del = rle(al.df$call.pos=='-')
  del.df = tibble(x=cumsum(c(1,del$lengths[-length(del$lengths)])), width=del$lengths, gap=del$values)
  if(del.df$gap[1]){
    del.df = del.df[-1]
  }
  if(del.df$gap[nrow(del.df)]){
    del.df = del.df[-nrow(del.df)]
  }
  del.df = subset(del.df, gap)
  if(nrow(del.df)){
    del.pos = as.numeric(al.df$svpop.pos[del.df$x-1])
    alts = sapply(del.pos, function(pos) as.character(svpop.seq[pos]))
    refs = sapply(1:nrow(del.df), function(ii){
      s = as.numeric(al.df$svpop.pos[del.df$x[ii]-1])
      e = as.numeric(al.df$svpop.pos[del.df$x[ii]-1]) + del.df$width[ii]
      as.character(svpop.seq[s:e])
    })
    vars = rbind(vars, tibble(pos=del.pos, ref=refs, alt=alts, type='DEL'))
  }
  if(nrow(vars)==0){
    return(NULL)
  }
  ## How well were the sequence aligned (potentially for post-filtering)
  vars$align.prop.match = nmatch(pas) / max(ol.df$svpop.size[ii], ol.df$call.size[ii])
  ## Annotate variants with original insertion information
  vars$ins.chr = as.character(seqnames(svpop))[ol.df$queryHits[ii]]
  vars$ins.pos = start(svpop[ol.df$queryHits[ii]])
  vars$ins.size = nchar(svpop.seq)
  vars$id = names(svpop)[ol.df$queryHits[ii]]
  return(vars)
}, mc.cores=nb.cores)
vars.df = do.call(rbind, vars.df)

## Output variants
write.table(vars.df, file=output.file, sep='\t', quote=FALSE, row.names=FALSE)
