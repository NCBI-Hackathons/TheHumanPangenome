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

## Optional: Keep only variants around the subset insersions (+- 100)
if(!is.null(regions.bed)){
  bed = read.table(regions.bed, as.is=TRUE, sep='\t')
  colnames(bed)[1:3] = c('chr','start','end')
  bed = makeGRangesFromDataFrame(bed)
  svpop = subsetByOverlaps(svpop, bed, maxgap=100)
  call = subsetByOverlaps(call, bed, maxgap=100)
}

## Remove insertions that are too large to align (by the pairwiseAlignment function)
## Max 45kbp (only a couple of insertions)
message(sum(svpop$size >= 45000), ' large insertions filtereds because too large to align.')
svpop = svpop[which(svpop$size < 45000)]

## Make sure ALT sequences are in a consistent format (oh the joys of R types...)
if(class(call$ALT)=='DNAStringSet'){
  call$ALT = split(call$ALT, 1:length(call$ALT))
}
if(class(svpop$ALT)=='DNAStringSet'){
  svpop$ALT = split(svpop$ALT, 1:length(svpop$ALT))
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
  ## message(ii)
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
  ## Looking at gaps to get indels
  ## Match positions in alignment with positions in input sequences
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
  ## Insertions: consecutive gaps in the query
  ## merge consecutive gap status, i.e. segment of consecutive gaps (or sequence)
  ins = rle(al.df$svpop.pos=='-')
  ins.df = tibble(x=cumsum(c(1,ins$lengths[-length(ins$lengths)])), width=ins$lengths, gap=ins$values)
  ## Remove gaps at the beginning/end of the alignment
  if(ins.df$gap[1]){
    ins.df = ins.df[-1]
  }
  if(ins.df$gap[nrow(ins.df)]){
    ins.df = ins.df[-nrow(ins.df)]
  }
  ## Keep only the gaps
  ins.df = subset(ins.df, gap)
  if(nrow(ins.df)){
    ## If any, get position and alleles
    ins.pos = as.numeric(al.df$svpop.pos[ins.df$x-1])
    refs = sapply(ins.pos, function(pos) as.character(svpop.seq[pos]))
    alts = sapply(1:nrow(ins.df), function(ii){
      s = as.numeric(al.df$call.pos[ins.df$x[ii]-1])
      e = as.numeric(al.df$call.pos[ins.df$x[ii]-1]) + ins.df$width[ii]
      as.character(calls.seq[s:e])
    })
    vars = rbind(vars, tibble(pos=ins.pos, ref=refs, alt=alts, type='INS'))
  }
  ## Deletions: gaps in the subject
  ## merge consecutive gap status, i.e. segment of consecutive gaps (or sequence)
  del = rle(al.df$call.pos=='-')
  del.df = tibble(x=cumsum(c(1,del$lengths[-length(del$lengths)])), width=del$lengths, gap=del$values)
  ## Remove gaps at the beginning/end of the alignment
  if(del.df$gap[1]){
    del.df = del.df[-1]
  }
  if(del.df$gap[nrow(del.df)]){
    del.df = del.df[-nrow(del.df)]
  }
  ## Keep only the gaps
  del.df = subset(del.df, gap)
  if(nrow(del.df)){
    ## If any, get position and alleles
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
  vars$align.prop.match = round(vars$align.prop.match, 4)
  ## Annotate variants with original insertion information
  vars$ins.chr = as.character(seqnames(svpop))[ol.df$queryHits[ii]]
  vars$ins.pos = start(svpop[ol.df$queryHits[ii]])
  vars$ins.size = nchar(svpop.seq)
  vars$id = names(svpop)[ol.df$queryHits[ii]]
  vars$ins.gt = call$GT[ol.df$subjectHits[ii]]
  return(vars)
}, mc.cores=nb.cores)

vars.df = do.call(rbind, vars.df)

## Output variants
write.table(vars.df, file=output.file, sep='\t', quote=FALSE, row.names=FALSE)

## Columns in output
## pos: position of the variant relative to the inserted sequence (in the original VCF)
## ref: reference allele (like in VCF)
## alt: alternate allele (like in VCF)
## type: type of variant (SNV, DEL, INS)
## id: id of the insertion, e.g. NA19240_chr1-136935-INS-244
## ins.chr/ins.pos: location of the insertion in the genome
## ins.gt: genotype of the called insertion.
## align.prop.match: proportion of matches when aligning inserted sequence (for filtering maybe)
