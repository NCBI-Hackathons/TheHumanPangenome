library(methods)
library(GenomicRanges)
library(VariantAnnotation)
library(sveval)
library(dplyr)

## Rscript callVariantsInInsertedSeq.R ORIG_VCF CALL_VCF OUTPUT_TSV

## Make sure both VCFs are split and normalized
## E.g. "bcftools norm -m -both -f REF.fa"

arg = commandArgs(TRUE)
svpop.vcf = args[1] # e.g. 'SVPOP-explicit.vcf.gz'
call.vcf = args[2] # 'vg-call-HG002-illumina.vcf.gz'
output.file = args[3]
nb.cores = 3


## Read original VCF
svpop = readVcf(svpop.vcf)
svpop = rowRanges(svpop)
## Keep insertions defined as REF=1bp + ALT>1bp
svpop$size = Biostrings::nchar(svpop$ALT)
ref.size = Biostrings::nchar(svpop$REF)
svpop = svpop[which(ref.size==1 & svpop$size>1)]

## Read VCF with called variants
call = readSVvcf(call.vcf, keep.ins.seq=TRUE)
## Maybe some quality filtering here ??????????????
## Keep insertions defined as REF=1bp + ALT>1bp
call$size = Biostrings::nchar(call$ALT)
ref.size = Biostrings::nchar(call$REF)
call = call[which(ref.size==1 & call$size>1)]

## Match original insertions with called insertions by location and size
max.ins.gap = 30 # Maximum distance to cluster insertions
rsim.ins.size = .8 # Minimum reciprocal size similarity to match insertions
ol = findOverlaps(svpop, call, maxgap=max.ins.gap)
ol.df = as.data.frame(ol) %>%
  mutate(svpop.size=svpop$size[queryHits],
         call.size=call$size[subjectHits],
         ol.size=width(pintersect(svpop[queryHits], call[subjectHits]))) %>%
  filter(ol.size / call.size > rsim.ins.size, ol.size / svpop.size > rsim.ins.size)

## Align sequence with Smith-Waterman
svpop.seq = svpop$ALT[ol.df$queryHits]
call.seq = call$ALT[ol.df$subjectHits]
vars.df = parallel::mclapply(1:nrow(ol.df), function(ii){
  ## Alignment
  svpop.seq = DNAString(svpop.seq[ol.df$queryHits[ii]])
  call.seq = DNAString(calls.seq[ol.df$subjectHits[ii]])
  pas = pairwiseAlignment(svpop.seq,
                          call.seq,
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
  al.df$svpop.pos[gregexpr('-', al.svpop, fixed=TRUE)[[1]]] = '-'
  al.df$svpop.pos[which(is.na(al.df$svpop.pos))] = 1:nchar(svpop.seq)
  al.df$call.pos[gregexpr('-', al.call, fixed=TRUE)[[1]]] = '-'
  al.df$call.pos[which(is.na(al.df$call.pos))] = 1:nchar(call.seq)
  ## Insertions
  ins = rle(al.df$svpop.pos=='-')
  ins.df = tibble(x=cumsum(c(1,ins$lengths[-length(ins$lengths)])), width=ins$lengths, gap=ins$values) %>% filter(gap)
  ins.pos = as.numeric(al.df$svpop.pos[ins.df$x-1])
  refs = sapply(ins.pos, function(pos) as.character(svpop.seq[pos]))
  alts = sapply(1:nrow(ins.df), function(ii){
    s = as.numeric(al.df$call.pos[ins.df$x[ii]-1])
    e = as.numeric(al.df$call.pos[ins.df$x[ii]-1]) + ins.df$width[ii]
    as.character(call.seq[s:e])
  })
  vars = rbind(vars, tibble(pos=ins.pos, ref=refs, alt=alts, type='INS'))
  ## Deletions
  del = rle(al.df$call.pos=='-')
  del.df = tibble(x=cumsum(c(1,del$lengths[-length(del$lengths)])), width=del$lengths, gap=del$values) %>% filter(gap)
  del.pos = as.numeric(al.df$svpop.pos[del.df$x-1])
  alts = sapply(del.pos, function(pos) as.character(svpop.seq[pos]))
  refs = sapply(1:nrow(del.df), function(ii){
    s = as.numeric(al.df$svpop.pos[del.df$x[ii]-1])
    e = as.numeric(al.df$svpop.pos[del.df$x[ii]-1]) + del.df$width[ii]
    as.character(svpop.seq[s:e])
  })
  vars = rbind(vars, tibble(pos=del.pos, ref=refs, alt=alts, type='DEL'))
  ## How well were the sequence aligned (potentially for post-filtering)
  vars$prop.match = nmatch(pas) / max(ol.df$svpop.size[ii], ol.df$call.size[ii])
  ## Annotate variants with original insertion information
  vars$ins.pos = start(svpop[ol.df$queryHits])
  vars$ins.size = nchar(svpop.seq)
  return(vars)
}, mc.cores=nb.cores)
vars.df = do.call(rbind, vars.df)


## Output variants
write.table(vars.df, file=output.file, sep='\t', quote=FALSE, row.names=FALSE)
