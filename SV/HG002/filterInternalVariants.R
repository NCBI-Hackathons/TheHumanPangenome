library(methods)

## Rscript filterInternalVariants.R INPUT_TSV OUTPUT_TSV

args = commandArgs(TRUE)
vars.in = args[1]
vars.out = args[2]

vars = read.table(vars.in, as.is=TRUE, header=TRUE)

vars = subset(vars, align.prop.match>.7)

write.table(vars, file=vars.out, sep='\t', row.names=FALSE, quote=FALSE)
