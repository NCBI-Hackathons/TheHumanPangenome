var1 = read.table('variantsInInsertions.1.tsv', as.is=TRUE, header=TRUE)
var2 = read.table('variantsInInsertions.2.tsv', as.is=TRUE, header=TRUE)

## Create ids and overlap for Venn diagram
var1$var.id = paste(var1$pos, var1$ref, var1$alt, var1$id, var1$ins.gt, sep='_')
var2$var.id = paste(var2$pos, var2$ref, var2$alt, var2$id, var2$ins.gt, sep='_')

length(intersect(var1$var.id, var2$var.id))
length(setdiff(var1$var.id, var2$var.id))
length(setdiff(var2$var.id, var1$var.id))

library(VennDiagram)

venn.diagram(list(var1=var1$var.id, var2=var2$var.id), filename='test.tiff',
             col='transparent', fill=c('steelblue', 'indianred2'))


## Merge
var12 = merge(var1, var2, by=c('pos', 'ref', 'alt', 'id', 'ins.gt', 'type', 'ins.chr', 'ins.pos', 'ins.size'), all=TRUE,
              suffixes=c(1,2))
var12$in1 = !is.na(var12$align.prop.match1)
var12$in2 = !is.na(var12$align.prop.match2)
var12$n = as.numeric(var12$in1) + as.numeric(var12$in2)
var12$size = ifelse(var12$type=='INS', nchar(var12$alt) - 1, nchar(var12$ref) - 1)

library(ggplot2)
ggplot(var12, aes(x=type, fill=factor(n))) + geom_bar(position='fill') + theme_bw()

ggplot(subset(var12, type!='SNV'), aes(x=size, fill=factor(n))) + geom_histogram() + theme_bw()


## Examples
## Unique to 1
head(subset(var12, n==1 & in1))
## Unique to 2
head(subset(var12, n==1 & in2))
