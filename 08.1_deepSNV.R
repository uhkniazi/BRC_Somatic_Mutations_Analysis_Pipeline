# File: 08.1_deepSNV.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: somatic mutations analysis
#       https://bioconductor.org/packages/release/bioc/vignettes/deepSNV/inst/doc/deepSNV.pdf
# Date: 28/09/2021

source('header.R')
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# select some genes of interest
cvSymbols = c('NOTCH1', 'TP53', 'TP63')
dfGenes = AnnotationDbi::select(org.Hs.eg.db, keys = cvSymbols, 
                                columns = c('ENTREZID'),
                                keytype = 'SYMBOL')

dfGenes

# get the genes into GRanges object
oGRgenes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
oGRgenes = oGRgenes[as.character(dfGenes$ENTREZID[1])]

dfMeta = read.csv(file.choose(), header=T)
dfMeta[, c('Run', 'tissue')]

library(deepSNV)
csFiles = list.files('dataExternal/subsample/', pattern = '*.bam$',full.names = T)
csFiles

oSnv = deepSNV(test = csFiles[1], control = csFiles[2], regions=oGRgenes, q=30, model='betabin')
summary(oSnv)
plot(oSnv)
show(oSnv)
c = control(oSnv)
t = test(oSnv)
dim(c)
rownames(c) = start(oGRgenes):end(oGRgenes)
rownames(t) = start(oGRgenes):end(oGRgenes)
s = summary(oSnv)
c[as.character(s$pos),]
t[as.character(s$pos),]
temp = summary(oSnv, value='data.frame')
temp2 = myDeepSNVSummary(oSnv, value='data.frame')
oVcf = myDeepSNVSummary(oSnv, value='VCF')
writeVcf(oVcf, file='results/vcf_results.vcf')

