# File: 08.1_deepSNV.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: somatic mutations analysis
#       https://bioconductor.org/packages/release/bioc/vignettes/deepSNV/inst/doc/deepSNV.pdf
# Date: 28/10/2021

source('header.R')
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# select some genes of interest
dfGenes = read.csv('dataExternal/pilotProjectGeneList.csv', stringsAsFactors = F,
                   header=F)
# select some genes of interest
cvSymbols = dfGenes$V1[1:10]
dfGenes = AnnotationDbi::select(org.Hs.eg.db, keys = cvSymbols, 
                                columns = c('ENTREZID'),
                                keytype = 'SYMBOL')

dfGenes
rownames(dfGenes) = dfGenes$ENTREZID
# get the genes into GRanges object
oGRgenes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
oGRgenes = oGRgenes[as.character(dfGenes$ENTREZID)]

# load meta data
dfMeta = read.csv('dataExternal/fileList.csv', header=T, stringsAsFactors = F)
str(dfMeta)
table(dfMeta$Group)
# match group names with available files
csFiles = list.files('dataExternal/subsample/', '*.bam$')
f = gsub('_q30_.+', '', csFiles)
dfMeta = dfMeta[dfMeta$SampleID %in% f,]
i = match(dfMeta$SampleID, f)
# sanity check
identical(dfMeta$SampleID, f[i])
dfMeta$files = csFiles[i]
setwd('dataExternal/subsample/')

library(deepSNV)
loSnv = lapply(seq_along(oGRgenes), function(x){
  r = deepSNV(test = dfMeta$files[1], control = dfMeta$files[2], 
              regions=oGRgenes[x], q=30, model='betabin')
})
setwd(gcswd)
save(loSnv, file='results/loSnv.rds')
names(loSnv) = dfGenes[names(oGRgenes),'SYMBOL']

lResults = lapply(loSnv, summary, value='data.frame')
## create a data frame
dfResults = do.call(rbind, lResults)

loVcf = lapply(loSnv, myDeepSNVSummary, value='VCF')
# drop the null values
loVcf[sapply(loVcf, is.null)] = NULL
oVcf = do.call(rbind, loVcf)

######################### 
##### test code
#########################
# oSnv = deepSNV(test = dfMeta$files[1], 
#                control = dfMeta$files[2], regions=oGRgenes[1], q=30, model='betabin')
# summary(oSnv)
# plot(oSnv)
# show(oSnv)
# c = control(oSnv)
# t = test(oSnv)
# dim(c)
# rownames(c) = start(oGRgenes):end(oGRgenes)
# rownames(t) = start(oGRgenes):end(oGRgenes)
# s = summary(oSnv)
# c[as.character(s$pos),]
# t[as.character(s$pos),]
# temp = summary(oSnv, value='data.frame')
# # use the modified version of the deepsnv summary function 
# # this is done because the original function has an error when producing
# # VCF output. 
# temp2 = myDeepSNVSummary(oSnv, value='data.frame')
# # sanity check
# identical(temp, temp2)
# oVcf = myDeepSNVSummary(oSnv, value='VCF')
# setwd(gcswd)
# writeVcf(oVcf, file='results/vcf_results.vcf')
#########################
