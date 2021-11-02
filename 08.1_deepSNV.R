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
cvSymbols = dfGenes$V1#[1:10]
# remove duplicated symbols and whitespace
cvSymbols = gsub(' ', '', cvSymbols)
length(cvSymbols)
cvSymbols = unique(cvSymbols); length(cvSymbols)
dfGenes = AnnotationDbi::select(org.Hs.eg.db, keys = cvSymbols, 
                                columns = c('ENTREZID'),
                                keytype = 'SYMBOL')

dim(dfGenes)
dfGenes = na.omit(dfGenes); dim(dfGenes)
rownames(dfGenes) = dfGenes$ENTREZID
# get the genes into GRanges object
oGRgenes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                 filter=list('gene_id'=dfGenes$ENTREZID))
table(rownames(dfGenes) %in% names(oGRgenes))
## synchronize both objects w.r.t. gene ids
dfGenes = dfGenes[names(oGRgenes),]
dim(dfGenes)
identical(rownames(dfGenes), names(oGRgenes))

### this chunk may vary w.r.t. file paths if run on 
### rosalind instead of local machine
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
  r = deepSNV(test = csFiles[2], control = csFiles[1], 
              regions=oGRgenes[x], q=30, model='betabin')
})
setwd(gcswd)
names(loSnv) = dfGenes[names(oGRgenes),'SYMBOL']
## saved after running on rosalind
save(loSnv, file='results/loSnv.rds')

## load from here in case this was run on a different machine
## or continue analysis from here
load('results/loSnv.rds')

lResults = lapply(loSnv, myDeepSNVSummary, value='data.frame')
## create a data frame
dfResults = do.call(rbind, lResults)
i = match(rownames(dfResults), dfGenes$SYMBOL)
identical(rownames(dfResults), dfGenes$SYMBOL[i])
dfResults$ENTREZID = dfGenes$ENTREZID[i]
write.csv(dfResults, file='results/pilotDataVariantCalls.csv')

## convert to VCF
loVcf = lapply(loSnv, myDeepSNVSummary, value='VCF')
# drop the null values
loVcf[sapply(loVcf, is.null)] = NULL
oVcf = do.call(rbind, loVcf)
writeVcf(oVcf, file='results/pilotResults.vcf')
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
