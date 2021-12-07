# File: 08_shearwaterML.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: somatic mutations analysis
# Date: 27/09/2021

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
csFiles = list.files('dataExternal/remote/Aligned/', '*rd.bam$')
f = gsub('_q30_.+', '', csFiles)
dfMeta = dfMeta[dfMeta$SampleID %in% f,]
i = match(dfMeta$SampleID, f)
# sanity check
identical(dfMeta$SampleID, f[i])
dfMeta$files = csFiles[i]
setwd('dataExternal/remote/Aligned/')

## start with the 1st test sample and others as controls
csFiles = c(test=dfMeta$files[8], control=dfMeta$files[-8])
library(deepSNV)
loVCF = lapply(seq_along(oGRgenes), function(x){
  arCounts = loadAllData(csFiles, oGRgenes[x], q=30)
})
setwd(gcswd)
names(loSnv) = dfGenes[names(oGRgenes),'SYMBOL']
csSave = paste0('results/', csFiles['test'], '_loSnv.rds')
## saved after running on rosalind
save(loSnv, file=csSave)

## load from here in case this was run on a different machine
## or continue analysis from here
load('results/WTCHG_894903_73075379_q30_fixm_sortc_rd.bam_loSnv.rds')
iMetaIndex = 7

lResults = lapply(loSnv, myDeepSNVSummary, value='data.frame')
## create a data frame
dfResults = do.call(rbind, lResults)
i = match(rownames(dfResults), dfGenes$SYMBOL)
identical(rownames(dfResults), dfGenes$SYMBOL[i])
## if this is false, it is probably due to multiple hits on single gene
## otherwise skip this matching section
s = rownames(dfResults)
s = gsub('\\.\\d*', '', s)
i = match(s, dfGenes$SYMBOL)
identical(s, dfGenes$SYMBOL[i])

dfResults$ENTREZID = dfGenes$ENTREZID[i]
cvName = paste0('results/varCall', dfMeta$SampleID[iMetaIndex], '.csv')
write.csv(dfResults, file=cvName)

## convert to VCF
loVcf = lapply(loSnv, myDeepSNVSummary, value='VCF')
# drop the null values
loVcf[sapply(loVcf, is.null)] = NULL
oVcf = do.call(rbind, loVcf)
cvName = paste0('results/varCall', dfMeta$SampleID[iMetaIndex], '.vcf')
writeVcf(oVcf, file=cvName)


######### old test code
##################################################################
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
oGRgenes = oGRgenes[as.character(dfGenes$ENTREZID)]

library(deepSNV)
# see here for usage
# https://www.bioconductor.org/packages/devel/bioc/vignettes/deepSNV/inst/doc/shearwaterML.html
csFiles = list.files('dataExternal/subsample/', pattern = '*.bam$',full.names = T)

arCounts = loadAllData(csFiles, oGRgenes, q=10)
pvals = betabinLRT(arCounts)$pvals
qvals = p.adjust(pvals, method = 'BH')
dim(qvals) = dim(pvals)
oVcf = qvals2Vcf(qvals, arCounts, oGRgenes, samples = csFiles, mvcf = TRUE)
dir.create('results')
save(oVcf, file='results/oVcf.rds')
#################################################################