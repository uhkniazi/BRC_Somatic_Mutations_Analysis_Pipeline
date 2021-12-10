# File: 07_bamSubsampleForGenes.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: filter a set of bam files for regions of interest
# Date: 27/09/2021

source('header.R')

require(Rsamtools)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(org.Hs.eg.db)

dfGenes = read.csv('dataExternal/pilotProjectGeneList.csv', stringsAsFactors = F,
                   header=F)
# select some genes of interest
cvSymbols = dfGenes$V1[1:10]
dfGenes = AnnotationDbi::select(org.Hs.eg.db, keys = cvSymbols, 
                                columns = c('ENTREZID'),
                                keytype = 'SYMBOL')

dfGenes

# get the genes into GRanges object
oGRgenes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
oGRgenes = oGRgenes[as.character(dfGenes$ENTREZID)]

# expand the query region size
oGRquery = resize(oGRgenes, width = width(oGRgenes)+10000, fix='center')
strand(oGRquery) = '*'

#### set working directory to appropriate location with bam files
setwd('dataExternal/remote/')
setwd('Aligned/')
csFiles = list.files('.', pattern = '*rd.bam$', recursive = T)
#csFiles = csFiles[c(1, 2)]

dir.create('subsample')

# ### testing reading one file
# param = ScanBamParam(which = oGRquery)
# d = filterBam(dfSample$name[2], destination = 'output/testing_2.bam', param=param)
# dir.create('merged')
# d2 = mergeBam(c('output/testing.bam', 'output/testing_2.bam'),
#               'merged/merged.bam', indexDestination=T, overwrite=T)
## filter the bam files 
param = ScanBamParam(which = oGRquery)
d = sapply(csFiles, function(x){
  filterBam(x, destination = paste0('subsample/', x, '.bam'), param=param)
  indexBam(paste0('subsample/', x, '.bam'))
})

list.files('subsample')
setwd(gcswd)
