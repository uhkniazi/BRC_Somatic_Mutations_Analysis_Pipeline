# File: 04_fastQC_trimmed.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the fastq files after trimming
# Date: 18/10/2021


## set variables and source libraries
source('header.R')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CFastqQuality/experimental/CFastqQuality.R'
download(url, 'CFastqQuality.R')

# load the required packages
source('CFastqQuality.R')
# delete the file after source
unlink('CFastqQuality.R')

# meta data
dfSample = read.csv('dataExternal/fileList.csv', header=T, stringsAsFactors = F)

#### get the names of the fastq files 
setwd('dataExternal/remote/Trimmed/')
csFiles = list.files('.', pattern = '*.gz', recursive = F)

cvMatch = gsub('_\\d\\.fastq.gz', '', as.character(csFiles))
cvMatch = gsub('trim_', '', cvMatch)
i = match(cvMatch, dfSample$SampleID)
dfFiles = data.frame(name=csFiles, id=dfSample$SampleID[i], group1=dfSample$Group[i],
                     stringsAsFactors = F)

# sanity check if all files present in directory
identical(dfFiles$name, csFiles)

## create factors and variables for fastqqualitybatch object
i = grep('_1.fastq.gz', dfFiles$name)
fReadDirection = rep('2', times=nrow(dfFiles))
fReadDirection[i] = '1'
fReadDirection = factor(fReadDirection)
# unique names
i = dfFiles$group1
table(i)
table(i, fReadDirection)

# unique names for each file
cNames = paste(i, dfFiles$id, as.character(fReadDirection),  sep='_')
lMetaData = list(files=dfFiles)

ob = CFastqQualityBatch(dfFiles$name, cNames, fReadDirection, lMetaData)

setwd(gcswd)
dir.create('results')
n = make.names(paste('CFastqQualityBatch george pilot data rds'))
n2 = paste0('results/', n)
save(ob, file=n2)

### create the plots of interest
getwd()

# set short names for samples for plotting etc
csShortNames = gsub('_WTCHG_894903', '', ob@csSampleNames)
iGetReadCount(ob)
barplot.readcount(ob)
plot.alphabetcycle(ob)
plot.qualitycycle(ob)
hist(iGetReadCount(ob), main='Read Count Distribution', xlab='')
######### some additional diagnostic plots on the data matrix
### some diagnostic plots
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

## extract the base quality matrix 
mBatch = mGetReadQualityByCycle(ob)
dim(mBatch)
mBatch[1:10, 1:4]
colnames(mBatch)
colnames(mBatch) = csShortNames
## creat an object of diagnostics class to make plots
oDiag = CDiagnosticPlots(mBatch, 'Base Quality')

## turning off automatic jitters
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

fBatch = ob@fReadDirection
str(ob@lMeta$files)
fBatch = factor(ob@lMeta$files$group1)
levels(fBatch)
## try some various factors to make the plots of low dimensional summaries
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
boxplot.median.summary(oDiag, fBatch)
plot.PCA(oDiag, fBatch, csLabels = ob@lMeta$files$idSample, cex=2)
plot.dendogram(oDiag, fBatch, labels_cex = 1)

## looking at alphabets 
## change direction and alphabet i.e. base as required
i = grep('1', ob@fReadDirection)

lAlphabets = lapply(i, function(x){
  m = t(mGetAlphabetByCycle(ob@lData[[x]]))
  m = m[,c('A', 'T', 'G', 'C')]
  r = rowSums(m)
  m = sweep(m, 1, r, '/')
  return(m)
})

mAlphabet = do.call(cbind, lapply(lAlphabets, function(x) return(x[,'C'])))
dim(mAlphabet)
i = grep('1', ob@fReadDirection)
colnames(mAlphabet) = csShortNames[i]
oDiag.2 = CDiagnosticPlots(mAlphabet, 'forward base C')

## turning off automatic jitters
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag.2)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)

i = grep('1', ob@fReadDirection)
fBatch = factor(ob@lMeta$files$group1[i])
length(fBatch)
levels(fBatch)
## try some various factors to make the plots of low dimensional summaries
plot.mean.summary(oDiag.2, fBatch)
plot.sigma.summary(oDiag.2, fBatch)
boxplot.median.summary(oDiag.2, fBatch)
plot.PCA(oDiag.2, fBatch, cex=2)
plot.PCA(oDiag.2, fBatch, xlim=c(-5, 5), cex=2, ylim=c(-1, 0.5))
plot.dendogram(oDiag.2, fBatch, labels_cex = 1)

############ make a figure for the read counts
mRC = as.matrix(iGetReadCount(ob))
dim(mRC)
rownames(mRC) = csShortNames
hc = hclust(dist(mRC))
plot(hc)
fBatch = factor(ob@lMeta$files$group1)
levels(fBatch)
# get a clustering variable
c = cutree(hc, k = 3)
table(c)
dfReads = data.frame(mRC, c, fBatch)
aggregate(dfReads$mRC, by=list(c), mean)
table(c)
aggregate(dfReads$mRC, by=list(fBatch, c), mean)

col.p = rainbow(nlevels(fBatch))
dend = as.dendrogram(hc)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) = col.p[as.numeric(fBatch)][order.dendrogram(dend)]
labels_cex(dend) = 1
# Plotting the new dendrogram
plot(dend, main=paste('Hierarchical clustering of distance matrix for read counts'),
     xlab='', sub='Coloured on Genotype')
abline(h = 2)
