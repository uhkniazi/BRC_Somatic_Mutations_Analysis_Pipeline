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

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'Sample')
# another way to get the query, preferred
g_did
dfSample = dbGetQuery(db, "select * from Sample where idData=46;")
# remove any whitespace from the names
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
head(dfSample)
## get the file names from the files table
dbListFields(db, 'File')

## write query to get file names
q = paste0('select File.*, Sample.id as sid, Sample.* from File, Sample where Sample.idData=46 AND Sample.id=File.idSample AND File.type="fastq"', ';')

# df = lapply(q, function(x) dbGetQuery(db, x))
# dfFiles = do.call(rbind, df)
dfFiles = dbGetQuery(db, q)
# close connection after getting data
dbDisconnect(db)
str(dfFiles)
# remove redundant columns
dfFiles = dfFiles[,-c(5,7)]
#### get the names of the fastq files for first sequencing run
setwd('dataExternal/remote/Trimmed/')
csFiles = list.files('.', pattern = '*.gz', recursive = F)

## drop the NA samples
dfFiles = na.omit(dfFiles)
dim(dfFiles)

# sanity check if all files present in directory
dfFiles$name = paste('trim_', dfFiles$name, sep ='')
table(dfFiles$name %in% csFiles)
identical(dfFiles$name, csFiles)
table(csFiles %in% dfFiles$name)
## create factors and variables for fastqqualitybatch object
i = grep('_1.fq.gz', dfFiles$name)
fReadDirection = rep('2', times=nrow(dfFiles))
fReadDirection[i] = '1'
fReadDirection = factor(fReadDirection)
i = dfFiles$group1
table(i)
i.2 = dfFiles$group2
table(i.2)
# unique names for each file
cNames = paste(i, i.2, dfFiles$idSample, as.character(fReadDirection),  sep='_')
lMetaData = list(files=dfFiles)

ob = CFastqQualityBatch(dfFiles$name, cNames, fReadDirection, lMetaData)

setwd(gcswd)
n = make.names(paste('CFastqQualityBatch post trimming data id 46 john rds'))
n2 = paste0('~/Data/MetaData/', n)
save(ob, file=n2)

## note: comment out as this entry has been made in db
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='trim FASTQ quality checks on John FOG1 KO mouse data')
#dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

### create the plots of interest
getwd()

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

## creat an object of diagnostics class to make plots
oDiag = CDiagnosticPlots(mBatch, 'Base Quality')

## turning off automatic jitters
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag = CDiagnosticPlotsSetParameters(oDiag, l)

fBatch = ob@fReadDirection
str(ob@lMeta$files)
d = ob@lMeta$files$description
d = strsplit(d, ';')
fBatch = factor(sapply(d, function(x) x[2]))
d = strsplit(ob@lMeta$files$group3, '_')
fBatch = factor(ob@lMeta$files$group1)
fBatch = fBatch:factor(ob@lMeta$files$group2)
f1 = factor(ob@lMeta$files$group2)
f2 = factor(ob@lMeta$files$group3)
f3 = factor(ob@lMeta$files$group1)
fBatch = f1:f2
levels(fBatch)
## try some various factors to make the plots of low dimensional summaries
plot.mean.summary(oDiag, fBatch)
plot.sigma.summary(oDiag, fBatch)
boxplot.median.summary(oDiag, fBatch)
plot.PCA(oDiag, fBatch, csLabels = ob@lMeta$files$idSample, cex=2)
plot.dendogram(oDiag, fBatch, labels_cex = 0.4)

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
colnames(mAlphabet) = ob@lMeta$files$idSample[i]
oDiag.2 = CDiagnosticPlots(mAlphabet, 'forward base C')

## turning off automatic jitters
## we set jitter to FALSE for PCA, otherwise, in sparse matrix a random jitter is added to avoid divisions by zeros
l = CDiagnosticPlotsGetParameters(oDiag.2)
l$PCA.jitter = F; l$HC.jitter = F;
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)

i = grep('1', ob@fReadDirection)
fBatch = factor(ob@lMeta$files$group2[i])
length(fBatch)
f1 = factor(ob@lMeta$files$group2[i])
f2 = factor(ob@lMeta$files$group3[i])
f3 = factor(ob@lMeta$files$group1[i])
fBatch = f1:f2#:f3
levels(fBatch)
## try some various factors to make the plots of low dimensional summaries
plot.mean.summary(oDiag.2, fBatch)
plot.sigma.summary(oDiag.2, fBatch)
boxplot.median.summary(oDiag.2, fBatch)
plot.PCA(oDiag.2, fBatch, cex=2)
plot.PCA(oDiag.2, fBatch, xlim=c(-10, 10), cex=2)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8)

############ make a figure for the read counts
mRC = as.matrix(iGetReadCount(ob))
rownames(mRC) = ob@lMeta$files$id
hc = hclust(dist(mRC))
# lane
d = ob@lMeta$files$description
d = strsplit(d, ';')
fLane = factor(sapply(d, function(x) x[2]))
f1 = factor(ob@lMeta$files$group2)
f2 = factor(ob@lMeta$files$group3)
f3 = factor(ob@lMeta$files$group1)
fBatch = f1:f2
levels(fBatch)
# get a clustering variable
c = cutree(hc, h=2)
table(c)
dfReads = data.frame(mRC, c, fBatch)
aggregate(dfReads$mRC, by=list(c), mean)
table(c)

col.p = rainbow(nlevels(fBatch))
dend = as.dendrogram(hc)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) = col.p[as.numeric(fBatch)][order.dendrogram(dend)]
labels_cex(dend) = 0.5
# Plotting the new dendrogram
plot(dend, main=paste('Hierarchical clustering of distance matrix for read counts'),
     xlab='', sub='Coloured on Genotype:Replicate')
abline(h = 2)