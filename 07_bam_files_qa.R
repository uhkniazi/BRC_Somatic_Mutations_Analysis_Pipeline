# File: 07_bam_files_qa.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the bam files
# Date: 17/11/2021


## set variables and source libraries
source('header.R')
library(Rsamtools)

## load the metadata
dfMeta = read.csv('dataExternal/fileList.csv', header=T, stringsAsFactors = F)

#### set working directory to appropriate location with bam files
setwd('dataExternal/remote/Aligned/')
csFiles = list.files('.', pattern = '*sortc.bam$')
## match the file names to sample ids in metadata table
f = gsub('_q30_.+', '', csFiles)
table(dfMeta$SampleID %in% f)
i = match(dfMeta$SampleID, f)
identical(dfMeta$SampleID, f[i])
dfMeta$files = csFiles[i]

## load the list of files with duplicates removed
csFiles = list.files('.', pattern = '*sortc_rd.bam$')
f = gsub('_q30_.+', '', csFiles)
table(dfMeta$SampleID %in% f)
i = match(dfMeta$SampleID, f)
identical(dfMeta$SampleID, f[i])
dfMeta$files_rd = csFiles[i]

dim(dfMeta)
csFiles = dfMeta$files

## perform the analysis one sample at a time
lAllBams = lapply(csFiles, function(x){
  c = countBam(x)
  cat(paste(x, 'done', '\n'))
  return(c)
})
dfBam = do.call(rbind, lAllBams)

## samples with duplicates removed
csFiles = dfMeta$files_rd
lAllBams = lapply(csFiles, function(x){
  c = countBam(x)
  cat(paste(x, 'done', '\n'))
  return(c)
})
df = do.call(rbind, lAllBams)
dfBam = rbind(dfBam, df)

setwd(gcswd)
n2 = paste0('results/', 'bamFileCounts.csv')
write.csv(dfBam, file=n2)

### load the results if this part was run on cluster
dfBam = read.csv(n2, header=T, stringsAsFactors = F, row.names=1)

i = grep('sortc.bam', dfBam$file)
iReadCount.pre = round(dfBam$records[i] / 1e6, 2)
names(iReadCount.pre) = gsub('_q30_.+', '', dfBam$file[i])
iReadCount.post = round(dfBam$records[-i] / 1e6, 2)
names(iReadCount.post) = gsub('_q30_.+', '', dfBam$file[-i])
identical(names(iReadCount.pre), names(iReadCount.post))
## check if or otherwise sort names in same order as meta data
identical(dfMeta$SampleID, names(iReadCount.pre))
identical(dfMeta$SampleID, names(iReadCount.post))
names(iReadCount.pre) = dfMeta$Group
names(iReadCount.post) = dfMeta$Group

### create the plots of interest
setwd(gcswd)
pdf(file='results/bam.qa.pdf')

mReadCount = rbind(iReadCount.post, iReadCount.pre)
rownames(mReadCount) = c('Post', 'Pre')

barplot(mReadCount, beside=T, las=2, main='No. of reads aligned', ylab = 'No. of Reads in Millions', cex.names =1,
        ylim=c(0, max(mReadCount)+2))
legend('topright', legend = c('Dup Rem', 'Orig'), fill=c('black', 'grey'))

dev.off(dev.cur())
