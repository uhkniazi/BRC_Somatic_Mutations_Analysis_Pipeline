# File: 09_annotateVCF.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Annotate the identified somatic mutations
# Date: 2/11/2021

source('header.R')

library(VariantAnnotation)

dfDeepSNV = read.csv('results/pilotDataVariantCalls.csv', header=T, row.names=1)
oVcf = readVcf('results/pilotResults.vcf', genome = 'hg38')
header(oVcf)
## header slots
samples(header(oVcf))
geno(header(oVcf))

## genomic positions
rowRanges(oVcf)
rowRanges(oVcf)$REF
ref(oVcf)
alt(oVcf)

## genotype data
# unique to each sample and may be unique for each variant and sample
geno(header(oVcf))
geno(oVcf)
sapply(geno(oVcf), class)
sapply(geno(oVcf), dim)
geno(oVcf)[[1]]

## info data
# unique to each variant, same across samples
info(header(oVcf))
info(oVcf)

## use the snplocs database
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
oDBSnp = SNPlocs.Hsapiens.dbSNP151.GRCh38

oGRsnps = rowRanges(oVcf)
#oGRsnps = as(oVcf, 'VRanges')
seqnames(oGRsnps)
seqlevels(oGRsnps)
## change the style to match snplocs package
seqlevelsStyle(oGRsnps); seqlevelsStyle(oDBSnp)
seqlevelsStyle(oGRsnps) = 'NCBI'
## check
seqnames(oGRsnps)
seqlevels(oGRsnps)
## do some acrobatics with seqlevels to perform search
## as this will give error
oGPsnps = snpsByOverlaps(oDBSnp, oGRsnps)

seqlevels(oGRsnps) = seqlevels(oDBSnp)
seqinfo(oGRsnps) = seqinfo(oDBSnp) 
oGPsnps = snpsByOverlaps(oDBSnp, oGRsnps)

# manually check a few rs numbers if any hits at
# https://www.ncbi.nlm.nih.gov/snp/?term=rs3745946
#IUPAC_CODE_MAP
oGPsnps

f = findOverlaps(oGRsnps, oGPsnps)
oGRsnps$RefSNP_id = NA_character_
oGRsnps$RefSNP_id[queryHits(f)] = oGPsnps$RefSNP_id[subjectHits(f)]
oGRsnps.snploc = oGRsnps
## match the names with deepSNV results table
gr = GRanges(as.character(dfDeepSNV$chr), IRanges(start=dfDeepSNV$pos, width = 1))
seqlevelsStyle(oGRsnps.snploc) = 'UCSC'
i = findOverlaps(oGRsnps.snploc, gr)
identical(queryHits(i), subjectHits(i))
dfDeepSNV$RefSNP_id = oGRsnps.snploc$RefSNP_id

## location of the variants 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
oTxdb = TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevelsStyle(oTxdb)
oGRsnps = rowRanges(oVcf)
seqlevelsStyle(oGRsnps)
seqlevelsStyle(oTxdb) = 'UCSC'
genome(oGRsnps)
genome(oTxdb)

oGRlocations = locateVariants(oGRsnps, oTxdb, AllVariants())
library(org.Hs.eg.db)
## summarize each snp by location
dfSymbol = AnnotationDbi::select(org.Hs.eg.db, keys = oGRlocations$GENEID, columns = 'SYMBOL',
                                keytype = 'ENTREZID')
## check if names match with deepsnv results
table(dfDeepSNV$ENTREZID %in% oGRlocations$GENEID)
table(rownames(dfDeepSNV) %in% dfSymbol$SYMBOL)

#lSnpLoc = split(oGRlocations$LOCATION, oGRlocations$QUERYID)
lSnpLoc = split(oGRlocations$LOCATION, dfSymbol$SYMBOL)
length(lSnpLoc)
dfLocations = data.frame(sapply(lSnpLoc, function(x) table(x)))

write.csv(dfDeepSNV, file='results/pilotDataVariantCalls_RefSNP.csv')
write.csv(dfLocations, file='results/pilotDataVariantLocations.csv')
