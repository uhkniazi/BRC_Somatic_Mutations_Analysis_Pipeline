# File: 09.1_annotateVariantTable.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Annotate the identified somatic mutations using shearwaterML.R script
# Date: 14/12/2021

source('header.R')

## load the relevant files
dfDeepSNV = read.csv('temp/dfResults.adj.csv', header=T, row.names=1, stringsAsFactors = F)

## use the snplocs database
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
oDBSnp = SNPlocs.Hsapiens.dbSNP151.GRCh38
df = DataFrame(dfDeepSNV[,-c(2,3)]); rownames(df) = NULL
dim(df)
oGRsnps = GRanges(dfDeepSNV$chr, IRanges(start=dfDeepSNV$pos, width=1), 
                  strand = '*')
mcols(oGRsnps) = df
length(oGRsnps)

seqnames(oGRsnps)
seqlevels(oGRsnps)
## change the style to match snplocs package
seqlevelsStyle(oGRsnps); seqlevelsStyle(oDBSnp)
seqlevelsStyle(oGRsnps) = 'NCBI'
## check
seqnames(oGRsnps)
seqlevels(oGRsnps)
## do some acrobatics with seqlevels to perform search
## IF this will give error, otherwise ignore next few lines
oGPsnps = snpsByOverlaps(oDBSnp, oGRsnps)

seqlevels(oGRsnps) = seqlevels(oDBSnp)
seqinfo(oGRsnps) = seqinfo(oDBSnp) 
oGPsnps = snpsByOverlaps(oDBSnp, oGRsnps)

# manually check a few rs numbers if any hits at
# https://www.ncbi.nlm.nih.gov/snp/?term=rs3745946
#IUPAC_CODE_MAP
oGPsnps

## update the SNPs GRanges object
f = findOverlaps(oGRsnps, oGPsnps)
oGRsnps$RefSNP_id = NA_character_
oGRsnps$RefSNP_id[queryHits(f)] = oGPsnps$RefSNP_id[subjectHits(f)]
oGRsnps.snploc = oGRsnps
## match the names with deepSNV results table
gr = GRanges(as.character(dfDeepSNV$chr), IRanges(start=dfDeepSNV$pos, width = 1))
seqlevelsStyle(oGRsnps.snploc) = 'UCSC'
i = findOverlaps(oGRsnps.snploc, gr, type='end')
identical(queryHits(i), subjectHits(i))
w = which(queryHits(i) != subjectHits(i))
identical(dfDeepSNV$pos[subjectHits(i)], end(oGRsnps[queryHits(i)]))
dfDeepSNV$RefSNP_id = oGRsnps.snploc$RefSNP_id

## location of the variants 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
oTxdb = TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevelsStyle(oTxdb)
seqlevelsStyle(oGRsnps) = 'UCSC'
#seqlevelsStyle(oTxdb) = 'UCSC'
genome(oGRsnps)
genome(oTxdb)
genome(oGRsnps) = genome(oTxdb)

library(VariantAnnotation)
oGRlocations = locateVariants(oGRsnps, oTxdb, AllVariants())
## as comparison is unstranded, remove those matches that do not
## have the same gene id as the snp object, as we are reporting
## a match from the opposite strand as well
oGRlocations = oGRlocations[oGRlocations$GENEID %in% oGRsnps$ENTREZID]
library(org.Hs.eg.db)
## summarize each snp by location
dfSymbol = AnnotationDbi::select(org.Hs.eg.db, keys = oGRlocations$GENEID, columns = 'SYMBOL',
                                keytype = 'ENTREZID')
## check if names match with deepsnv results
table(dfDeepSNV$ENTREZID %in% oGRlocations$GENEID)
table(dfDeepSNV$SYMBOL %in% dfSymbol$SYMBOL)

#lSnpLoc = split(oGRlocations$LOCATION, oGRlocations$QUERYID)
lSnpLoc = split(oGRlocations$LOCATION, dfSymbol$SYMBOL)
length(lSnpLoc)
dfLocations = data.frame(sapply(lSnpLoc, function(x) table(x)))

write.csv(dfDeepSNV, file='temp/testing_RefSNP.csv')
write.csv(dfLocations, file='results/testing_Locations.csv')
