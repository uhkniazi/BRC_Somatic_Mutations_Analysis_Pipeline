# File: 09.1_annotateVariantTable.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Annotate the identified somatic mutations using shearwaterML.R script
# Date: 14/12/2021

source('header.R')
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)

cvFiles = list.files('results/shearwaterML/')

ldfSNPs = vector(mode = 'list', length = length(cvFiles))
names(ldfSNPs) = cvFiles
ldfSNPLocs = vector(mode = 'list', length = length(cvFiles))
names(ldfSNPLocs) = cvFiles
setwd('results/shearwaterML/')

for(iIndex in seq_along(cvFiles)){
  
  ## load the relevant files
  dfDeepSNV = read.csv(cvFiles[iIndex], header=T, row.names=1, stringsAsFactors = F)
  dim(dfDeepSNV)
  dfDeepSNV = dfDeepSNV[dfDeepSNV$p.adj < 0.01, ] 
  dim(dfDeepSNV)
  ## use the snplocs database
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
  
  # seqlevels(oGRsnps) = seqlevels(oDBSnp)
  # seqinfo(oGRsnps) = seqinfo(oDBSnp) 
  # oGPsnps = snpsByOverlaps(oDBSnp, oGRsnps)
  # 
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
  oTxdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  seqlevelsStyle(oTxdb)
  seqlevelsStyle(oGRsnps) = 'UCSC'
  #seqlevelsStyle(oTxdb) = 'UCSC'
  genome(oGRsnps)
  genome(oTxdb)
  genome(oGRsnps) = genome(oTxdb)
  
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
  
  ldfSNPs[[iIndex]] = dfDeepSNV
  ldfSNPLocs[[iIndex]] = dfLocations
  # write.csv(dfDeepSNV, file='temp/testing_RefSNP.csv')
  # write.csv(dfLocations, file='results/testing_Locations.csv')
}

sapply(ldfSNPs, dim)
lapply(ldfSNPs, head)

n = names(ldfSNPs)
setwd(gcswd)
sapply(seq_along(n), function(x){
  write.csv(ldfSNPs[[x]], paste0('results', '/SNPs', n[x]))
  write.csv(ldfSNPLocs[[x]], paste0('results', '/Locations', n[x]))
})
