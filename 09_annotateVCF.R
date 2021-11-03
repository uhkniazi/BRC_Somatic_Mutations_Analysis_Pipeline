# File: 09_annotateVCF.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Annotate the identified somatic mutations
# Date: 2/11/2021

source('header.R')

library(VariantAnnotation)

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
# oGRsnps = GRanges(as.character(seqnames(oGRsnps)),
#                   ranges(oGRsnps),
#                   )

seqlevels(oGRsnps) = seqlevels(oDBSnp)
seqinfo(oGRsnps) = seqinfo(oDBSnp) 
oGPsnps = snpsByOverlaps(oDBSnp, oGRsnps)

# manually check a few rs numbers if any hits at
# https://www.ncbi.nlm.nih.gov/snp/?term=rs3745946
#IUPAC_CODE_MAP
oGPsnps

## location of the variants 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
oTxdb = TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevelsStyle(oTxdb)
oGRsnps = rowRanges(oVcf)
seqlevelsStyle(oGRsnps)
seqlevelsStyle(oTxdb) = 'UCSC'
oGRlocations = locateVariants(oGRsnps, oTxdb, AllVariants())

## summarize each snp by location
lSnpLoc = split(oGRlocations$LOCATION, oGRlocations$QUERYID)
sapply(lSnpLoc, function(x) table(x))
