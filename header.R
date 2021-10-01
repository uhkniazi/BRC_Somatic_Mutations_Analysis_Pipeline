# File: header.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: global variables
# Date: 22/9/2021


## variables
#g_pid = 22
#g_did = 43
gcswd = getwd()
gcRemoteDir = "/run/user/1000/gvfs/sftp:host=10.202.64.29,user=k1625253/users/k1625253/brc_scratch/Data/ProjectsData/"

p.old = par()

###### utility functions

f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}


# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$adj.P.Val < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logFC, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

###################### modification of deepSNV summary/.significantSNV function
## as it gives error when generating VCF output
# summary(oSnv, value='VCF')
# Error in metadata(v)$header@header$META["date", 1] <- paste(Sys.time()) : 
#   incorrect number of subscripts on matrix
# use getMethod to find this function
# getMethod('summary', 'deepSNV')
# function source copied from
# https://rdrr.io/github/gerstung-lab/deepSNV/src/R/deepSNV-functions.R
# or see code here
# deepSNV:::.significantSNV
myDeepSNVSummary = function (deepSNV, sig.level = 0.05, adjust.method = "bonferroni", 
                             fold.change = 1, value = "data.frame") 
{
  if (is.null(adjust.method)) 
    q = deepSNV@p.val
  else q = p.adjust(deepSNV@p.val, method = adjust.method)
  CV.control <- consensusSequence(deepSNV, vector = TRUE)
  f <- RF(deepSNV@test, total = T)
  g <- RF(deepSNV@control, total = T)
  d <- f - g
  fc <- pmax(f, g)/pmin(f, g)
  var.d <- f/rowSums(deepSNV@test) + g/rowSums(deepSNV@control)
  if (is.null(fold.change)) 
    fccond <- TRUE
  else fccond <- fc > fold.change
  cond = q <= sig.level & !is.na(q) & fccond
  cond.idx <- which(cond, arr.ind = TRUE)
  cond.rows <- cond.idx[, 1]
  cond.cols <- cond.idx[, 2]
  lengths = deepSNV@regions$stop - deepSNV@regions$start + 
    1
  beg = cumsum(c(1, lengths[-length(lengths)]))
  end = cumsum(lengths)
  reg_ix = as.numeric(sapply(cond.rows, function(i) which(end >= 
                                                            i)[1]))
  pos = deepSNV@regions$start[reg_ix] + cond.rows - beg[reg_ix]
  table <- data.frame(chr = deepSNV@regions$chr[reg_ix], pos = pos, 
                      ref = CV.control[cond.rows], var = factor(cond.cols, 
                                                                levels = 1:5, labels = colnames(f)), p.val = q[cond], 
                      freq.var = d[cond], sigma2.freq.var = var.d[cond], n.tst.fw = deepSNV@test[, 
                                                                                                 1:5][cond], cov.tst.fw = rowSums(deepSNV@test[cond.rows, 
                                                                                                                                               1:5, drop = FALSE]), n.tst.bw = deepSNV@test[, 6:10][cond], 
                      cov.tst.bw = rowSums(deepSNV@test[cond.rows, 6:10, drop = FALSE]), 
                      n.ctrl.fw = deepSNV@control[, 1:5][cond], cov.ctrl.fw = rowSums(deepSNV@control[cond.rows, 
                                                                                                      1:5, drop = FALSE]), n.ctrl.bw = deepSNV@control[, 
                                                                                                                                                       6:10][cond], cov.ctrl.bw = rowSums(deepSNV@control[cond.rows, 
                                                                                                                                                                                                          6:10, drop = FALSE]))
  if (ncol(deepSNV@regions) > 3) {
    table <- cbind(table, deepSNV@regions[reg_ix, -c(1, 2, 
                                                     3)])
    colnames(table)[15 + 1:(ncol(deepSNV@regions) - 3)] = colnames(deepSNV@regions)[-c(1, 
                                                                                       2, 3)]
  }
  if (!is.null(adjust.method)) {
    table$raw.p.val <- deepSNV@p.val[cond]
  }
  rownames(table) <- NULL
  if (value == "data.frame") {
    o <- order(table$chr, table$pos)
    table <- table[o, ]
    return(table)
  }
  else {
    if (nrow(table) == 0) 
      return(NULL)
    isDel <- table$var == "-"
    isIns <- table$ref == "-"
    if (length(deepSNV@files) == 2) 
      samples <- sub("/.+/", "", deepSNV@files)
    else samples <- c("test", "control")
    v = VCF(rowRanges = GRanges(table$chr, 
                                IRanges(table$pos - (isDel | isIns), width = 1 + (isDel | isIns)), ), 
            fixed = DataFrame(REF = DNAStringSet(paste(ifelse(isDel | 
                                                                isIns, as.character(CV.control[cond.rows - 1]), 
                                                              ""), sub("-", "", table$ref), sep = "")), ALT = do.call(DNAStringSetList, 
                                                                                                                      as.list(paste(ifelse(isDel | isIns, as.character(CV.control[cond.rows - 
                                                                                                                                                                                    1]), ""), sub("-", "", table$var), sep = ""))), 
                              QUAL = round(-10 * log10(table$raw.p.val)), FILTER = "PASS"), 
            info = DataFrame(VF = table$freq.var, VFV = table$sigma2.freq.var), 
            geno = SimpleList(FW = cbind(table$n.tst.fw, table$n.ctrl.fw), 
                              BW = cbind(table$n.tst.bw, table$n.ctrl.bw), 
                              DFW = cbind(table$cov.tst.fw, table$cov.ctrl.fw), 
                              DBW = cbind(table$cov.tst.bw, table$cov.ctrl.bw)), 
            metadata = list(header = scanVcfHeader(system.file("extdata", 
                                                               "deepSNV.vcf", package = "deepSNV"))), colData = DataFrame(samples = 1:length(samples), 
                                                                                                                          row.names = samples), collapsed = TRUE)
    metadata(v)$header@samples <- samples
    #metadata(v)$header@header$META["date"] <- paste(Sys.time())
    return(sort(v))
  }
}
