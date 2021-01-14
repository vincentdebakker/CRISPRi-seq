# Vincent de Bakker
# 2020 NOV 10
#
# Adapted from DESeq2::makeExampleDESeqDataSet()
# Love et al. (2014) Genome Biology 15(12):550

makeExampleDESeqDataSet_vdb <- function (n = 1000, m = 12, betaMU = 0, betaSD = 0, interceptMean = 4, interceptSD = 2, 
          dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1, 
                                                                 m)) 
{
  beta <- cbind(rnorm(n, interceptMean, interceptSD), rnorm(n, 
                                                            betaMU, betaSD))
  dispersion <- dispMeanRel(2^(beta[, 1]))
  colData <- DataFrame(condition = factor(rep(c("A", "B"), 
                                              times = c(ceiling(m/2), floor(m/2)))))
  x <- if (m > 1) {
    stats::model.matrix.default(~colData$condition)
  }
  else {
    cbind(rep(1, m), rep(0, m))
  }
  mu <- t(2^(x %*% t(beta)) * sizeFactors)
  countData <- matrix(rnbinom(m * n, mu = mu, size = 1/dispersion), 
                      ncol = m)
  mode(countData) <- "integer"
  colnames(countData) <- paste("sample", 1:m, sep = "")
  rowRanges <- GRanges("1", IRanges(start = (1:n - 1) * 100 + 
                                      1, width = 100))
  names(rowRanges) <- paste0("gene", 1:n)
  design <- if (m > 1) {
    as.formula("~ condition", env = .GlobalEnv)
  }
  else {
    as.formula("~ 1", env = .GlobalEnv)
  }
  object <- DESeqDataSetFromMatrix(countData = countData, 
                                   colData = colData, design = design, rowRanges = rowRanges)
  trueVals <- DataFrame(trueIntercept = beta[, 1], trueBeta = beta[, 
                                                                   2], trueDisp = dispersion)
  mcols(trueVals) <- DataFrame(type = rep("input", ncol(trueVals)), 
                               description = c("simulated intercept values", "simulated beta values", 
                                               "simulated dispersion values"))
  mcols(object) <- cbind(mcols(object), trueVals)
  return(object)
}