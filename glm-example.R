options(width=75)
options(digits=4)
rawCounts <- data.frame(
  mrna_normox_a = c( 195, 200, 310, rep(200, 4) ),
  mrna_normox_b = c( 200, 210, 295, rep(200, 4) ),
  mrna_hypox_a  = c( 210,  95, 295, rep(200, 4) ),
  mrna_hypox_b  = c( 205, 100, 300, rep(200, 4) ),
  ribo_normox_a = c( 305, 295, 105, rep(200, 4) ),
  ribo_normox_b = c( 310, 310, 100, rep(200, 4) ),
  ribo_hypox_a  = c( 290, 150, 395, rep(200, 4) ),
  ribo_hypox_b  = c( 315, 155, 410, rep(200, 4) ),
  row.names=c("A", "B", "C", lapply(1:4, function(n) { sprintf("etc%d", n) })))

rawCounts[1:3,]

library("DESeq")

conditions <- data.frame(
  biol     = factor(c("mrna", "mrna", "mrna", "mrna", "ribo", "ribo", "ribo", "ribo"), levels=c("mrna", "ribo")),
  oxia     = factor(c("norm", "norm", "hypo", "hypo", "norm", "norm", "hypo", "hypo"), levels=c("norm", "hypo")),
  row.names= colnames(rawCounts)
  )
conditions

countData <- newCountDataSet(rawCounts, conditions)
countData <- estimateSizeFactors(countData)
countData <- estimateDispersions(countData, method = "blind", sharingMode = "fit-only" )

glmNoChg  <- fitNbinomGLMs( countData, count ~ biol - 1 )
format(glmNoChg[1:3,])
glmNoChg$mrnaCounts <- 2**(glmNoChg$biolmrna)
glmNoChg$riboCounts <- 2**(glmNoChg$biolribo)
format(glmNoChg[1:3,])

glmNoTrl <- fitNbinomGLMs( countData, count ~ biol + oxia - 1 )
glmNoTrl$mrnaNormC <- 2**(glmNoTrl$biolmrna)
glmNoTrl$mrnaHypoC <- 2**(glmNoTrl$biolmrna + glmNoTrl$oxiahypo)
glmNoTrl$riboNormC <- 2**(glmNoTrl$biolribo)
glmNoTrl$riboHypoC <- 2**(glmNoTrl$biolribo + glmNoTrl$oxiahypo)
glmNoTrl$hypoxChange <- 2**(glmNoTrl$oxiahypo)

format(glmNoTrl[1:3,])

pNoTrlVsNoChg <- nbinomGLMTest( glmNoTrl, glmNoChg )
pNoTrlVsNoChg[1:3]

glmFull  <- fitNbinomGLMs( countData, count ~ biol * oxia - 1 )
glmFull$mrnaNormC <- 2**(glmFull$biolmrna)
glmFull$mrnaHypoC <- 2**(glmFull$biolmrna + glmFull$oxiahypo)
glmFull$riboNormC <- 2**(glmFull$biolribo)
glmFull$riboHypoC <- 2**(glmFull$biolribo + glmFull$oxiahypo + glmFull$"biolribo:oxiahypo")
glmFull$hypoxMrnaChg <- 2**(glmFull$oxiahypo)
glmFull$hypoxTEChg <- 2**(glmFull$"biolribo:oxiahypo")

format(glmFull[1:3,])

pFullVsNoChg  <- nbinomGLMTest( glmFull, glmNoChg )
pFullVsNoTrl  <- nbinomGLMTest( glmFull, glmNoTrl )

pFullVsNoChg[1:3]
pFullVsNoTrl[1:3]
