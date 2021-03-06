options(width=75)
options(digits=4)
rawCounts <- data.frame(
  mrna_normox_a = c( 195, 200, 310, 100, 150 + rbinom(17, 100, 0.5) ),
  mrna_normox_b = c( 200, 210, 295, 105, 150 + rbinom(17, 100, 0.5)),
  mrna_hypox_a  = c( 210,  95, 295, 195, 150 + rbinom(17, 100, 0.5) ),
  mrna_hypox_b  = c( 205, 100, 300, 205, 150 + rbinom(17, 100, 0.5) ),
  ribo_normox_a = c( 305, 295, 105, 145, 150 + rbinom(17, 100, 0.5) ),
  ribo_normox_b = c( 310, 310, 100, 150, 150 + rbinom(17, 100, 0.5) ),
  ribo_hypox_a  = c( 290, 150, 395, 440, 150 + rbinom(17, 100, 0.5) ),
  ribo_hypox_b  = c( 315, 155, 410, 455, 150 + rbinom(17, 100, 0.5)),
  row.names=c("A", "B", "C", "D", lapply(1:17, function(n) { sprintf("etc%d", n) })))

rawCounts[1:4,]

library("DESeq")

conditions <- data.frame(
  biol     = factor(c("mrna", "mrna", "mrna", "mrna", "ribo", "ribo", "ribo", "ribo"), levels=c("mrna", "ribo")),
  oxia     = factor(c("norm", "norm", "hypo", "hypo", "norm", "norm", "hypo", "hypo"), levels=c("norm", "hypo")),
  bioloxia = factor(c("norm", "norm", "norm", "norm", "norm", "norm", "ribohypo", "ribohypo"), levels=c("norm", "ribohypo")),
  row.names= colnames(rawCounts)
  )
conditions

countData <- newCountDataSet(rawCounts, conditions)
countData <- estimateSizeFactors(countData)
countData <- estimateDispersions(countData, method = "pooled", sharingMode="fit-only", fitType="local")

glmNoChg  <- fitNbinomGLMs( countData, count ~ biol - 1 )
format(glmNoChg[1:4,])
glmNoChg$mrnaCounts <- 2**(glmNoChg$biolmrna)
glmNoChg$riboCounts <- 2**(glmNoChg$biolribo)
format(glmNoChg[1:4,])

glmNoTrl <- fitNbinomGLMs( countData, count ~ biol + oxia - 1 )

glmNoTrl$mrnaNormC <- 2**(glmNoTrl$biolmrna)
glmNoTrl$mrnaHypoC <- 2**(glmNoTrl$biolmrna + glmNoTrl$oxiahypo)
glmNoTrl$riboNormC <- 2**(glmNoTrl$biolribo)
glmNoTrl$riboHypoC <- 2**(glmNoTrl$biolribo + glmNoTrl$oxiahypo)
glmNoTrl$hypoxChange <- 2**(glmNoTrl$oxiahypo)
format(glmNoTrl[1:4,])

pNoTrlVsNoChg <- nbinomGLMTest( glmNoTrl, glmNoChg )
pNoTrlVsNoChg[1:4]

glmFull  <- fitNbinomGLMs( countData, count ~ biol * oxia - 1 )
glmFull$mrnaNormC <- 2**(glmFull$biolmrna)
glmFull$mrnaHypoC <- 2**(glmFull$biolmrna + glmFull$oxiahypo)
glmFull$riboNormC <- 2**(glmFull$biolribo)
glmFull$riboHypoC <- 2**(glmFull$biolribo + glmFull$oxiahypo + glmFull$"biolribo:oxiahypo")
glmFull$hypoxMrnaChg <- 2**(glmFull$oxiahypo)
glmFull$hypoxTEChg <- 2**(glmFull$"biolribo:oxiahypo")

format(glmFull[1:4,])

pFullVsNoChg  <- nbinomGLMTest( glmFull, glmNoChg )
pFullVsNoTrl  <- nbinomGLMTest( glmFull, glmNoTrl )

pFullVsNoChg[1:4]
pFullVsNoTrl[1:4]

glmNoTrx <- fitNbinomGLMs( countData, count ~ biol + bioloxia - 1 )

glmNoTrx$mrnaNormC <- 2**(glmNoTrx$biolmrna)
glmNoTrx$mrnaHypoC <- 2**(glmNoTrx$biolmrna)
glmNoTrx$riboNormC <- 2**(glmNoTrx$biolribo)
glmNoTrx$riboHypoC <- 2**(glmNoTrx$biolribo + glmNoTrx$bioloxiaribohypo)
glmNoTrx$hypoxChange <- 2**(glmNoTrx$bioloxiaribohypo)
options(width=60)
format(glmNoTrx[1:4,])

pFullVsNoTrx <- nbinomGLMTest( glmFull, glmNoTrx )
pNoTrxVsNoChg <- nbinomGLMTest( glmNoTrx, glmNoChg )

pFullVsNoTrx[1:4]
pNoTrxVsNoChg[1:4]

glmFull$pChg <- pFullVsNoChg
glmFull$pTrx <- pFullVsNoTrx
glmFull$pTrl <- pFullVsNoTrl
glmFull$chg <- glmFull$pChg < 0.01
glmFull$trx <- glmFull$pTrx < 0.01
glmFull$trl <- glmFull$pTrl < 0.01
options(width=95)
glmFull[1:4,]
        