## PREPARE GENE EXPRESSION DATA FOR ccRCC DATSET
## =============================================

## load libraries
library(limma)
library(GEOquery)

## define results folder
resFolder <- "results"
plotsFolder <- "results/plots"
dir.create(plotsFolder, recursive = TRUE)

## load GSE16441 gene expression data
fList <- list.files("data/raw/", full.names = TRUE)[1:34]
raw <- read.maimages(fList,
                     source = "agilent",
                     green.only = TRUE,
                     other.columns = "gIsWellAboveBG")

## explore expression distributions in samples
boxplot(log2(raw$E))

## background correction and quantile normalization
data <- backgroundCorrect(raw, method = "normexp")
data <- normalizeBetweenArrays(data, method = "quantile")

## filter control probes and not annotated probes
control <- data$genes$ControlType == 1L
noID <- is.na(data$genes$GeneName)

## filter probes not above background on at least 17 arrays
aboveBG <- rowSums(data$other$gIsWellAboveBG > 0) >= 17

## maintain the desired probes
exprMatrix <- data[!control & !noID & aboveBG, ]
dim(exprMatrix)
exprMatrix$genes <- exprMatrix$genes[, c("ProbeName", "GeneName")]

## average the probes with the same targets
exprMatrix <- avereps(exprMatrix, ID = exprMatrix$genes$GeneName)


## DEFINE SAMPLES AND COVARIATES
## =============================

## extract samples metadata from GEO
geo <- getGEO("GSE16441")
geo <- geo[[grep("GPL6480", names(geo))]]

## define samples metadata
meta <- pData(geo)
meta <- meta[, c(1, 2, 57, 58)]

## rename samples
colnames(exprMatrix$E) <- meta$geo_accession
rownames(meta) <- meta$geo_accession

## save metadata
write.csv(meta, "data/gene_metadata.csv")


## EXPRESSION FILTERING
## ====================

## calculate median intensities for each probe
rMed <- rowMedians(exprMatrix$E)

## plot the distribution of median intensities
hist(rMed, 100, col = "lightblue", border = "royalblue", freq = FALSE,
     main = "Distribution of probe intensities", xlab = "Median intensities")

## define a suitable threshold
thr <- 8
abline(v = thr, col = "darkgreen", lwd = 2)

## save the plot
intPlot <- recordPlot()
saveRDS(intPlot, paste(plotsFolder, "gene_probe_intensities.RDS", sep = "/"))

## remove probes where the intensity is under the threshold in more than 17 samples
aboveThr <- apply(exprMatrix$E, 1, function(x) sum(x > thr) >= 17)
table(aboveThr)
exprFlt <- exprMatrix$E[aboveThr, ]
saveRDS(exprFlt, "data/gene_expression_matrix.RDS")

