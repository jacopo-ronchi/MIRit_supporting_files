## PREPROCESSING OF MICROARRAY DATA
## ================================

## load required libraries
library(oligo)
library(affycoretools)
library(hta20transcriptcluster.db)
library(limma)

## define project folders
dataDir <- "data/GSE150696"
plotsFolder <- "results/plots"
dir.create(plotsFolder, recursive = TRUE)

## load sample metadata
mesMeta <- read.csv(paste(dataDir, "metadata.csv", sep = "/"),
                    row.names = "X")

## import CEL files
celpath <- paste(dataDir, "GSE150696_RAW", sep = "/")
celFiles <- list.files(celpath, full.names = TRUE, recursive = TRUE)
celFiles <- celFiles[endsWith(celFiles, ".cel")]
data <- read.celfiles(celFiles)

## rename samples
colnames(data) <- mesMeta$sample

## normalize data with RMA algorithm
dataNorm <- rma(data, target = "core")

## annotate probes at the transcript level
eset <- getMainProbes(dataNorm, level = "core")
eset <- annotateEset(eset, hta20transcriptcluster.db)

## remove probes without proper annotation
eset <- eset[!is.na(featureData(eset)@data$SYMBOL)]

## filter lowly expressed probes
rMed <- rowMedians(exprs(eset))
hist(rMed, 100, col = "lightblue", border = "royalblue", freq = FALSE,
     main = "Distribution of probe intensities", xlab = "Median intensities")
thr <- 6
abline(v = thr, col = "darkgreen", lwd = 2)
intPlot <- recordPlot()
saveRDS(intPlot, paste(plotsFolder, "gene_probe_intensities.RDS", sep = "/"))
aboveThr <- apply(eset, 1, function(x) sum(x > thr) >= min(table(mesMeta$condition)))
table(aboveThr)
eset <- subset(eset, aboveThr)

## average probes targeting the same gene
geneConv <- featureData(eset)@data$SYMBOL[match(featureData(eset)@data$PROBEID,
                                                   rownames(eset))]
eset <- avereps(eset, ID = geneConv)
write.csv(eset, "normalized_expression.csv")

