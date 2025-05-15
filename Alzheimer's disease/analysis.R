## DATA PREPARATION
## ================

## load libraries
library(MIRit)
library(qvalue)
library(patchwork)
library(ggpubr)

## define results folder
resFolder <- "results"
plotsFolder <- "results/plots"

## set seed for reproducible results
set.seed(1234)

## load data
genExpr <- read.csv("data/GSE150696/normalized_expression.csv", row.names = "X")
mirExpr <- read.table("data/GSE63501/counts.csv", header = TRUE, sep = ",", row.names = "X")

## define samples metadata
genMeta <- read.csv("data/GSE150696/metadata.csv", row.names = "X")
mirMeta <- read.csv("data/GSE63501/metadata.csv", row.names = "X")
genMeta$sample <- gsub("/", ".", genMeta$sample)
genMeta$condition[genMeta$condition == "CTRL"] <- "HC"
mirMeta$classification[mirMeta$classification == "Control"] <- "HC"
mirMeta$primary <- mirMeta$sample
genMeta$primary <- genMeta$sample
colnames(mirMeta)[colnames(mirMeta) == "sample"] <- "mirnaCol"
colnames(mirMeta)[colnames(mirMeta) == "classification"] <- "condition"
colnames(genMeta)[colnames(genMeta) == "sample"] <- "geneCol"
meta <- merge(mirMeta, genMeta, by = c("primary", "condition"), all = TRUE)

## filter miRNA counts
keep <- rowSums(mirExpr >= 10) >= min(table(mirMeta$condition))
table(keep)
mirExpr <- mirExpr[keep, ]

## create a MirnaExperiment object
experiment <- MirnaExperiment(mirnaExpr = mirExpr,
                              geneExpr = genExpr,
                              samplesMetadata = meta,
                              pairedSamples = FALSE)


## DIFFERENTIAL EXPRESSION ANALYSIS
## ================================

## display expression variability
plotDimensions(experiment, "microRNA", condition = "condition")
plotDimensions(experiment, "genes", condition = "condition")

## perform miRNA and gene differential expression
experiment <- performGeneDE(experiment,
                            group = "condition",
                            contrast = "AD-HC",
                            design = ~ condition + sex,
                            method = "limma",
                            useArrayWeights = TRUE)
experiment <- performMirnaDE(experiment,
                             group = "condition",
                             contrast = "AD-HC",
                             design = ~ condition,
                             method = "voom",
                             useVoomWithQualityWeights = TRUE,
                             pCutoff = 0.1)

## correct p-values with the Storey's q-value approach
experiment@mirnaDE$data$adj.P.Val <- qvalue(experiment@mirnaDE$data$P.Value)$qvalues
experiment@mirnaDE$significant <- experiment@mirnaDE$data$ID[experiment@mirnaDE$data$adj.P.Val < 0.1]

## save DE results
write.csv(geneDE(experiment), "results/deg.csv")
write.csv(mirnaDE(experiment), "results/dem.csv")

## display differential expression results
plotVolcano(experiment, "microRNA")
plotVolcano(experiment, "genes")


## FUNCTIONAL ENRICHMENT ANALYSIS AND MIRNA-SNPs
## =============================================

## perform disease-SNPs association
searchDisease("alzheimer")
snp <- findMirnaSNPs(experiment, "Alzheimer disease")
write.csv(snp, "results/snp.csv")

## perform gene functional enrichment
ora <- enrichGenes(experiment, method = "ORA", database = "GO")
gsea <- enrichGenes(experiment, method = "GSEA", database = "KEGG")


## DATA INTEGRATION
## ================

## load a local copy of miRTarBase v10
mtb <- read.csv("data/miRTarBase_MTI.csv")

## retrieve miRNA target genes
experiment <- getTargets(experiment, score = "High", local = mtb)

## integrate miRNA and gene expression levels
experiment <- mirnaIntegration(experiment, associationMethod = "boschloo", pCutoff = 0.1)
write.csv(integration(experiment), "results/association.csv")

## enrichment of integrated targets
netGO <- enrichTargets(experiment, database = "GO")

## save the MirnaExperiment object
saveRDS(experiment, "results/experiment.RDS")


## FIGURES
## =======

## Panel 1
## -------

## set the color scale
sc <- c(AD = "#1a80bb", HC = "#b8b8b8")

## generate the plots
mdsGene <- plot.MDS(experiment, "genes", colorScale = sc, pointAlpha = 0.8)
mdsMirna <- plot.MDS(experiment, "microRNA", colorScale = sc, pointAlpha = 0.8)
volGene <- plot.Volcano(experiment, "genes")
volMirna <- plot.Volcano(experiment, "microRNA")
trPlot <- readRDS("figures/trackplot.RDS")
netPlot <- plot.Enrichment(netGO$downregulated,
                           sel = c(4, 5, 6, 8, 9, 11, 14, 15, 19, 28,
                                   82, 151, 155, 185, 202, 230, 306, 317))

## extract the legend from the first plot
lg <- get_legend(mdsGene)

## adjust the plots
mdsGene <- mdsGene + theme(legend.position = "none")
mdsMirna <- mdsMirna + theme(legend.position = "none")
volGene <- volGene + theme(legend.position = "none")
volMirna <- volMirna + theme(legend.position = "none")
netPlot <- netPlot + scale_size_continuous(name = "Overlap",
                                           breaks = c(0.1, 0.2, 0.3),
                                           range = c(2, 5))

## compose the panel
pan1top <- (mdsGene + mdsMirna + volGene + volMirna) +
    plot_layout(nrow = 1) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))
pan1bot <- (free(trPlot) + free(netPlot)) +
    plot_layout(nrow = 1, widths = c(1, 2)) +
    plot_annotation(tag_levels = list(c("E", "F"))) &
    theme(plot.tag = element_text(size = 14, face = "bold"))

## add the custom legend
pan1 <- ggarrange(pan1top, pan1bot, lg, nrow = 3, heights = c(5, 10, 1))

## save the panel
ggsave("figures/AD_de.pdf", pan1, width = 180, height = 150, units = "mm",
       scale = 1.5, device = cairo_pdf)





