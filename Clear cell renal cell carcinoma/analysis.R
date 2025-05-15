## DATA PREPARATION
## ================

## load libraries
library(MIRit)
library(patchwork)

## define results folder
resFolder <- "results"
plotsFolder <- "results/plots"

## set seed for reproducible results
set.seed(1234)

## load data
genExpr <- readRDS("data/gene_expression_matrix.RDS")
mirExpr <- readRDS("data/miRNA_expression_matrix.RDS")

## define samples metadata
genMeta <- read.csv("data/gene_metadata.csv")
mirMeta <- read.csv("data/mirna_metadata.csv")
identical(genMeta$identifier.ch1, mirMeta$identifier.ch1)
identical(genMeta$histology.ch1, mirMeta$histology.ch1)
genMeta$condition <- "HC"
genMeta$condition[genMeta$histology.ch1 == "Clear Cell"] <- "ccRCC"
genMeta$ID <- paste(genMeta$identifier.ch1, genMeta$condition, sep = " - ")
meta <- data.frame(primary = genMeta$ID,
                   mirnaCol = mirMeta$geo_accession,
                   geneCol = genMeta$geo_accession,
                   condition = genMeta$condition,
                   individual = genMeta$identifier.ch1)

## create a MirnaExperiment object
experiment <- MirnaExperiment(mirnaExpr = mirExpr,
                              geneExpr = genExpr,
                              samplesMetadata = meta)


## DIFFERENTIAL EXPRESSION ANALYSIS
## ================================

## display expression variability
plotDimensions(experiment, "microRNA", condition = "condition")
plotDimensions(experiment, "genes", condition = "condition")

## perform miRNA and gene differential expression
experiment <- performGeneDE(experiment,
                            group = "condition",
                            contrast = "ccRCC-HC",
                            design = ~ condition + individual,
                            method = "limma",
                            logFC = log2(1.5),
                            useArrayWeights = TRUE,
                            eBayes.args = list(robust = TRUE))
experiment <- performMirnaDE(experiment,
                             group = "condition",
                             contrast = "ccRCC-HC",
                             design = ~ condition + individual,
                             method = "limma",
                             logFC = log2(1.5),
                             useArrayWeights = TRUE,
                             eBayes.args = list(robust = TRUE))

## save DE results
write.csv(geneDE(experiment), "results/deg.csv")
write.csv(mirnaDE(experiment), "results/dem.csv")

## display differential expression results
plotVolcano(experiment, "microRNA")
plotVolcano(experiment, "genes")


## FUNCTIONAL ENRICHMENT ANALYSIS
## ==============================

## perform gene functional enrichment
ora <- enrichGenes(experiment, method = "ORA", database = "GO", pCutoff = 0.1)
gsea <- enrichGenes(experiment, method = "GSEA", database = "GO", pCutoff = 0.1)


## DATA INTEGRATION
## ================

## load a local copy of miRTarBase v10
mtb <- read.csv("data/miRTarBase_MTI.csv")

## retrieve miRNA target genes
experiment <- getTargets(experiment, score = "High", local = mtb)

## integrate miRNA and gene expression levels
experiment <- mirnaIntegration(experiment, pCutoff = 0.1)
write.csv(integration(experiment), "results/correlation.csv")

## perform the enrichment of integrated targets
netGO <- enrichTargets(experiment, database = "GO", pCutoff = 0.1, minSize = 20, maxSize = 1000)
netDO <- enrichTargets(experiment, database = "DO", pCutoff = 0.1)

## save the MirnaExperiment object
saveRDS(experiment, "results/experiment.RDS")


## FIGURES
## =======

## Panel 1
## -------

## set the color scale
sc <- c(ccRCC = "#1a80bb", HC = "#b8b8b8")

## generate the plots
mdsGene <- plot.MDS(experiment, "genes", colorScale = sc, pointAlpha = 0.8)
mdsMirna <- plot.MDS(experiment, "microRNA", colorScale = sc, pointAlpha = 0.8)
volGene <- plot.Volcano(experiment, "genes")
volMirna <- plot.Volcano(experiment, "microRNA")
mir1 <- plot.Barplot(experiment, "hsa-miR-142-3p", "microRNA", colorScale = sc,
                     brackets = list(c("ccRCC", "HC", "***")),
                     jitter = TRUE, title = "miR-142-3p")
mir2 <- plot.Barplot(experiment, "hsa-miR-200c", "microRNA", colorScale = sc,
                     brackets = list(c("ccRCC", "HC", "****")),
                     jitter = TRUE, title = "miR-200c")
mir3 <- plot.Barplot(experiment, "hsa-miR-210", "microRNA", colorScale = sc,
                     brackets = list(c("ccRCC", "HC", "****")),
                     jitter = TRUE, title = "miR-210")
mir4 <- plot.Barplot(experiment, "hsa-miR-141", "microRNA", colorScale = sc,
                     brackets = list(c("ccRCC", "HC", "****")),
                     jitter = TRUE, title = "miR-141")
cor1 <- plot.Correlation(experiment, "hsa-miR-142-3p", "HOXA7", colorScale = sc,
                         corMethod = "spearman", lineCol = "grey20", pointAlpha = 0.8)
cor2 <- plot.Correlation(experiment, "hsa-miR-200c", "F3", colorScale = sc,
                         corMethod = "spearman", lineCol = "grey20", pointAlpha = 0.8)
cor3 <- plot.Correlation(experiment, "hsa-miR-27a", "SLC13A3", colorScale = sc,
                         corMethod = "spearman", lineCol = "grey20", pointAlpha = 0.8)

## compose the panel
pan1top <- (mdsGene + mdsMirna + volGene + volMirna + plot_layout(nrow = 1)) &
    theme(legend.position = "none")
pan1mid <- (mir1 + mir2 + mir3 + mir4 + plot_layout(nrow = 1)) &
    theme(legend.position = "none")
pan1mid <- pan1mid & ylab(expression(log[2]~expression))
pan1bot <- (cor1 + cor2 + cor3 + plot_layout(nrow = 1, guides = "collect")) &
    theme(legend.position = "bottom")
pan1bot <- pan1bot & guides(color = guide_legend(override.aes = list(alpha = 1) ) )
pan1 <- pan1top / pan1mid / pan1bot

## adjust the panel
pan1 <- pan1 + plot_layout(heights = c(1, 1, 1.25))
pan1 <- pan1 +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

## save the panel
ggsave("figures/ccRCC_de.pdf", pan1, width = 180, height = 150, units = "mm",
       scale = 1.5, device = cairo_pdf)


## Panel 2
## -------

## create the plots
net <- plot.Network(experiment, minCorr = -0.7) +
    theme(legend.position = "none")
netEnr <- plot.DirEnrichment(netGO, sel = list(c(1:4),
                                               c(1, 2, 4, 5, 6, 7, 8, 15, 17,
                                                 18, 22, 23, 29, 30, 37, 43)))

## create the panel
pan2 <- net + netEnr + plot_layout(widths = c(1, 0.6))
pan2 <- pan2 +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

## save the panel
ggsave("figures/ccRCC_net.pdf", pan2, width = 180, height = 130, units = "mm",
       scale = 1.5, device = cairo_pdf)


