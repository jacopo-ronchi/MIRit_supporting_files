## DATA PREPARATION
## ================

## load libraries
library(MIRit)
library(readxl)
library(qvalue)
library(patchwork)
library(ggpubr)

## define results folder
resFolder <- "results"
dir.create(resFolder)

## set seed for reproducible results
set.seed(1234)

## load data
mirnaCounts <- read_xlsx("data/mirna_counts.xlsx")
geneCounts <- read_xlsx("data/gene_counts.xlsx")

## prepare matrices
mirnaCounts <- as.data.frame(mirnaCounts)
rownames(mirnaCounts) <- mirnaCounts$Gene_symbol
mirnaCounts <- mirnaCounts[, 1:7]
geneCounts <- as.data.frame(geneCounts)
rownames(geneCounts) <- geneCounts$Gene_symbol
geneCounts <- geneCounts[, 1:8]

## exclude lowly expressed miRNAs and genes
keepMir <- rowSums(mirnaCounts >= 10) >= 3
table(keepMir)
mirnaCounts <- mirnaCounts[keepMir, ]
keepGen <- rowSums(geneCounts >= 10) >= 3
table(keepGen)
geneCounts <- geneCounts[keepGen, ]

## define samples metadata
meta <- data.frame(primary = colnames(geneCounts),
                   mirnaCol = c(colnames(mirnaCounts)[1:4], NA,
                                colnames(mirnaCounts)[5:7]),
                   geneCol = colnames(geneCounts),
                   disease = c(rep("VCM", 4), rep("ICM", 4)))

## create a MirnaExperiment object
experiment <- MirnaExperiment(mirnaExpr = mirnaCounts,
                              geneExpr = geneCounts,
                              samplesMetadata = meta)


## DIFFERENTIAL EXPRESSION ANALYSIS
## ================================

## display expression variability
plotDimensions(experiment, "microRNA", condition = "disease")
plotDimensions(experiment, "genes", condition = "disease")

## perform miRNA and gene differential expression
experiment <- performGeneDE(experiment,
                            group = "disease",
                            contrast = "VCM-ICM",
                            design = ~ disease,
                            pCutoff = 0.05,
                            method = "DESeq2")
experiment <- performMirnaDE(experiment,
                             group = "disease",
                             contrast = "VCM-ICM",
                             design = ~ disease,
                             pCutoff = 0.05,
                             method = "DESeq2")

## correct p-values with the Storey's q-value approach
experiment@mirnaDE$data$adj.P.Val <- qvalue(experiment@mirnaDE$data$P.Value)$qvalues
experiment@mirnaDE$significant <- experiment@mirnaDE$data$ID[experiment@mirnaDE$data$adj.P.Val < 0.05]
experiment@geneDE$data$adj.P.Val <- qvalue(experiment@geneDE$data$P.Value)$qvalues
experiment@geneDE$significant <- experiment@geneDE$data$ID[experiment@geneDE$data$adj.P.Val < 0.05]

## fix the order of DE results
experiment@mirnaDE$data <- experiment@mirnaDE$data[order(experiment@mirnaDE$data$P.Value), ]
experiment@geneDE$data <- experiment@geneDE$data[order(experiment@geneDE$data$P.Value), ]

## save DE results
write.csv(geneDE(experiment), "results/deg.csv")
write.csv(mirnaDE(experiment), "results/dem.csv")

## display differential expression results
plotVolcano(experiment, "microRNA")
plotVolcano(experiment, "genes")


## FUNCTIONAL ENRICHMENT ANALYSIS
## ==============================

## perform over-representation analysis
ora <- enrichGenes(experiment, method = "ORA", database = "GO", pCutoff = 0.1)

## set seed
set.seed(1234)

## perform GSEA
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

## save the MirnaExperiment object
saveRDS(experiment, "results/experiment.RDS")


## FIGURES
## =======

## Panel 1
## -------

## set the color scale
sc <- c(VCM = "#1a80bb", ICM = "#b8b8b8")

## generate the plots
volGene <- plot.Volcano(experiment, "genes")
volMirna <- plot.Volcano(experiment, "microRNA")
sel <- c(4, 5, 9, 11, 18, 23, 26, 100, 153, 304, 355, 605)
pRid <- plot.Ridgeplot(gsea, sel = sel)
pNet <- plot.Network(experiment, minCorr = -0.75)
cor1 <- plot.Correlation(experiment, "hsa-miR-218-5p", "DDX6", "disease",
                         colorScale = sc, corMethod = "spearman",
                         lineCol = "grey20", pointAlpha = 0.8) +
    scale_color_manual(values = sc, breaks = names(sc))
cor2 <- plot.Correlation(experiment, "hsa-miR-218-5p", "SEMA4A", "disease",
                         colorScale = sc, corMethod = "spearman",
                         lineCol = "grey20", pointAlpha = 0.8)
cor3 <- plot.Correlation(experiment, "hsa-miR-218-5p", "TTC39C", "disease",
                         colorScale = sc, corMethod = "spearman",
                         lineCol = "grey20", pointAlpha = 0.8)
cor4 <- plot.Correlation(experiment, "hsa-miR-218-5p", "NUP210", "disease",
                         colorScale = sc, corMethod = "spearman",
                         lineCol = "grey20", pointAlpha = 0.8)

## extract the legend from the first plot
lg <- get_legend(cor1)

## adjust the plots
pRid <- pRid + coord_cartesian(xlim = c(-5, 5))
volGene <- volGene + theme(legend.position = "none")
volMirna <- volMirna + theme(legend.position = "none")
pNet <- pNet + theme(legend.position = "none")
pNet <- ggplotify::as.ggplot(pNet)
cor1 <- cor1 +
    scale_y_continuous(breaks = c(3000, 4000, 5000)) +
    scale_x_continuous(breaks = c(3000, 5000, 7000)) +
    theme(legend.position = "none")
cor2 <- cor2 +
    scale_x_continuous(breaks = c(3000, 5000, 7000)) +
    theme(legend.position = "none")
cor3 <- cor3 +
    scale_x_continuous(breaks = c(3000, 5000, 7000)) +
    theme(legend.position = "none")
cor4 <- cor4 +
    scale_x_continuous(breaks = c(3000, 5000, 7000)) +
    theme(legend.position = "none")

## compose the panel
lay <- "AAACCCCCCCCCCCC
        BBBCCCCCCCCCCCC
        DDDDDDDEEEEFFFF
        DDDDDDDGGGGHHHH"
pan1 <- volGene + volMirna + free(pRid) + free(pNet) + cor1 + cor2 + cor3 + cor4 +
    plot_layout(design = lay, heights = c(1, 1, 1.1, 1.1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

## add the custom legend
pan1 <- ggarrange(pan1, lg, nrow = 2, heights = c(20, 1))

## save the panel
ggsave("figures/DCM_de.pdf", pan1, width = 180, height = 160, units = "mm",
       scale = 1.5, device = cairo_pdf)

