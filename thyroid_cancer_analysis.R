## Integrative miRNA-mRNA analysis of thyroid cancer with MIRit
## ============================================================

## load MIRit
library(MIRit)

## define results folder
resFolder <- "results"
dir.create(resFolder)

## set seed for reproducible results
set.seed(1234)

## load data
data("mirnaCounts")
data("geneCounts")

## define samples metadata
meta <- data.frame(primary = colnames(mirnaCounts),
                   mirnaCol = colnames(mirnaCounts),
                   geneCol = colnames(geneCounts),
                   disease = c(rep("PTC", 8), rep("NTH", 8)),
                   patient = c(rep(paste("Sample_", seq(8), sep = ""), 2)))

## create a MirnaExperiment object
experiment <- MirnaExperiment(mirnaExpr = mirnaCounts,
                              geneExpr = geneCounts,
                              samplesMetadata = meta)

## perform miRNA and gene differential expression
experiment <- performGeneDE(experiment,
                            group = "disease",
                            contrast = "PTC-NTH",
                            design = ~ disease + patient,
                            method = "edgeR",
                            logFC = log2(2))
experiment <- performMirnaDE(experiment,
                             group = "disease",
                             contrast = "PTC-NTH",
                             design = ~ disease + patient,
                             method = "edgeR",
                             logFC = log2(2))

## perform gene functional enrichment
ora <- enrichGenes(experiment, method = "ORA", database = "GO")
gsea <- enrichGenes(experiment, method = "GSEA", database = "KEGG")

## retrieve miRNA target genes
experiment <- getTargets(experiment, score = "High")

## integrate miRNA and gene expression levels
experiment <- mirnaIntegration(experiment)

## perform ORA of integrated genes
oraInt <- enrichTargets(experiment, database = "DO")

## create miRNA-augmented pathways
pathways <- preparePathways(experiment, database = "KEGG")

## perform a topology-aware integrative pathway analysis (TAIPA)
taipa <- topologicalAnalysis(experiment,
                             pathways = pathways,
                             pCutoff = 0.05,
                             pAdjustment = "max-T",
                             nPerm = 10000)

## save resulting objects
saveRDS(experiment, paste(resFolder, "experiment.RDS", sep = "/"))
saveRDS(taipa, paste(resFolder, "taipa.RDS", sep = "/"))


## Supplementary tables generation
## ===============================

## extract and save differential expression results
write.csv(mirnaDE(experiment),
          file = paste(resFolder, "mirna_de.csv", sep = "/"),
          row.names = FALSE)
write.csv(geneDE(experiment),
          file = paste(resFolder, "gene_de.csv", sep = "/"),
          row.names = FALSE)

## extract and save miRNA-mRNA correlation analysis
write.csv(integration(experiment),
          file = paste(resFolder, "integration.csv", sep = "/"),
          row.names = FALSE)


## Figure generation
## =================

## load required libraries
library(ggpubr)

## define plots folder
pDir <- "plots"
dir.create(pDir)

## differential expression and functional enrichment panel
panel1.row1 <- ggarrange(
  plotDimensions(experiment, "genes", condition = "disease",
                 legend = "none", title = "Genes"),
  plotDimensions(experiment, "microRNA", condition = "disease",
                 legend = "none", title = "miRNAs"),
  plotVolcano(experiment, "genes", title = "DEGs"),
  plotVolcano(experiment, "microRNA", title = "DE-miRNAs"),
  nrow = 1, labels = "auto"
)
panel1.row2 <- ggarrange(
  plotDE(experiment, "TG", nameAsTitle = TRUE),
  plotDE(experiment, "TPO", nameAsTitle = TRUE),
  plotDE(experiment, "DIO2", nameAsTitle = TRUE),
  plotDE(experiment, "PAX8", nameAsTitle = TRUE),
  nrow = 1, labels = c("e", "f", "g", "h"), common.legend = TRUE
)
panel1.row3 <- ggarrange(
  enrichmentDotplot(ora$downregulated, title = "Depleted terms"),
  gseaPlot(gsea, "Thyroid hormone synthesis", rankingMetric = TRUE),
  nrow = 1, labels = c("i", "j"), widths = c(2, 1)
)
panel1 <- ggarrange(
  ggarrange(panel1.row1, panel1.row2, nrow = 2, align = "hv"),
  panel1.row3,
  nrow = 2, heights = c(1, 1)
)
ggsave(paste(pDir, "panel1.pdf", sep = "/"), panel1,
       width = 170, height = 178, units = "mm", scale = 1.5)

## integration panel
panel2.row1 <- ggarrange(
  plotCorrelation(experiment, "hsa-miR-146b-5p", "DIO2"),
  plotCorrelation(experiment, "hsa-miR-146b-5p", "PAX8"),
  plotCorrelation(experiment, "hsa-miR-146b-3p", "PAX8"),
  nrow = 1, labels = "auto", common.legend = TRUE, legend = "bottom"
)
panel2.row2 <- ggarrange(
  enrichmentDotplot(oraInt$downregulated, title = "Depleted diseases",
                    showTerms = 8, showTermsParam = "padj"),
  integrationDotplot(taipa, title = "Integrated pathways"),
  nrow = 1, labels = c("d", "e"), widths = c(1.2, 1)
)
panel2 <- ggarrange(
  panel2.row1, panel2.row2,
  nrow = 2, heights = c(1, 1)
)
ggsave(paste(pDir, "panel2.pdf", sep = "/"), panel2,
       width = 170, height = 128, units = "mm", scale = 1.5)

## create the network plot
pdf(paste(pDir, "network.pdf", sep = "/"), width = 7, height = 4)
visualizeNetwork(taipa, pathway = "Thyroid hormone synthesis")
dev.off()
