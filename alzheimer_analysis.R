## Integrative miRNA-mRNA analysis of Alzheimer's disease with MIRit
## =================================================================

## load MIRit
library(MIRit)

## define results folder
resFolder <- "results"
dir.create(resFolder)

## set seed for reproducible results
set.seed(1234)

## read miRNA count matrix
mir <- read.table("data/GSE63501/counts.csv",
                  header = TRUE, sep = ",", row.names = "X")

## read miRNA sample metadata
mirMeta <- read.csv("data/GSE63501/metadata.csv",
                    row.names = "X")

## read mRNA expression matrix
mes <- read.csv("data/GSE150696/normalized_expression.csv",
                row.names = "X")

## read mRNA sample metadata
mesMeta <- read.csv("data/GSE150696/metadata.csv",
                    row.names = "X")
mesMeta$sample <- gsub("/", ".", mesMeta$sample)
mesMeta$condition[mesMeta$condition == "CTRL"] <- "Control"

## create required metadata table
mirMeta$primary <- mirMeta$sample
mesMeta$primary <- mesMeta$sample
colnames(mirMeta)[colnames(mirMeta) == "sample"] <- "mirnaCol"
colnames(mirMeta)[colnames(mirMeta) == "classification"] <- "condition"
colnames(mesMeta)[colnames(mesMeta) == "sample"] <- "geneCol"
meta <- merge(mirMeta, mesMeta, by = c("primary", "condition"),
              all = TRUE)

## create MirnaExperiment object
experiment <- MirnaExperiment(mirnaExpr = mir,
                              geneExpr = mes,
                              samplesMetadata = meta,
                              pairedSamples = FALSE)

## perform miRNA and gene differential expression
experiment <- performMirnaDE(experiment,
                             group = "condition",
                             contrast = "AD-Control",
                             design = ~0 + condition,
                             method = "voom",
                             pCutoff = 0.05,
                             pAdjustment = "none")
experiment <- performGeneDE(experiment,
                            group = "condition",
                            contrast = "AD-Control",
                            design = ~ 0 + condition + sex + age,
                            method = "limma",
                            logFC = log2(1.5),
                            useArrayWeights = TRUE)

## perform disease-SNPs association
searchDisease("alzheimer")
snp <- findMirnaSNPs(experiment, "Alzheimer disease")

## retrieve miRNA target genes
experiment <- getTargets(experiment, score = "Medium")

## integrate miRNA and gene expression levels
experiment <- mirnaIntegration(experiment, associationMethod = "boschloo")

## create miRNA-augmented pathways
paths <- preparePathways(experiment, database = "WikiPathways")

## perform a topology-aware integrative pathway analysis (TAIPA)
taipa <- topologicalAnalysis(experiment,
                             paths,
                             pAdjustment = "BH",
                             nPerm = 20000)

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
library(gridGraphics)

## define plots folder
pDir <- "plots"
dir.create(pDir)

## alzheimer's disease panel
panel1.row1 <- ggarrange(
  plotDimensions(experiment, "genes", condition = "condition", legend = "left",
                 colorScale = c(AD = scales::hue_pal()(2)[2],
                                Control = scales::hue_pal()(2)[1]),
                 title = "Genes"),
  plotDimensions(experiment, "microRNA", condition = "condition", legend = "none",
                 colorScale = c(AD = scales::hue_pal()(2)[2],
                                Control = scales::hue_pal()(2)[1]),
                 title = "miRNAs"),
  plotVolcano(experiment, "genes", title = "DEGs"),
  plotVolcano(experiment, "microRNA", title = "DE-miRNAs"),
  nrow = 1, common.legend = TRUE, legend = "bottom", labels = "auto"
)
mirVariantPlot("rs2632516", snp, showContext = TRUE, showSequence = FALSE,
               cex = 1, cex.group = 0.8,
               from = 58329000, to = 58333000)
trackPlot <- grid.grab()
dev.off()
panel1.row2 <- ggarrange(
  trackPlot,
  integrationDotplot(taipa, title = "Integrated pathways"),
  nrow = 1, labels = c("e", "f"), widths = c(1, 2.5)
)
panel1 <- ggarrange(
  panel1.row1,
  panel1.row2,
  nrow = 2, heights = c(1, 1.4)
)
ggsave(paste(pDir, "panel3.pdf", sep = "/"), panel1,
       width = 170, height = 123, units = "mm", scale = 1.5)
