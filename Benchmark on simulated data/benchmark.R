## BENCHMARK OF DIFFERENT STATISTICAL TESTS FOR MIRNA INTEGRATIVE ANALYSES
## =======================================================================
##
## This simple script evaluates the performance of simple and partial Pearson
## and Spearman correlation analysis for the identification of influential
## miRNA-target relationships. In addition, the performance of Fisher's exact
## test, with or without Lancaster's midP correction, Boschloo's exact test,
## and rotation gene-set tests in identifying relevant active miRNAs is also
## evaluated through multiple performance metrics, including AUC, AUPRC,
## Precision, and F1.
##
## --------------------------------------------------

## load required libraries
library(limma)
library(ROCR)
library(caret)
library(ppcor)
library(exact2x2)
library(BiocParallel)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

## load helper functions
source("helper_functions.R")

## set the backend for parallel computation
BP <- MulticoreParam(workers = 28)


## RUN THE BENCHMARK
## -----------------

## set the seed for reproducibility
set.seed(1234)

## perform the simulation n times
nPerm <- 500
res <- bplapply(1:nPerm, function(p) {
    
    ## simulate a miRNA-mRNA experiment
    sim <- simulateExperiment()
    
    ## identify differentially expressed miRNAs and genes
    resMiR <- runLimma(sim$miRNA, sim$group)
    resGene <- runLimma(sim$gene, sim$group)
    upMir <- rownames(resMiR)[resMiR$adj.P.Val < 0.05 & resMiR$logFC > 0]
    downMir <- rownames(resMiR)[resMiR$adj.P.Val < 0.05 & resMiR$logFC < 0]
    upGene <- rownames(resGene)[resGene$adj.P.Val < 0.05 & resGene$logFC > 0]
    downGene <- rownames(resGene)[resGene$adj.P.Val < 0.05 & resGene$logFC < 0]
    
    ## identify all miRNA-target pairs
    mirnas <- names(sim$targets)
    allPairs <- do.call(rbind, lapply(seq_along(sim$targets), function(i) {
        if (length(sim$targets[[i]]) > 0) {
            data.frame(
                miRNA = rep(mirnas[i], length(sim$targets[[i]])),
                gene  = sim$targets[[i]],
                stringsAsFactors = FALSE
            )
        } else {
            NULL
        }
    }))
    
    ## filter meaningful pairs
    allPairs <- allPairs[(allPairs$miRNA %in% upMir & allPairs$gene %in% downGene) |
                             (allPairs$miRNA %in% downMir & allPairs$gene %in% upGene), ]
    
    ## define TRUE and FALSE pairs
    allPairs$label <- paste(allPairs$miRNA, allPairs$gene) %in%
        paste(sim$true_pairs$miRNA, sim$true_pairs$gene)
    
    ## perform simple correlation analyses
    simPearson <- simpleCorrelation(sim, allPairs, method = "pearson")
    simSpearman <- simpleCorrelation(sim, allPairs, method = "spearman")
    
    ## perform partial correlation analyses
    parPearson <- partialCorrelation(sim, allPairs, method = "pearson")
    parSpearman <- partialCorrelation(sim, allPairs, method = "spearman")
    
    ## estimate the performance of the different correlation approaches
    corMethods <- list(simPearson, simSpearman,
                       parPearson, parSpearman)
    pfTable <- as.data.frame(t(sapply(corMethods, calculatePerformance)))
    pfTable$method <- c("Simple Pearson", "Simple Spearman",
                        "Partial Pearson", "Partial Spearman")
    
    ## estimate the AUC and AUCPR of the different correlation approaches
    auMet <- as.data.frame(t(sapply(corMethods, calculateCurves)))
    colnames(auMet) <- c("AUC", "AUPRC")
    auMet$method <- pfTable$method
    
    ## perform the integration using unpaired methods
    pFisher <- fisherTest(upMir, downMir, upGene, downGene, sim)
    pMidp <- fisherMidp(upMir, downMir, upGene, downGene, sim)
    pBoschloo <- boschlooTest(upMir, downMir, upGene, downGene, sim)
    pFry <- fryTest(upMir, downMir, upGene, downGene, sim)
    pCamera <- cameraTest(upMir, downMir, upGene, downGene, sim)
    
    ## evaluate the performance of unpaired approaches in finding relevant miRNAs
    unpMethods <- list(pFisher, pMidp, pFry, pCamera, pBoschloo)
    unpMirs <- as.data.frame(t(sapply(unpMethods, function(x) {
        colnames(x)[4] <- "label"
        calculatePerformance(x)
    })))
    auMirs <- as.data.frame(t(sapply(unpMethods, function(x) {
        colnames(x)[colnames(x) == "active"] <- "label"
        x$cor <- log10(x$pval)
        calculateCurves(x)
    })))
    colnames(auMirs) <- c("AUC", "AUPRC")
    auMirs$method <- unpMirs$method <- c("Fisher's exact test",
                                         "Fisher's test with midP adjustment",
                                         "Fry test", "CAMERA test",
                                         "Boschloo's exact test")
    
    ## combine performance metrics
    intTable <- merge(pfTable, auMet, by = "method", all.x = TRUE)
    mirTable <- merge(unpMirs, auMirs, by = "method", all.x = TRUE)
    
    ## mark groups
    intTable$type <- "interactions"
    mirTable$type <- "miRNAs"
    
    ## combine correlation + unpaired methods
    resTab <- rbind(intTable, mirTable)
    
    ## add iteration index
    resTab$iteration <- p
    
    ## prepare a list with the resulting identifications
    names(corMethods) <- pfTable$method
    names(unpMethods) <- unpMirs$method
    results <- list(interactions = corMethods,
                    miRNAs = unpMethods)
    
    ## return results and performance
    return(list(performance = resTab,
                results = results))
    
}, BPPARAM = BP)

## create the final performance table
pfList <- lapply(res, function(x) x$performance)
benchmark <- do.call(rbind, pfList)

## save the benchmark table
write.csv(benchmark, "output/benchmark_performance.csv")

## create a list with predictions for each iteration
resList <- lapply(res, function(x) x$results)

## save raw simulation results
saveRDS(resList, "output/benchmark.RDS")

## create a performance summary
pfSummary <- do.call(rbind, by(benchmark[, 1:22],
                               benchmark$method,
                               function(x) {
                                   
        ## select numeric data
        numData <- x[, sapply(x, is.numeric), drop = FALSE]
        
        ## compute mean and sd
        means <- colMeans(numData, na.rm = TRUE)
        sds <- sapply(numData, sd, na.rm = TRUE)
        
        ## create mean ± sd strings
        metrics <- paste0(round(means, 3), " ± ", round(sds, 3))
        
        ## return as data.frame with method
        data.frame(method = unique(x$method), t(metrics), row.names = NULL)
        
    }
))
colnames(pfSummary) <- colnames(benchmark)[1:22]
pfSummary$type <- "interactions"
pfSummary$type[1:5] <- "miRNAs"
write.csv(pfSummary, "output/performance_summary.csv")


## COMPARE F1 AT VARYING SAMPLE SIZES
## ----------------------------------

## set the seed for reproducibility
set.seed(1234)

## define sample sizes
sizes <- c(10, 15, 20, 25, 30)

## monitor performance at varying sample sizes
nRep <- 100
f1SS <- bplapply(1:nRep, function(p) {
    
    ## define a data.frame for storing F1 scores
    f1forPerm <- matrix(NA, nrow = 5, ncol = 4)
    colnames(f1forPerm) <- c("Simple Pearson", "Simple Spearman",
                             "Partial Pearson", "Partial Spearman")
    
    ## loop through various sample sizes
    for (i in 1:length(sizes)) {
        
        ## simulate a miRNA-mRNA experiment
        ss <- sizes[i]
        sim <- simulateExperiment(nA = ss, nB = ss)
        
        ## identify differentially expressed miRNAs and genes
        resMiR <- runLimma(sim$miRNA, sim$group)
        resGene <- runLimma(sim$gene, sim$group)
        upMir <- rownames(resMiR)[resMiR$adj.P.Val < 0.05 & resMiR$logFC > 0]
        downMir <- rownames(resMiR)[resMiR$adj.P.Val < 0.05 & resMiR$logFC < 0]
        upGene <- rownames(resGene)[resGene$adj.P.Val < 0.05 & resGene$logFC > 0]
        downGene <- rownames(resGene)[resGene$adj.P.Val < 0.05 & resGene$logFC < 0]
        
        ## identify all miRNA-target pairs
        mirnas <- names(sim$targets)
        allPairs <- do.call(rbind, lapply(seq_along(sim$targets), function(i) {
            if (length(sim$targets[[i]]) > 0) {
                data.frame(
                    miRNA = rep(mirnas[i], length(sim$targets[[i]])),
                    gene  = sim$targets[[i]],
                    stringsAsFactors = FALSE
                )
            } else {
                NULL
            }
        }))
        
        ## filter meaningful pairs
        allPairs <- allPairs[(allPairs$miRNA %in% upMir & allPairs$gene %in% downGene) |
                                 (allPairs$miRNA %in% downMir & allPairs$gene %in% upGene), ]
        
        ## define TRUE and FALSE pairs
        allPairs$label <- paste(allPairs$miRNA, allPairs$gene) %in%
            paste(sim$true_pairs$miRNA, sim$true_pairs$gene)
        
        ## perform simple correlation analyses
        simPearson <- simpleCorrelation(sim, allPairs, method = "pearson")
        simSpearman <- simpleCorrelation(sim, allPairs, method = "spearman")
        
        ## perform partial correlation analyses
        parPearson <- partialCorrelation(sim, allPairs, method = "pearson")
        parSpearman <- partialCorrelation(sim, allPairs, method = "spearman")
        
        ## estimate the performance of the different correlation approaches
        corMethods <- list(simPearson, simSpearman,
                           parPearson, parSpearman)
        pfTable <- as.data.frame(t(sapply(corMethods, calculatePerformance)))
        pfTable$method <- c("Simple Pearson", "Simple Spearman",
                            "Partial Pearson", "Partial Spearman")
        
        ## add F1 scores
        f1 <- pfTable$F1
        f1forPerm[i, ] <- f1
        
    }
    
    ## convert to a suitable data.frame
    resIter <- as.data.frame(f1forPerm)
    resIter$sample.size <- sizes
    resIter$iteration <- p
    
    ## return results
    return(resIter)
    
}, BPPARAM = BP)

## create a single table
f1Res <- do.call(rbind, f1SS)

## save results
write.csv(f1Res, "output/F1_varying_sizes.csv")

