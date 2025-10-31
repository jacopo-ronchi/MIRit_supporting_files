## helper function for simulating a paired miRNA-mRNA experiment where
## influential miRNA-target pairs show differential expression in opposite
## direction and also present anti-correlated expression values
simulateExperiment <- function(G = 20000,
                               M = 2500,
                               nA = 50,
                               nB = 50,
                               mir_mu = 2,
                               gene_mu = 3,
                               mir_sd = 0.4,
                               gene_sd = 0.4,
                               samp_sd = 0.2,
                               n_DE_miRNAs = 50,
                               frac_active = 0.3,
                               mean_targets_per_miRNA = 300,
                               prop_targets_affected = 0.3,
                               n_DE_genes = 1000,
                               group_effect = 1,
                               sigma_effect = 0.4,
                               hill_range = c(2, 4),
                               K_range = c(4, 12),
                               eff_range = c(0.2, 0.5),
                               alpha = 4,
                               noise = 0.1,
                               verbose = FALSE) {
    
    ## define the total number of samples and features
    N <- nA + nB
    group <- c(rep(0, nA), rep(1, nB))
    sample_names <- paste0("S", seq_len(N))
    mir_names <- paste0("miR_", seq_len(M))
    gene_names <- paste0("gene_", seq_len(G))
    
    ## simulate log-mean expression for each gene and miRNA
    mirMean  <- rnorm(M, mean = mir_mu, sd = mir_sd)
    geneMean <- rnorm(G, mean = gene_mu, sd = gene_sd)
    
    ## define baseline expression on the linear scale
    mirExpr <- t(sapply(mirMean, function(j) rlnorm(N, meanlog = j, sdlog = samp_sd)))
    geneExpr <- t(sapply(geneMean, function(j) rlnorm(N, meanlog = j, sdlog = samp_sd)))
    rownames(mirExpr) <- mir_names
    rownames(geneExpr) <- gene_names
    colnames(mirExpr) <- colnames(mirExpr) <- sample_names
    
    ## define the probabilities of being targeted by miRNAs
    gene_suscept <- rbeta(G, 2, 100)
    gene_prob <- gene_suscept / sum(gene_suscept)
    
    ## randomly define miRNA targets
    targets <- setNames(vector("list", M), mir_names)
    targets_count <- pmax(1, rnbinom(M, size = 5, mu = mean_targets_per_miRNA))
    for (i in seq_len(M)) {
        targets[[i]] <- sample(gene_names, size = min(targets_count[i], G),
                               replace = FALSE, prob = gene_prob)
    }
    
    ## choose DE-miRNAs and the subset of "active" ones
    deIdx <- sample(seq_len(M), n_DE_miRNAs)
    mirSigns <- sample(c(-1, 1), n_DE_miRNAs, replace = TRUE)
    names(mirSigns) <- mir_names[deIdx]
    nActive <- ceiling(n_DE_miRNAs * frac_active)
    activeIdx <- sample(deIdx, nActive)
    
    ## add group effect to DE miRNAs
    log2FC <- mirSigns * rnorm(n_DE_miRNAs, mean = group_effect, sd = sigma_effect)
    FC <- 2^log2FC
    for (j in seq_along(deIdx)) {
        i <- deIdx[j]
        mirExpr[i, group == 1] <- mirExpr[i, group == 1] * FC[j]
    }
    
    ## define the list of active miRNA-target pairs
    truePairs <- data.frame(miRNA = character(),
                            gene = character(),
                            stringsAsFactors = FALSE)
    affectedList <- list()
    for (i in activeIdx) {
        mir <- mir_names[i]
        tgs <- targets[[i]]
        nAff <- max(1, round(length(tgs) * prop_targets_affected *
                                 runif(1, 0.7, 1.3)))
        affected <- sample(tgs, nAff)
        affectedList[[mir]] <- affected
        truePairs <- rbind(truePairs, data.frame(miRNA = mir,
                                                 gene = affected,
                                                 stringsAsFactors = FALSE))
    }
    
    ## calculate perturbed gene expression
    geneExpr <- computePerturbation(mirExpr, geneExpr, affectedList,
                                    K_range, hill_range, eff_range, alpha)
    
    ## add extra independent DEGs
    remaining <- setdiff(gene_names, unique(truePairs$gene))
    degs <- sample(remaining, n_DE_genes)
    signs <- sample(c(-1, 1), length(degs), replace = TRUE)
    effects <- signs * rnorm(length(degs), mean = group_effect, sd = sigma_effect)
    effects <- 2^effects
    geneExpr[degs, group == 1] <- geneExpr[degs, group == 1] * effects
    
    ## add a small random noise
    mirExpr <- mirExpr * matrix(rlnorm(length(mirExpr), 0, noise), nrow = M)
    geneExpr <- geneExpr * matrix(rlnorm(length(geneExpr), 0, noise), nrow = G)
    
    ## print diagnostics
    if (verbose) {
        mdNonTg <- median(geneExpr[!rownames(geneExpr) %in% truePairs$gene, ])
        mdTg <- median(geneExpr[rownames(geneExpr) %in% truePairs$gene, ])
        sp <- sapply(1:nrow(truePairs), function(i) {
            cor(mirExpr[truePairs$miRNA[i], ],
                geneExpr[truePairs$gene[i], ],
                method = "spearman")
        })
        pr <- sapply(1:nrow(truePairs), function(i) {
            cor(mirExpr[truePairs$miRNA[i], ],
                geneExpr[truePairs$gene[i], ],
                method = "pearson")
        })
        cat("Average Pearson correlation:", round(mean(pr), 3), "\n",
            "Range: ", paste(round(range(pr), 2), collapse = ":"), "\n",
            "Average Spearman correlation:", round(mean(sp), 3), "\n",
            "Range: ", paste(round(range(sp), 2), collapse = ":"), "\n",
            "Median miRNA expression:", round(median(mirExpr), 2), "\n",
            "Median expression of non-target genes:", round(mdNonTg, 2), "\n",
            "Median expression of target genes:", round(mdTg, 2))
    }
    
    ## return results
    res <- list(miRNA = mirExpr,
                gene = geneExpr,
                group = group,
                targets = targets,
                DE_miRNAs = mir_names[deIdx],
                active_miRNAs = mir_names[activeIdx],
                DE_genes = degs,
                true_pairs = truePairs,
                affected_genes = unique(truePairs$gene))
    return(res)
    
}



## helper function for calculating gene expression after miRNA regulation
computePerturbation <- function(mi_mat,
                                gene_baseline,
                                miRNA_targets,
                                hill_range,
                                K_range,
                                eff_range,
                                alpha) {
    
    ## define the number of features and samples
    n_genes <- nrow(gene_baseline)
    n_samples <- ncol(gene_baseline)
    
    ## create a list with regulators for each gene
    gene2regulators <- vector("list", n_genes)
    names(gene2regulators) <- rownames(gene_baseline)
    for (m in names(miRNA_targets)) {
        targets <- miRNA_targets[[m]]
        gene2regulators[targets] <- Map(c, gene2regulators[targets], list(m))
    }
    
    ## define a matrix to store perturbation factors
    totP <- matrix(1, nrow = nrow(gene_baseline), ncol = n_samples,
                   dimnames = list(rownames(gene_baseline),
                                   colnames(gene_baseline)))
    
    ## loop over each gene
    for (g in rownames(gene_baseline)) {
        
        ## extract all the miRNAs regulating this gene
        regulators <- gene2regulators[[g]]
        if (length(regulators) == 0L) next
        
        ## define a vector to store miRNA contributions
        s <- numeric(n_samples)
        
        ## iterate over each regulator
        for (m in regulators) {
            
            ## extract miRNA expression
            mi_exp <- mi_mat[m, ]
            
            ## simulate repression efficiency, K and Hill parameters
            h <- runif(1, hill_range[1], hill_range[2])
            K <- runif(1, K_range[1], K_range[2])
            eff <- runif(1, eff_range[1], eff_range[2])
            
            ## calculate miRNA repression contribution
            rf <- mi_exp^h / (K^h + mi_exp^h)
            s <- s + eff * rf
            
        }
        
        ## calculate overall perturbation
        totP[g, ] <- 1 - (s / (1 + s))
        
    }
    
    ## return the final gene expression
    geneNorm <- log2(gene_baseline + 1e-6)
    totNorm <- log2(totP)
    geNorm <- geneNorm + alpha * totNorm
    ge <- 2^geNorm - 1e-6
    return(ge)
    
}



## very simple function for detecting DE features using the limma pipeline
runLimma <- function(expr, group) {
    expr <- log2(expr + 1e-6)
    design <- model.matrix(~ group)
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    top <- topTable(fit, coef = 2, number = nrow(expr))
    top
}



## helper function for performing a simple miRNA-mRNA correlation analysis
simpleCorrelation <- function(sim, pairs, method = "pearson") {
    
    ## calculate the correlation
    cors <- t(mapply(function(m, g) {
        ct <- cor.test(sim$miRNA[m,], sim$gene[g,], method = method)
        c(ct$estimate, ct$p.value)
    }, pairs$miRNA, pairs$gene))
    colnames(cors) <- c("cor", "pvalue")
    
    ## format the results
    pairs <- cbind(pairs, cors)
    
    ## adjust p-values
    pairs$FDR <- p.adjust(pairs$pvalue, method = "BH")
    
    ## define significant pairs
    pairs$sig <- FALSE
    pairs$sig[pairs$FDR < 0.05] <- TRUE
    
    ## return results
    return(pairs)
    
}



## helper function for performing a partial miRNA-mRNA correlation analysis
partialCorrelation <- function(sim, pairs, method = "pearson") {
    
    ## calculate the partial correlation
    cors <- t(mapply(function(m, g) {
        pc <- pcor.test(sim$miRNA[m, ], sim$gene[g, ],
                        sim$group, method = method)
        c(pc$estimate, pc$p.value)
    }, pairs$miRNA, pairs$gene))
    colnames(cors) <- c("cor", "pvalue")
    
    ## format the results
    pairs <- cbind(pairs, cors)
    
    ## adjust p-values
    pairs$FDR <- p.adjust(pairs$pvalue, method = "BH")
    
    ## define significant pairs
    pairs$sig <- FALSE
    pairs$sig[pairs$FDR < 0.05] <- TRUE
    
    ## return results
    return(pairs)
    
}



## helper function for identifying influential miRNAs using Fisher's exact test
fisherTest <- function(upMir, downMir, upGene, downGene, sim) {
    
    ## calculate p-values using Fisher's exact test
    fPU <- sapply(upMir, function(m){
        tgs <- sim$targets[[m]]
        tbl <- matrix(c(
            sum(tgs %in% downGene),
            sum(!tgs %in% downGene),
            sum(setdiff(rownames(sim$gene), tgs) %in% downGene),
            sum(!setdiff(rownames(sim$gene), tgs) %in% downGene)
        ), nrow = 2)
        fisher.test(tbl, alternative = "greater")$p.value
    })
    fPD <- sapply(downMir, function(m){
        tgs <- sim$targets[[m]]
        tbl <- matrix(c(
            sum(tgs %in% upGene),
            sum(!tgs %in% upGene),
            sum(setdiff(rownames(sim$gene), tgs) %in% upGene),
            sum(!setdiff(rownames(sim$gene), tgs) %in% upGene)
        ), nrow = 2)
        fisher.test(tbl, alternative = "greater")$p.value
    })
    
    ## create a table with the results
    res <- data.frame(miRNA = c(upMir, downMir), pval = c(fPU, fPD))
    res$FDR <- p.adjust(res$pval)
    res$active <- res$miRNA %in% sim$active_miRNAs
    
    ## define significant pairs
    res$sig <- FALSE
    res$sig[res$FDR < 0.05] <- TRUE
    
    ## return results
    return(res)
    
}



## helper function for identifying influential miRNAs using Fisher's exact test
## with Lancaster's midP correction
fisherMidp <- function(upMir, downMir, upGene, downGene, sim) {
    
    ## calculate p-values using Fisher's exact test (upregulated miRNAs)
    fPU <- sapply(upMir, function(m){
        tgs <- sim$targets[[m]]
        tbl <- matrix(c(
            sum(tgs %in% downGene),
            sum(!tgs %in% downGene),
            sum(setdiff(rownames(sim$gene), tgs) %in% downGene),
            sum(!setdiff(rownames(sim$gene), tgs) %in% downGene)
        ), nrow = 2)
        fisher.exact(tbl, alternative = "greater", midp = TRUE)$p.value
    })
    fPD <- sapply(downMir, function(m){
        tgs <- sim$targets[[m]]
        tbl <- matrix(c(
            sum(tgs %in% upGene),
            sum(!tgs %in% upGene),
            sum(setdiff(rownames(sim$gene), tgs) %in% upGene),
            sum(!setdiff(rownames(sim$gene), tgs) %in% upGene)
        ), nrow = 2)
        fisher.exact(tbl, alternative = "greater", midp = TRUE)$p.value
    })
    
    ## create a table with the results
    res <- data.frame(miRNA = c(upMir, downMir), pval = c(fPU, fPD))
    res$FDR <- p.adjust(res$pval)
    res$active <- res$miRNA %in% sim$active_miRNAs
    
    ## define significant pairs
    res$sig <- FALSE
    res$sig[res$FDR < 0.05] <- TRUE
    
    ## return results
    return(res)
    
}



## helper function for identifying influential miRNAs using Boschloo's test
boschlooTest <- function(upMir, downMir, upGene, downGene, sim) {
    
    ## calculate p-values using Boschloo's exact test
    bPU <- sapply(upMir, function(m){
        tgs <- sim$targets[[m]]
        tbl <- matrix(c(
            sum(tgs %in% downGene),
            sum(!tgs %in% downGene),
            sum(setdiff(rownames(sim$gene), tgs) %in% downGene),
            sum(!setdiff(rownames(sim$gene), tgs) %in% downGene)
        ), nrow = 2)
        MIRit:::boshloo.test(tbl[2:1, ], 100)
    })
    bPD <- sapply(downMir, function(m){
        tgs <- sim$targets[[m]]
        tbl <- matrix(c(
            sum(tgs %in% upGene),
            sum(!tgs %in% upGene),
            sum(setdiff(rownames(sim$gene), tgs) %in% upGene),
            sum(!setdiff(rownames(sim$gene), tgs) %in% upGene)
        ), nrow = 2)
        MIRit:::boshloo.test(tbl[2:1, ], 100)
    })
    
    ## create a table with the results
    res <- data.frame(miRNA = c(upMir, downMir), pval = c(bPU, bPD))
    res$FDR <- p.adjust(res$pval)
    res$active <- res$miRNA %in% sim$active_miRNAs
    
    ## define significant pairs
    res$sig <- FALSE
    res$sig[res$FDR < 0.05] <- TRUE
    
    ## return results
    return(res)
    
}



## helper function for identifying influential miRNAs using the 'fry' method
fryTest <- function(upMir, downMir, upGene, downGene, sim) {
    
    ## perform the 'fry' method
    gene <- log2(sim$gene + 1e-6)
    fr <- fry(gene, index = sim$targets[c(upMir, downMir)],
              design = model.matrix(~ sim$group))
    
    ## create a table with the results
    fr$miRNA <- rownames(fr)
    fr$active <- fr$miRNA %in% sim$active_miRNAs
    fr <- fr[, c(8, 7, 2, 3, 4)]
    colnames(fr)[4] <- "pval"
    
    ## define significant pairs
    fr$sig <- FALSE
    fr$sig[fr$FDR < 0.05 &
               (fr$miRNA %in% upMir & fr$Direction == "Down") |
               (fr$miRNA %in% downMir & fr$Direction == "Up")] <- TRUE
    
    ## reorder the data.frame
    fr <- fr[, c(2, 4, 5, 1, 6)]
    
    ## return results
    return(fr)
    
}



## helper function for identifying influential miRNAs using the 'CAMERA' method
cameraTest <- function(upMir, downMir, upGene, downGene, sim) {
    
    ## perform the 'CAMERA' method
    gene <- log2(sim$gene + 1e-6)
    cm <- camera(gene, index = sim$targets[c(upMir, downMir)],
                 design = model.matrix(~ sim$group))
    
    ## create a table with the results
    cm$miRNA <- rownames(cm)
    cm$active <- cm$miRNA %in% sim$active_miRNAs
    colnames(cm)[3] <- "pval"
    
    ## define significant pairs
    cm$sig <- FALSE
    cm$sig[cm$FDR < 0.05 &
               (cm$miRNA %in% upMir & cm$Direction == "Down") |
               (cm$miRNA %in% downMir & cm$Direction == "Up")] <- TRUE
    
    ## reorder the data.frame
    cm <- cm[, c(5, 3, 4, 6, 7)]
    
    ## return results
    return(cm)
    
}



## helper function for estimating the performance of different tests
calculatePerformance <- function(pairs) {
    
    ## make sure that the levels are the same
    pairs$label <- factor(pairs$label, levels = c(FALSE, TRUE))
    pairs$sig   <- factor(pairs$sig, levels = c(FALSE, TRUE))
    
    ## estimate performance with the caret package
    cm <- confusionMatrix(data = pairs$sig,
                          reference = pairs$label,
                          positive = "TRUE")
    
    ## calculate the FDR
    fdr <- sum(pairs$sig == TRUE & pairs$label == FALSE) /
        sum(pairs$sig == TRUE)
    
    ## create a vector with the performance metrics
    met <- c(cm$overall, cm$byClass, FDR = fdr)
    
    ## return performance metrics
    return(met)
    
}



## helper function for estimating AUC and AUPRC
calculateCurves<- function(pairs) {
    
    ## calculate AUC and AUPRC
    pred <- prediction(-pairs$cor, pairs$label)
    auc <- as.numeric(performance(pred, measure = "auc")@y.values)
    auprc <- as.numeric(performance(pred, measure = "aucpr")@y.values)
    
    ## return performance metrics
    return(c(auc, auprc))
    
}



## helper function for computing precision-recall curves
getPrCurve <- function(df) {
    pred <- prediction(-df$cor, df$label)
    perf <- performance(pred, "prec", "rec")
    data.frame(
        recall = perf@x.values[[1]],
        precision = perf@y.values[[1]]
    )
}



## helper function for comparing paired performance metrics through
## Wilcoxon test followed by Holm correction
testMetrics <- function(data, metric, comparisons) {
    
    ## prepare the comparison
    data <- data[, c(metric, "method", "iteration")]
    wide <- reshape(data,
                    idvar = "iteration",
                    timevar = "method",
                    direction = "wide")
    colnames(wide)[2:(length(unique(data$method)) + 1)] <- unique(data$method)
    
    ## compare methods
    res <- as.data.frame(t(sapply(comparisons, function(cmp) {
        tst <- wilcox.test(wide[, cmp[1]], wide[, cmp[2]], paired = TRUE)
        c(paste(cmp, collapse = " vs "), tst$p.value)
    })))
    colnames(res) <- c("comparison", "P.Value")
    res$adj.P.Val <- p.adjust(res$P.Value, method = "holm")
    res$star <- sapply(res$adj.P.Val, mapSignificance)
    
    ## return results
    return(res)
    
}



## helper function fo rmapping significance to star
mapSignificance <- function(p) {
    
    ## map p-value to star
    if (p > 0.05) {
        return("ns")
    } else if (p <= 0.05 & p > 0.01) {
        return("*")
    } else if (p <= 0.01 & p > 0.001) {
        return("**")
    } else if (p <= 0.001 & p > 0.0001) {
        return("***")
    } else if (p <= 0.0001) {
        return("****")
    }
    
}

