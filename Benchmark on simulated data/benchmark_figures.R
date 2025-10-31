## GENERATION OF BENCHMARK FIGURES
## ===============================
##
## This script is used to generate the plots for reporting benchmarking results
## for the comparison of both correlation-based and association-based tests in
## finding influential miRNA-target interactions.
##
## --------------------------------------------------

## load required libraries
library(ROCR)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(patchwork)
library(ggsci)
library(tidyr)
library(dplyr)

## load helper functions
source("helper_functions.R")

## load benchmark results
resList <- readRDS("output/benchmark.RDS")
benchmark <- read.csv("output/benchmark_performance.csv")
benchmark$method <- gsub(" adjustment", "", benchmark$method)

## split results for interactions and miRNAs
intRes <- lapply(resList, function(x) x$interactions)
mirRes <- lapply(resList, function(x) x$miRNAs)

## define a colorscale for each unique method
cScale <- pal_npg("nrc")(length(unique(benchmark$method)))
names(cScale) <- unique(benchmark$method)


## METRICS FOR INTERACTIONS
## ------------------------

## restrict to interactions
intTab <- benchmark[benchmark$type == "interactions", ]

## create the plot for the F1 in finding influential pairs
cmp <- list(c("Partial Pearson", "Partial Spearman"))
f1Corr <- testMetrics(intTab, "F1", cmp)
f1Corr
f1CorPlot <- ggviolin(intTab, x = "method", y = "F1", fill = "method",
                      add = "median_iqr", add.params = list(color = "black")) +
    xlab(NULL) +
    geom_signif(comparisons = cmp, annotations = f1Corr$star, textsize = 6, margin_top = 0.2) +
    scale_fill_manual(values = cScale) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                       expand = expansion(mult = c(0, 0.15))) +
    theme_pubr(x.text.angle = 45)

## create the plot for the AUPRC in finding influential pairs
cmp <- list(c("Partial Pearson", "Partial Spearman"),
            c("Simple Pearson", "Simple Spearman"))
auprcCorr <- testMetrics(intTab, "AUPRC", cmp)
auprcCorr
auprcCorPlot <- ggviolin(intTab, x = "method", y = "AUPRC", fill = "method",
                         add = "median_iqr", add.params = list(color = "black")) +
    xlab(NULL) +
    geom_signif(comparisons = cmp, annotations = auprcCorr$star, textsize = 6, margin_top = 0.2) +
    scale_fill_manual(values = cScale) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                       expand = expansion(mult = c(0, 0.15))) +
    theme_pubr(x.text.angle = 45)


## METRICS FOR INFLUENTIAL MIRNAS
## ------------------------------

## restrict to influential miRNAs
mirTab <- benchmark[benchmark$type == "miRNAs", ]

## create the plot for the F1 in finding influential miRNAs
f1UnpPlot <- ggviolin(mirTab, x = "method", y = "F1", fill = "method",
                      add = "median_iqr", add.params = list(color = "black")) +
    xlab(NULL) +
    scale_fill_manual(values = cScale) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_pubr(x.text.angle = 45)

## create the plot for the AUPRC in finding influential miRNAs
auprcUnpPlot <- ggviolin(mirTab, x = "method", y = "AUPRC", fill = "method",
                         add = "median_iqr", add.params = list(color = "black")) +
    xlab(NULL) +
    scale_fill_manual(values = cScale) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_pubr(x.text.angle = 45)


## AVERAGE PRECISION-RECALL CURVES
## -------------------------------

## aggregate results from all iterations
methods <- names(intRes[[1]])
allCurves <- list()
for (m in methods) {
    
    ## compute the curves
    crvs <- lapply(intRes, function(x) getPrCurve(x[[m]]))
    
    ## interpolate all curves on a common recall grid
    recallGrid <- seq(0, 1, length.out = 200)
    interp <- lapply(crvs, function(df) {
        approx(df$recall, df$precision, xout = recallGrid, rule = 2)$y
    })
    
    ## define the mean and sd
    mat <- do.call(rbind, interp)
    meanPrec <- colMeans(mat, na.rm = TRUE)
    sdPrec   <- apply(mat, 2, sd, na.rm = TRUE)
    se <- sdPrec / sqrt(length(intRes))
    lower <- meanPrec - 1.96 * se
    upper <- meanPrec + 1.96 * se
    
    ## add results for this method
    allCurves[[m]] <- data.frame(recall = recallGrid,
                                 mean_precision = meanPrec,
                                 lower = lower,
                                 upper = upper,
                                 Method = m)
    
}

## combine into a single data.frame for plotting
crvDf <- do.call(rbind, allCurves)

## create the plot
avgAucpr <- ggplot(crvDf, aes(x = recall, y = mean_precision,
                              color = Method, fill = Method)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    scale_color_manual(values = cScale) +
    scale_fill_manual(values = cScale) +
    labs(x = "Recall",
         y = "Precision") +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_pubr(border = TRUE)


## DENSITY PLOTS FOR SPEARMAN VS PEARSON
## -------------------------------------

## extract correlation coefficients
spVals <- lapply(intRes, function(x) x$`Simple Pearson`$cor)
ssVals <- lapply(intRes, function(x) x$`Simple Spearman`$cor)
ppVals <- lapply(intRes, function(x) x$`Partial Pearson`$cor)
psVals <- lapply(intRes, function(x) x$`Partial Spearman`$cor)

## define true and false interactions
intLab <- unlist(lapply(intRes, function(x) x$`Simple Pearson`$label))
lbl <- rep("Not Affected", length(intLab))
lbl[intLab] <- "Affected"

## create a data.frame for plotting
spVals <- data.frame(Coefficient = unlist(spVals),
                     Type = "Pearson",
                     Label = lbl,
                     Correlation = "Simple")
ssVals <- data.frame(Coefficient = unlist(ssVals),
                     Type = "Spearman",
                     Label = lbl,
                     Correlation = "Simple")
ppVals <- data.frame(Coefficient = unlist(ppVals),
                     Type = "Pearson",
                     Label = lbl,
                     Correlation = "Partial")
psVals <- data.frame(Coefficient = unlist(psVals),
                     Type = "Spearman",
                     Label = lbl,
                     Correlation = "Partial")
corDf <- rbind(spVals, ssVals, ppVals, psVals)

## create the density plot
densPlot <- ggplot(corDf, aes(x = Coefficient, group = Type, fill = Type)) +
    geom_density(adjust = 1.5, alpha = 0.4) +
    facet_grid(Correlation ~ Label) +
    ylab("Density") +
    theme_pubr(border = TRUE)


## CREATE THE PANEL FOR THE MAIN TEXT
## ----------------------------------

## create the panel
bmPan <- free(densPlot, type = "label") + f1CorPlot + auprcCorPlot +
    free(avgAucpr, type = "label") + f1UnpPlot + auprcUnpPlot &
    theme(legend.position = "none")

## create a ghost plot only for defining a specific legend
mtLeg <- cowplot::get_legend(avgAucpr)
dnLeg <- cowplot::get_legend(densPlot)
leg <- (wrap_plots(dnLeg) + wrap_plots(mtLeg)) + plot_layout(widths = c(1, 3))

## add the new legend
bmPan <- bmPan / leg +
    plot_layout(heights = c(1, 0.05))

## adjust the panel
bmPan <- bmPan +
    plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "", ""))) &
    theme(plot.tag = element_text(size = 14, face = "bold"))

## save the panel
ggsave("output/benchmark_figure.pdf", bmPan, width = 180, height = 150, units = "mm",
       scale = 1.5, device = cairo_pdf)


## CREATE THE SUPPLEMENTARY PANEL
## ------------------------------

## create the FDR plot
fdr <- ggviolin(benchmark, x = "method", y = "FDR", fill = "method",
                add = "median_iqr", add.params = list(color = "black")) +
    xlab(NULL) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    scale_fill_manual(values = cScale) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_pubr(x.text.angle = 45, legend = "none")

## load F1 results at various sample sizes
f1Res <- read.csv("output/F1_varying_sizes.csv")

## reformat the data.frame
colnames(f1Res)[2:6] <- c("Simple Pearson", "Simple Spearman",
                          "Partial Pearson", "Partial Spearman", "Sample Size")
f1Res <- f1Res |>
    select(-X) |>
    pivot_longer(
        cols = c("Simple Pearson", "Simple Spearman",
                 "Partial Pearson", "Partial Spearman"),
        names_to = "method",
        values_to = "F1")

## summarize data
f1Summary <- f1Res %>%
    group_by(`Sample Size`, method) %>%
    summarise(
        mean = mean(F1, na.rm = TRUE),
        se = sd(F1, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    )

## create the plot
f1overSize <- ggplot(f1Summary, aes(x = `Sample Size`, y = mean,
                                    color = method, group = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
    scale_color_manual(values = cScale) +
    labs(x = "Sample Size",
         y = "F1",
         color = "Method") +
    theme_pubr(legend = "bottom")

## create the panel
supPan <- ((fdr + guides(fill = "none")) + free(f1overSize, type = "label")) +
    plot_layout(guides = "collect", widths = c(1.3, 1)) &
    theme(legend.position = "bottom")
supPan <- supPan +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

## save the panel
ggsave("output/benchmark_supplementary.pdf", supPan, width = 180,
       height = 80, units = "mm",
       scale = 1.5, device = cairo_pdf)

