## ----preamble
# Differential analysis CyTOF workflow with CATALYST
# done by Helena Crowell, mods here by Mark Robinson
# setwd("/Users/mark/projects/cytof/Neutrophils/CATALYST_differential_analysis")
# devtools::install_github("HelenaLC/CATALYST", ref="f1000")
# knitr::spin("run_Samanta_01-25-2018.R")
# knitr::spin("run_Samanta_01-25-2018_20feb2018.R")
# file.symlink("01-25-2018/script.R", "run_Samanta_01-25-2018_24jun2018.R")
# knitr::spin("run_Samanta_01-25-2018_24jun2018.R")
# file.symlink("01-25-2018/script.R", "run_Samanta_01-25-2018_24jun2018_with_pdf_05apr2019.R")
# knitr::spin("run_Samanta_01-25-2018_24jun2018_with_pdf_05apr2019.R")
# knitr::spin("run_Samanta_01-25-2018_24jun2018_with_pdf_05apr2019.R", format="Rnw")

# file.symlink("01-25-2018/script.R", "run_Samanta_01-25-2018_24jun2018_with_pdf_09apr2019.R")
# knitr::spin("run_Samanta_01-25-2018_24jun2018_with_pdf_09apr2019.R")
# knitr::spin("run_Samanta_01-25-2018_24jun2018_with_pdf_09apr2019.R", format="Rnw")


# file.symlink("01-25-2018/script.R", "run_Samanta_01-25-2018_24jun2018_with_pdf_23apr2019.R")
# knitr::spin("run_Samanta_01-25-2018_24jun2018_with_pdf_23apr2019.R")
# knitr::spin("run_Samanta_01-25-2018_24jun2018_with_pdf_23apr2019.R", format="Rnw")

## ----libraries
suppressPackageStartupMessages({
  library(CATALYST)
  library(flowCore)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
})


## ----read in data
path <- "01-25-2018"

fcs_path <- file.path(path, "001_cleanfcs")
fcs <- list.files(fcs_path, ".fcs", full.names=TRUE)
fs <- read.flowSet(fcs)

## ----guess panel
panel <- CATALYST:::guessPanel(fs)  # so far, not exported
panel$use_channel[panel$antigen=="live"] <- 0

## ----read in metadata
md_path <- file.path(path, "metadata_samanta_CK_2018-01-25.xlsx")
(md <- readxl::read_xlsx(md_path))
md$filename <- file.path(fcs_path, md$filename)

m <- match(basename(fcs), basename(md$filename))
md <- md[m,]


stopifnot( all(basename(fcs) == basename(md$filename)) )

md_cols = list(file = "filename", id = "sample_id", 
               factors = c("condition1", "condition2"))

md$condition1 <- factor(md$condition1)

stopifnot( basename(fcs) == basename(md$filename) )

## ----construct daFrame
(cols_to_use <- as.logical(panel$use_channel))
se <- daFrame(fs, panel, md, cols_to_use, 
              md_cols = md_cols, cofactor=5)

o <- c(18,15,7,6,17,16,2,21,20,9,4,10,5,3,1,11,13,14,8,19,12)
se <- se[,o]

# little hack to recover text of condition factor
rowData(se)$condition1 <- factor(levels(md$condition1)[rowData(se)$condition1])
# rowData(se)$condition2 <- factor(levels(md$condition2)[rowData(se)$condition2])


## ----general data plotting, fig.height = 8, fig.width = 12----
# plotCounts(se, color_by="condition1")
# plotMDS(se, color_by="condition1")
(nrs <- plotNRS(se, color_by="condition1"))

## ----expr-heatmap, fig.height = 4, fig.width = 12----
plotExprHeatmap(se, anno = TRUE, color_by = "condition1")
plotExprHeatmap(se, anno = TRUE, scale = FALSE, color_by = "condition1")

## ----FlowSOM
nrs_df <- nrs$data %>% group_by(antigen) %>% summarise(nrs_mean=mean(NRS))
(cluster_on <- as.character(nrs_df$antigen[1:12]))
set.seed(1234)
# se <- cluster(se, cols_to_use=cluster_on, maxK = 40)
se <- cluster(se, cols_to_use=cluster_on)

# (t <- table(cluster_codes(se)[,"meta0"][cluster_ids(se)], sample_ids(se)))
(t <- table(cluster_codes(se)[,"meta20"][cluster_ids(se)], sample_ids(se)))

# round(100*t(t(t)/colSums(t)),2)
# round(rowSums(t)/sum(t)*100, 2)


## ----manual merging
# merging_table <- readxl::read_xlsx(file.path(path, "CK_2018-01-25merge.xlsx"))  
merging_table <- readxl::read_xlsx(file.path(path, "CK_2018-06-23merge3.xlsx"))
sem <- mergeClusters(se, k="meta20", table=merging_table, id="merging3")
# 
(t <- table(cluster_codes(sem)[,"merging3"][cluster_ids(sem)], sample_ids(sem)))
round(100*t(t(t)/colSums(t)),2)
round(rowSums(t)/sum(t)*100, 2)

# clust <- cluster_codes(sem)[,"merging3"][cluster_ids(sem)]
# med_exprs <- data.frame(CATALYST:::scale_exprs(CATALYST::exprs(sem)), cluster_ids=clust) %>%
#   group_by(cluster_ids) %>%
#   summarize_all(funs(median))
# 
# barplot(as.matrix(df["population 6",type_markers(sem)]))
# 
# df <- med_exprs %>% column_to_rownames("cluster_ids") %>% as.data.frame()f

# med_exprs <- data.frame(CATALYST::exprs(sem), cluster_ids=clust) %>% 
#   group_by(cluster_ids) %>% 
#   summarize_all(funs(median))


## ----plot-abundances, fig.height = 12, fig.width = 18----
plotAbundances(se, k="meta20", group = "condition2")
# plotAbundances(sem, k="meta40", group = "condition2")
plotAbundances(sem, k="merging3", group="condition2")


## ----plot-cluster-heatmap, fig.height = 9, fig.width = 12----
plotClusterHeatmap(x=sem, scale=TRUE, k="meta20", m="merging3",
                   hm2="abundances", draw_freqs = TRUE)


## ----plot-merged-heatmap, fig.height = 6, fig.width = 12----
plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="abundances", draw_freqs = TRUE)
plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="state_markers", draw_freqs = TRUE)

# plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="abundances", 
#                    draw_freqs = TRUE, split_by = "condition1")
# plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="abundances", 
#                      draw_freqs = TRUE, split_by = "condition2")
plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="abundances",
                     draw_freqs = TRUE, split_by = "sample_id")
# 
# # plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="state_markers", 
# #                    draw_freqs = TRUE, split_by = "condition1")
# # plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="state_markers", 
# #                    draw_freqs = TRUE, split_by = "condition2")
plotClusterHeatmap(x=sem, scale=TRUE, k="merging3", hm2="state_markers",
                   draw_freqs = TRUE, split_by = "sample_id")

## ----tSNE
# se2 <- tSNE(sem, n=1000)
se2 <- runDR(sem, dr = "TSNE", rows_to_use = 1000)
# se2 <- runDR(sem, dr = "TSNE", rows_to_use = 1000)
se2 <- runDR(se2, dr = "UMAP", rows_to_use = 1000)


## ----tSNE-facet, fig.height=8, fig.width=12
# plotDR(se2, color="meta30", facet = "sample_id")
# plotDR(se2, dr = "UMAP", color="meta30", facet = "sample_id")
# plotDR(se2, color="meta20", facet = "sample_id")
plotDR(se2, color="merging3", facet = "sample_id")
# plotDR(se2, dr = "UMAP", color="meta20", facet = "sample_id")
plotDR(se2, dr = "UMAP", color="merging3", facet = "sample_id")

sessionInfo()
