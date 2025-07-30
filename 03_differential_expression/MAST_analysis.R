##### library loading #####
library(Seurat)
library(apeglm)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(MAST)
library(tidyverse)
library(data.table)
library(lme4)
library(R.utils)
library(glue)
library(optparse)
library(qs)

################################################################################
# parse parameters
################################################################################

options(stringsAsFactors = F)
option_list <- list(
  make_option("--seurat_obj", default = "https://github.com/bcbio/bcbioR-test-data/raw/refs/heads/main/singlecell/tiny.rds"),
  make_option("--resolution_column", default = "integrated_snn_res.0.4"),
  make_option("--cluster_name", default = "2"),
  make_option("--contrast", default = "age"),
  make_option("--outputDir", default = "out")
)
args <- parse_args(OptionParser(option_list = option_list))

invisible(list2env(args, environment()))
column <- contrast
system(glue("mkdir -p {outputDir}"))


message("[Preparing inputs for MAST modeling]")

################################################################################
# Read in provided seurat object
################################################################################

if (isUrl(seurat_obj)) {
  seurat <- readRDS(url(seurat_obj))
} else {
  seurat <- qread(seurat_obj)
}

message("Input seurat object: ", seurat_obj)
DefaultAssay(seurat) <- "RNA"
message("RNA is set as the default assay")
message("Column name of clustering to use: ", resolution_column)
Idents(object = seurat) <- resolution_column
message("Subset original seurat to be only cluster ", cluster_name, " for faster computing!")

################################################################################
# Prepare MAST inputs
################################################################################

data_subset <- subset(x = seurat, idents = cluster_name)

data_subset <- CreateSeuratObject(
  counts = data_subset[["RNA"]]$counts,
  project = "subset",
  meta.data = data_subset@meta.data
)


message("Natural log of raw counts with pseudobulk 1 used for MAST modeling")
sce <- as.SingleCellExperiment(data_subset)

# ##### Log-Normalize Seurat for visualization later #####
message("Total counts normalization and log1p transformation done to raw counts")
message("New layer Data added to the seurat object for visualization later!")
data_subset <- NormalizeData(
  data_subset,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1
)
qsave(data_subset, glue("{outputDir}/processed_seurat.qs"))

# log-normalize SCE
assay(sce, "log") <- log(counts(sce) + 1)

# scaling ngenes
cdr <- colSums(assay(sce, "log") > 0)
colData(sce)$cngeneson <- scale(cdr)

# Create new sce object (only 'log' count data)
sce.1 <- SingleCellExperiment(assays = list(log = assay(sce, "log")))
colData(sce.1) <- colData(sce)

# change to sca
sca <- SceToSingleCellAssay(sce.1)

message("Subset genes observed in at least 10% of cells")

expressed_genes <- freq(sca) > 0.1
sca_filtered <- sca[expressed_genes, ]

cdr2 <- colSums(SummarizedExperiment::assay(sca_filtered) > 0)

SummarizedExperiment::colData(sca_filtered)$ngeneson <- scale(cdr2)
SummarizedExperiment::colData(sca_filtered)$orig.ident <-
  factor(SummarizedExperiment::colData(sca_filtered)$orig.ident)
SummarizedExperiment::colData(sca_filtered)[[column]] <-
  factor(SummarizedExperiment::colData(sca_filtered)[[column]])

################################################################################
# MAST modeling
################################################################################

message("[MAST modeling with supplied contrasts]")

message("Note: this step is time-consuming!")

# define model coefficients to examine
SummarizedExperiment::colData(sca_filtered)[[contrast]] <- factor(SummarizedExperiment::colData(sca_filtered)[[contrast]])
comp_name <- levels(SummarizedExperiment::colData(sca_filtered)[[contrast]])[2]
lrt_names <- paste0(contrast, comp_name)

# fit model
formula_touse <- as.formula(paste0("~ ngeneson + (1 | orig.ident) + ", contrast))
zlmCond <- suppressMessages(MAST::zlm(formula_touse, sca_filtered,
  method = "glmer",
  ebayes = F, strictConvergence = FALSE
))

# examine coefficients
summaryCond_column <- suppressMessages(MAST::summary(zlmCond, doLRT = lrt_names))

message("[Main MAST computation done, result outputs]")

summary_cond_file <- paste0(outputDir, "/MAST_RESULTS_", cluster_name, "_", contrast, ".rds")
saveRDS(summaryCond_column, file = summary_cond_file)

message("Full MAST object saved to file ", summary_cond_file)

################################################################################
# collect MAST outputs
################################################################################

summaryDt_column <- summaryCond_column$datatable

fcHurdle_column <- summaryDt_column %>%
  filter(component %in% c("H", "logFC"), contrast %in% lrt_names) %>%
  dplyr::select(-component) %>%
  pivot_longer(!primerid & !contrast, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value), metric != "z") %>%
  pivot_wider(names_from = "metric", values_from = "value") %>%
  mutate(cluster_name = cluster_name) %>%
  relocate(cluster_name)

fcHurdle_column <- stats::na.omit(as.data.frame(fcHurdle_column))

fcHurdle_column <- lapply(lrt_names, function(curr_contrast) {
  fcHurdle_column %>%
    filter(contrast == curr_contrast) %>%
    mutate(fdr = p.adjust(`Pr(>Chisq)`, "fdr"))
}) %>% bind_rows()


full_res_file <- paste0(outputDir, "/FULL_MAST_RESULTS_", cluster_name, "_", contrast, ".csv")
write.table(fcHurdle_column, file = full_res_file, row.names = FALSE, sep = ",")

message("MAST summary results output to csv files")


fcHurdleSig_column <- fcHurdle_column %>%
  filter(fdr < 0.05) %>%
  arrange(fdr)

sig_res_file <- paste0(outputDir, "/SIG_MAST_RESULTS_padj_less_0.05_", cluster_name, "_", contrast, ".csv")
write.table(fcHurdleSig_column, file = sig_res_file, row.names = FALSE, sep = ",")

message("Significant MAST summary results output to csv files")
