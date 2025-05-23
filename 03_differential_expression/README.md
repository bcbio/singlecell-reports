# Differential expression

## pseudobulk approach

![](https://img.shields.io/badge/status-alpha-yellow) [`02_differential_expression/scRNA_pseudobulk.Rmd`](02_differential_expression/scRNA_pseudobulk.Rmd) is a template that performs pseudobulk differential expression analysis using DESeq2. It takes as input:

-   information.R, containing basic facts about the project, analysis, and scientific question
-   Seurat object containing both raw counts and log normalized data
-   a metadata variable of interest and two levels of that variable to be compared
-   a column in the seurat object metadata containing cluster names/numbers
-   the name/number of the cluster of interest

The template aggregates single cell counts to the sample x cluster x factor_of_interest pseudobulk level, subsets to the cluster of interest, performs the desired statistical comparisons, and visualizes the results.

## MAST - single cell approach

![](https://img.shields.io/badge/status-alpha-yellow) [`02_differential_expression/scRNA_MAST.Rmd`](02_differential_expression/scRNA_MAST.Rmd) is a template to visualize differentially expressed genes (DEG) results generated from MAST analysis. Main visualizations include:

-   Group-level mean expression shown in heatmap;
-   Volcano plots highlighting top DEGs;
-   An adapted `seurat` Dotplot to show % of cells expressing top DEGs;
-   An integrated Violin-Box-Scatter (VBS) plot displaying the normalized expression of top DEGs per single cell across contrast groups.

We separately prepare the Rscript `02_differential_expression/MAST_analysis.R` to pre-compute the DEG results since MAST computation is relatively time-consuming. For our demo dataset, \~1000 cells and two group comparison, it took around 10-15 minutes.

To run this Rscript, you should input below parameters:

-   `--seurat_obj`: the seurat object containing your scRNA-seq data, raw counts stored in the layer `counts` of assay `RNA` is required
-   `--resolution_column`: the column name of your choice of clustering method + resolution chosen, for example `pca_res0.4`, this is only required if you want to subset your data
-   `--cluster_name`: the cluster you want to subset for MAST analysis
-   `--contrast`: the column name in your `seurat` cell metadata indicating the group you want to do contrast for
-   `--outputDir`: the output directory you will find your intermediate results for running `scRNA_MAST.Rmd` Default: `out`.

If none of the parameters are supplied, just run `Rscript 02_differential_expression/MAST_analysis.R` will run it on our demo dataset.

To obtain more informative console logs from MAST analysis, please consider running:

`Rscript --no-save --no-restore --verbose differential_expression/MAST_analysis.R > MAST.Rout 2>&1`

You will expect to have three main outputs in your specified output folder:

-   `processed_seurat.rds`: the processed seurat object containing log-normalized data
-   `MAST_RESULTS*`: the MAST modeling object to allow for further follow-up analysis of your own choice
-   `FULL_MAST_RESULTS_*`: full MAST differential expression analysis results in csv format
-   `SIG_MAST_RESULTS_padj<0.05*`: significant DEGs using 0.05 as the threshold for FDR
