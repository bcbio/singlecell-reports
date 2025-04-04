Templates with ![](https://img.shields.io/badge/status-stable-green) revision indicates that the components or processes have undergone comprehensive parameterization and testing.

Templates with ![](https://img.shields.io/badge/status-alpha-yellow) revision indicates that the components or processes are currently being tested. There is some test data available, but there are parameters that need to be set up manually within the code.

Templates with ![](https://img.shields.io/badge/status-draft-grey) revision indicates that the components or processes are not fully tested. There is no test data available, parameters need to be set up manually within the code, and specific code changes are required based on the data used.

Read [main page](https://github.com/bcbio) to know how to collaborate with us.

# Guideline for scRNAseq analysis

Make sure there is a valid project name, and modify [`information.R`](information.R) with the right information for your project. You can use this file with any other Rmd to include the project/analysis information.

-   Set the working directory to this file level. We recommend to use **Projects** in Rstudio.
-   Use [`install_dependencies.R`](install_dependencies.R) to install all packages used in these reports.

## nf-core

`cellranger` outputs are in your [output directory](https://nf-co.re/scrnaseq/4.0.0/docs/output/#cellranger) under `results/cellranger`.

## running cell-ranger by yourself

[`pre-process-w-cellranger.md`](pre-process-w-cellranger.md) contains step by step guidelines on how to run cellranger.

Then, the [`scripts/seurat_init.R`](scripts/seurat_init.R) script contains all the pieces to go from cellranger output to Seurat obj. It is assuming a mouse genome. Alternatively, if you have an especially large single cell dataset, you may wish to construct a Seurat object where the counts matrix remains stored on-disk and is not loaded into R memory. The [`scripts/seurat_ondisk.R`](scripts/seurat_ondisk.R) script can assist with this.

# Quality Assessment

## scATAC

The Rmd that helps to visualize ATAC metrics is [`01_quality_assessment/scATAC_QC.Rmd`](01_quality_assessment/scATAC_QC.Rmd).

## scRNA

Currently we are working on deploying a shiny app to inspect the single cell object and find the best cut-offs for filtering. The Rmd that helps to visualize the before and after is [`01_quality_assessment/scRNA_QC.Rmd`](01_quality_assessment/scRNA_QC.Rmd).

# Integration

[`02_integration/norm_integration.rmd`](02_integration/norm_integration.rmd) is a template with guidelines on how to work with multiple samples. It compares log2norm vs SCT, work with SCT by samples to remove batch biases better, provide options for integration between CCA and Harmony. As last step, it contains cell type clustering and visualization to help decide the best parameters.

# Differential expression

## pseudobulk approach

[`02_differential_expression/scRNA_pseudobulk.Rmd`](02_differential_expression/scRNA_pseudobulk.Rmd) is a template that performs pseudobulk differential expression analysis using DESeq2. It takes as input:

-   information.R, containing basic facts about the project, analysis, and scientific question
-   Seurat object containing both raw counts and log normalized data
-   a metadata variable of interest and two levels of that variable to be compared
-   a column in the seurat object metadata containing cluster names/numbers
-   the name/number of the cluster of interest

The template aggregates single cell counts to the sample x cluster x factor_of_interest pseudobulk level, subsets to the cluster of interest, performs the desired statistical comparisons, and visualizes the results.

## MAST - single cell approach

[`02_differential_expression/scRNA_MAST.Rmd`](02_differential_expression/scRNA_MAST.Rmd) is a template to visualize differentially expressed genes (DEG) results generated from MAST analysis. Main visualizations include:

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
