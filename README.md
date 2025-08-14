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

The Rmd that helps to visualize ATAC metrics is ![](https://img.shields.io/badge/status-alpha-yellow) [`01_quality_assessment/scATAC_QC.qmd`](01_quality_assessment/scATAC_QC.qmd).

## ![](https://img.shields.io/badge/status-stable-green) [scRNA quality assessment](01_quality_assessment/scRNA_QC.qmd) 👀 [example](https://bcbio.github.io/singlecell-reports/01_quality_assessment/scRNA_QC.html)

Currently we are working on deploying a shiny app to inspect the single cell object and find the best cut-offs for filtering. This tempalte helps to visualize the before and after is ![](https://img.shields.io/badge/status-alpha-yellow).

# ![](https://img.shields.io/badge/status-stable-green) [scRNA integration](02_integration/norm_integration.qmd) 👀 [example](https://bcbio.github.io/singlecell-reports/02_integration/norm_integration.html)

This template has guidelines on how to work with multiple samples. It compares log2norm vs SCT, work with SCT by samples to remove batch biases better, provide options for integration between CCA and Harmony. As last step, it contains cell type clustering and visualization to help decide the best parameters.

# Differential Expression

Read full documentation at [03_differential_expression/README.md](03_differential_expression/README.md).

-   ![](https://img.shields.io/badge/status-stable-green) [scRNA_pseudobulk](03_differential_expression/scRNA_pseudobulk.qmd) is a template that performs pseudobulk differential expression analysis using DESeq2. 👀 [See an example](https://bcbio.github.io/singlecell-reports/03_differential_expression/scRNA_pseudobulk.html)

- ![](https://img.shields.io/badge/status-stable-green) [MAST scRNA](03_differential_expression/scRNA_MAST.Rmd) is a template to visualize differentially expressed genes (DEG) results generated from MAST analysis. 👀 [See an example](https://bcbio.github.io/singlecell-reports/03_differential_expression/scRNA_MAST.html)

# Compositional Analysis

![](https://img.shields.io/badge/status-draft-grey) [`04_compositional/propeller.Rmd`](04_compositional/propeller.Rmd) and [`04_compositional/sscomp.Rmd`](04_compositional/sccomp.Rmd) are templates to run compositional analysis with two different methods. `comp.png` is an example of `sccomp.Rmd` analysis.
