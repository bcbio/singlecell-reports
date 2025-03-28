library(Seurat)
library(qs)
library(BPCells)

# Goal: use BPCells package to work with very large single cell datasets. Allows counts matrix to
# be accessed on-disk in bit-packing compressed format rather than loaded into R memory, 
# decreasing the size of the Seurat object

# For more info, see Seurat vignette "Using BPCells with Seurat Objects"

################################### Scenario 1: ############################################
# From cellranger output of a single sample, create a Seurat object with counts matrix stored on disk
############################################################################################

# Load counts matrix
counts_matrix_h5 <- open_matrix_10x_hdf5(
  path = "/path/to/cellranger/output/raw_feature_bc_matrix.h5"
)
# Write the bitpacked matrix to a directory
counts_dir_name <- 'path/to/directory/where/you/want/to/save/bitpacked/counts/matrices'
write_matrix_dir(
    mat = counts_matrix_h5,
    dir = counts_dir_name
)

# Now that we have the bp matrix on disk, we can load it
counts_mat_bp <- open_matrix_dir(dir = counts_dir_name)

# Create Seurat Object
seurat <- CreateSeuratObject(counts = counts_mat_bp)
qsave(seurat, '/path/where/you/want/to/save/seurat/object/_.qs')

################################# Scenario 2: #############################################
# From an existing Seurat object, convert the counts and data layers (if present) to 
# be stored on disk
###########################################################################################

seurat <- qread('path/to/existing/seurat/object/_.qs')
assay_name <- 'RNA' # could be different if working with spatial data
counts_dir_name <- 'path/to/directory/where/you/want/to/save/bitpacked/counts/matrices'
data_dir_name <- 'path/to/directory/where/you/want/to/save/bitpacked/data/matrices'

# write existing matrices to bp format
write_matrix_dir(mat = seurat[[assay_name]]$counts, dir = counts_dir_name)
write_matrix_dir(mat = seurat[[assay_name]]$data, dir = data_dir_name)

# load bp matrices
counts.mat <- open_matrix_dir(dir = counts_dir_name)
data.mat <- open_matrix_dir(dir = data_dir_name)

# replace existing matrices with paths to bp matrices
seurat[[assay_name]]$counts <- counts.mat
seurat[[assay_name]]$data <- data.mat

# save newly bitpacked object
qsave(seurat, '/path/where/you/want/to/save/seurat/object/_.qs')
