#!/bin/bash

#SBATCH --job-name=NK_MAST      # Job name
#SBATCH --partition=short            # Partition name
#SBATCH --time=0-04:00                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --mem=16G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)
#SBATCH --array=1-2%2

module load gcc/9.2.0 imageMagick/7.1.0 geos/3.10.2 cmake/3.22.2 R/4.3.1 fftw/3.3.10 gdal/3.1.4 udunits/2.2.28 boost/1.75.0 python/3.9.14 hdf5/1.14.0

clusters=(2 3)
cluster=${clusters[$SLURM_ARRAY_TASK_ID-1]}

Rscript MAST_analysis.R \
  --seurat_obj=https://github.com/bcbio/bcbioR-test-data/raw/refs/heads/main/singlecell/tiny.rds \
  --resolution_column=integrated_snn_res.0.4 \
  --cluster_name="${cluster}" \
  --contrast=age \
  --outputDir=out