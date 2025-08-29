install.packages("BiocManager")
BiocManager::install("renv")
BiocManager::install("immunogenomics/presto")
BiocManager::install("sqjin/CellChat")
BiocManager::install(renv::dependencies(path = ".")[["Package"]])

create_snapshot <- function(folder, extra_libs = NULL) {
  pkgs <- renv::dependencies(path = folder)[["Package"]]
  if (!is.null(extra_libs)) {
    pkgs <- c(pkgs, renv::dependencies(path = extra_libs)[["Package"]])
  }
  renv::snapshot(lockfile = file.path(folder, "renv.lock"), packages = pkgs)
}

# Snapshots for each folder
create_snapshot("01_quality_assessment")
create_snapshot("02_integration")
create_snapshot("03_differential_expression")
create_snapshot("04_compositional")
create_snapshot("05_signaling")
create_snapshot("06_imputation")
