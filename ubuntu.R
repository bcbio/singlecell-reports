# Configure BiocManager to use Posit Package Manager
options(BioC_mirror = "https://packagemanager.posit.co/bioconductor/latest")

# Configure BiocManager to load its configuration from Package Manager
options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/latest/config.yaml")

# Set the Bioconductor version to prevent defaulting to a newer version
Sys.setenv("R_BIOC_VERSION" = "3.22")

# Configure a CRAN snapshot compatible with Bioconductor 3.22
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
