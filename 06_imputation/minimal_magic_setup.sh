#!/bin/bash

# Minimal MAGIC Setup - Just the essentials
# This installs only what's absolutely needed for MAGIC to work

echo "=== Minimal MAGIC Setup ==="

# Navigate to 06_imputation directory
cd "06_imputation"

# Step 1: Remove existing environment if it exists and create fresh one
echo "Creating fresh minimal conda environment..."
conda env remove -n magic_env -y 2>/dev/null || true
conda create -n magic_env python=3.9 -y

# Activate environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate magic_env

# Step 2: Install compatible packages from conda-forge (avoids binary compatibility issues)
echo "Installing compatible packages from conda-forge..."
conda install -y -c conda-forge numpy=1.24.3 pandas=2.0.3 scipy matplotlib

# Step 3: Install MAGIC dependencies via conda/pip carefully
echo "Installing MAGIC dependencies..."
pip install scprep future tasklogger graphtools

# Step 4: Install MAGIC package directly
echo "Installing MAGIC package..."
cd MAGIC/python
pip install .
cd ../..

# Step 5: Test MAGIC (with error handling)
echo "Testing MAGIC..."
python -c "
try:
    import magic
    print('✓ MAGIC version:', magic.__version__)
    print('✓ MAGIC import successful!')
except Exception as e:
    print('✗ MAGIC import failed:', e)
    exit(1)
"

# Step 5: Install R MAGIC package
echo "Installing R MAGIC package..."
R --slave << 'EOF'
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}
install.packages("./MAGIC/Rmagic", repos = NULL, type = "source")
library(Rmagic)
cat("✓ R MAGIC installed and loaded\n")
renv::snapshot()
EOF

# Step 6: Create minimal R setup
cat > "minimal_magic_setup.R" << 'EOF'
# Minimal MAGIC Setup for R
reticulate::use_condaenv("magic_env", required = TRUE)
library(Rmagic)
cat("✓ Minimal MAGIC environment ready\n")
EOF

# Update .Rprofile
echo "" >> .Rprofile
echo "# Minimal MAGIC Setup" >> .Rprofile
echo "source('minimal_magic_setup.R')" >> .Rprofile

echo ""
echo "=== Minimal Setup Complete! ==="
echo "To use MAGIC:"
echo "1. conda activate magic_env"
echo "2. Start R in this directory"
echo ""
echo "Test with: python -c 'import magic; print(magic.__version__)'"
