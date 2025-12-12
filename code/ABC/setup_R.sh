#!/bin/bash
# Setup R 4.1.2 for RedSwampCrayfish_MISGP project with x86 library

# Get the current directory (should be rscABC)
PROJECT_DIR="/mnt/ufs18/rs-012/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/"
R_INSTALL="${PROJECT_DIR}/SRC/R-4.1.2-install"
R_LIBS="${PROJECT_DIR}/SRC/myR/x86_64-pc-linux-gnu-library/4.1"

# Set R environment to use our R 4.1.2
export PATH="${R_INSTALL}/bin:$PATH"
export LD_LIBRARY_PATH="${R_INSTALL}/lib/R/lib:$LD_LIBRARY_PATH"
export R_HOME="${R_INSTALL}/lib/R"
export R_LIBS_SITE="${R_LIBS}"
export R_LIBS_USER="${R_LIBS}"

echo "R 4.1.2 configured for RedSwampCrayfish_MISGP project"
echo "R installation: ${R_INSTALL}"
echo "R libraries: ${R_LIBS}"
echo "Using R version: $(R --version | head -1)"
echo "Available packages: $(ls ${R_LIBS} | wc -l)"
