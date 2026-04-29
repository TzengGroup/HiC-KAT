# =========================================================================
# Script Name: Example.R
# Description: An example of execution of HiC-KAT
# =========================================================================

# Set working directory to the current directory
setwd(getwd())

# Load simulated genotype data and phenotype/covariate data
load("Dummy_Genotype_201_2000.RData")
load("yx.RData")

# Define the focal locus (SNP/region of interest)
rslocus = 95055000

# Create a vector of loci spanning +/- 1 Megabase (100 bins of 10,000 bp) around the focal locus
loci <- seq((rslocus - 100*10000), (rslocus + 100*10000), 10000)

# load HICKAT Function 
source("../code/HiC-KAT.R")

# -------------------------------------------------------------------------
# Execute the Function
# -------------------------------------------------------------------------
HICKAT(sample=2000, q.therashold=0.05, HiC.data = "LCL", which.chr = "chr5", rslocus = 95055000, 
       input.genotype=simulated_data, loci.vec = loci, cv.vec = c(1/4, 1/5, 1/6, 1/7, 1/8, 1/100))