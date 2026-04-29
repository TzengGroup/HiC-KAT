# =========================================================================
# Script Name: snp_MAF.R
# Description: This script identifies Hi-C interacting regions for a specific 
#              target locus, extracts the corresponding genotype data, filters 
#              SNPs based on Minor Allele Frequency (MAF), and formats the 
#              resulting matrix for downstream kernel analysis.
# =========================================================================

# Source required variables and base parameters locally
source("../Code/HiC_Info.R", local = TRUE)

# loci <- seq(ifelse((rslocus - 100*10000) >= min(mids), (rslocus - 100*10000),  min(mids)),
#             ifelse((rslocus + 100*10000) <= max(mids), (rslocus + 100*10000),  max(mids)), 
#             10000)


# -------------------------------------------------------------------------
# 1. Identify Interacting Loci (Hi-C Data)
# -------------------------------------------------------------------------
# Extract all unique fragment midpoints that interact with the focal locus (rslocus).
# It checks both fragmentMid1 and fragmentMid2 to ensure all interactions are caught,
# and also explicitly includes the focal rslocus itself.
mids.hic.rslocus <- unique(c(TargetRegion$fragmentMid2[which(TargetRegion$fragmentMid1==rslocus)],
                             TargetRegion$fragmentMid1[which(TargetRegion$fragmentMid2==rslocus)],
                             rslocus))

# Filter interacting loci to only keep those within a 1 Megabase (1,000,000 bp) 
# distance from the focal rslocus (100 bins * 10,000 bp resolution).
mids.hic.rslocus <- mids.hic.rslocus[which(abs(mids.hic.rslocus - rslocus)/10000 <= 100)]

# Find the specific indices of these valid interacting midpoints within the 
# master vector of all theoretical loci (loci.vec)
index.hic <- which(loci.vec %in% mids.hic.rslocus)

# -------------------------------------------------------------------------
# 2. Extract and Filter Genotype Data
# -------------------------------------------------------------------------
# Initialize empty data frames to store aggregated SNPs and locus summaries.
# (Populated initially with NAs to establish row lengths).
snp <- data.frame(rep(NA,sample))
Zmatrix <- data.frame(rep(NA,length(index.hic)))

Zmatrix$locus <- loci.vec[index.hic]

# Loop through each valid interacting locus to extract and filter SNPs
for (i in 1:length(index.hic)){
  
  # Extract the genotype matrix for the current locus
  tmp_geno <- input.genotype[[i]]$genotypes
  
  # Calculate Minor Allele Frequency (MAF). 
  # Assuming genotypes are coded as 0, 1, 2, dividing colMeans by 2 gives the allele frequency.
  MAF_geno <- colMeans(tmp_geno) / 2.0
  
  # Filter SNPs: Keep only those with a MAF below the specified threshold (rare variants)
  # Append these filtered SNPs to the master 'snp' data frame
  snp <- cbind(snp, tmp_geno[, which(MAF_geno < q.therashold)])
  
  # Record summary statistics in the Zmatrix for this locus
  Zmatrix$snp[i] <- dim(tmp_geno)[2]                             # Total SNPs initially available
  Zmatrix$included.SNP[i] <- length(which(MAF_geno < q.therashold)) # Total SNPs kept after filtering
}

# Clean up the initialization columns (the first column of NAs)
snp <- snp[,-1]
Zmatrix <- Zmatrix[,-1]

# Free up system memory by removing the large, unfiltered genotype list
rm(input.genotype)


# -------------------------------------------------------------------------
# 3. Format Data for Downstream Analysis
# -------------------------------------------------------------------------

# Convert the final SNP data frame into a standard matrix
snp <- as.matrix(snp)

# Convert to a sparse matrix to heavily optimize memory usage and speed up 
# downstream algebraic operations (highly recommended for large genomic datasets)
snp <- Matrix::Matrix(data = snp, sparse = TRUE)

# Calculate final dataset dimensions and MAFs for the retained markers
nmarker <- ncol(x = snp)
MAF <- colMeans(x = snp) / 2.0

# Print the final dimensions of the sparse SNP matrix to the console
cat("dim of snp", dim(snp))