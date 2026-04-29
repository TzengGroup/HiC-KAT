# Analysis parameters needed
# q.therashold: FDR threshold for significant interactions (Note: 'therashold' is a likely typo for 'threshold')
# HiC.data: Cell line/tissue type identifier
# which.chr: Target chromosome for analysis

# Load and format chromosome sizes
chromsize <- read.table("../Additional Data/hg19.chrom.sizes.txt")
chromsize <- chromsize[1:24,] # Keep only primary chromosomes (chr1-22, X, Y)
chromsize <- chromsize[-which(chromsize[,1]%in%c("chrY")),] # Remove chrY
chrom <- levels(as.factor(as.character(chromsize[,1]))) # Get list of remaining chromosome names


# Extract the specific length for the target chromosome
size <- chromsize[which(chromsize[,1] == which.chr),2]

# Conditionally load the appropriate FitHiC workspace (.RData) based on the specified cell type
if (HiC.data == "LCL"){
  load("../FitHiC/LCL - 10k resolution - hg19/LCL_FitHiC.RData")
}else if (HiC.data == "LV"){
  load("../FitHiC/Left Ventricle - 10k resolution/LeftVentricle_FitHiC.RData")
}

# Extract the FitHiC data specifically for the target chromosome 
# (assumes a list object named 'chr' was loaded from the .RData file)
TargetRegion <- chr[[which.chr]]
rm(chr)

# Standardize the q-value column name if dealing with the 'LV' dataset
if (HiC.data == "LV"){
  TargetRegion$q.value <- TargetRegion$fithic_qvalue
}

# Filter interactions to keep only those that meet the significance threshold
TargetRegion <- TargetRegion[TargetRegion$q.value < q.therashold,]

# Convert q-values to Z-scores
qvalue <- TargetRegion$q.value
# Use the inverse normal cumulative distribution function to get Z-scores
z <- qnorm(qvalue/2, lower.tail = F)

# Handle infinite Z-scores (which occur when q-value is exactly 0) by capping them 
# at the maximum finite Z-score in the dataset
z[which(z == Inf)] <- max(z[which(z!=Inf)])

# Apply a square root transformation to the Z-scores
z <- z**(1/2)
summary(z)
TargetRegion$"sqrtZ" <- z

# Clean up memory by removing intermediate variables
rm(qvalue, z)

# Generate theoretical 10kb resolution bin midpoints for the entire chromosome
mids <- seq(1, size, 10000)+4999
mids <- c(mids[-length(mids)], size) # Ensure the last bin extends to the chromosome size
mids <- as.integer(mids)

# Extract all unique bin midpoints that actually appear in the filtered Hi-C data
mids.hic <- unique(c(TargetRegion$fragmentMid1, TargetRegion$fragmentMid2))

# Align the maximum theoretical midpoint with the empirical maximum if they are very close
if (abs(max(mids.hic) - max(mids))<5000){
  mids[max(mids)] <- max(mids.hic)
}

# Identify which theoretical bins actually have data in the significant Hi-C interactions
include <- which(mids%in%mids.hic)