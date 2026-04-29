# =========================================================================
# Script Name: Zvec.R
# Description: This script isolates the SNP indices corresponding to a focal 
#              locus and constructs a vector of Z-scores (Zvec) that maps 
#              locus-level Hi-C interaction weights to individual SNPs.
# =========================================================================

##############################################
# 1. Locate the SNPs for the Focal Locus
##############################################
msnplocation <- NULL

# Calculate the cumulative sum of SNPs up to the focal locus (rslocus).
tmptmp <- sum(Zmatrix$included.SNP[0:which(Zmatrix$locus == rslocus)])

# Generate a sequence of column indices corresponding to the focal locus.
msnplocation <- c(msnplocation, seq(from = tmptmp + 1, 
                                    to = tmptmp + Zmatrix$included.SNP[which(Zmatrix$locus == rslocus)], 
                                    by = 1))

# Clean up intermediate variable
rm(tmptmp)

# Create a backup of the Zmatrix
Zmatrix_tmp <- Zmatrix

################################################
# 2. Construct the Z_m Locus Vector (Zvec)
################################################
tmp <- TargetRegion

# Initialize a vector of zeros with a length equal to the number of loci in Zmatrix
loci.vec <- rep(0, dim(Zmatrix)[1])

# If there are significant Hi-C interactions.
if (dim(tmp)[1] != 0){
  for (i in 1:dim(tmp)[1]){
    if (tmp[i,"fragmentMid1"]==rslocus){
      loci.vec[which(Zmatrix$locus == tmp[i,"fragmentMid2"])] <- tmp[i,"sqrtZ"]}
    else if (tmp[i,"fragmentMid2"]==rslocus){
      loci.vec[which(Zmatrix$locus == tmp[i,"fragmentMid1"])] <- tmp[i,"sqrtZ"]}
    
    # Expand the locus-level Z-scores to the SNP-level.
    # The 'rep' function duplicates each locus's Z-score 'included.SNP' times, 
    # ensuring every individual SNP gets the interaction weight of its parent locus.
    Zvec <- rep(loci.vec,times=Zmatrix$included.SNP)
  }
}else{
  # If there are no interactions (dim(tmp)[1] == 0), just expand the vector of zeros 
  # so it matches the total number of SNPs.
  Zvec <- rep(loci.vec,times=Zmatrix$included.SNP)
}

# Clean up memory by removing the temporary dataframe and locus vector
rm(tmp,loci.vec)
