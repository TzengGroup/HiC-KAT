# =========================================================================
# Script Name: kernel_trait2.R
# Description: This script sets up the kernel parameters, handles marker 
#              weighting (e.g., up-weighting rare variants using a Beta 
#              distribution), validates the phenotype trait type, and 
#              computes the local kernel matrix factorization using Hi-C 
#              derived spatial weights.
# =========================================================================

# -------------------------------------------------------------------------
# 1. Marker Weighting Setup
# -------------------------------------------------------------------------
if (!weighted) {
  # If the weighted flag is FALSE, explicitly set weight to NULL
  weight <- NULL
  if (verbose) cat("\nno weighting of markers\n")
} else {
  # If weighted is TRUE but no specific weights are provided, 
  # assign default weights using a Beta(1, 25) distribution based on the 
  # Minor Allele Frequency (MAF). This is a standard approach in sequence 
  # kernel association tests (SKAT) to give rare variants more influence.
  if (is.null(x = weight)) {
    weight <- dbeta(MAF,1,25)
  } else if (nmarker != length(x = weight)) {
    # Guardrail: Ensure the length of a custom weight vector exactly matches 
    # the total number of genetic markers (SNPs) in the dataset.
    stop("if provided, dimension of weight must be the # of markers")
  }
}

# -------------------------------------------------------------------------
# 2. Trait Validation
# -------------------------------------------------------------------------
# Standardize the trait string to lowercase to avoid case-sensitivity errors
trait <- tolower(x = trait)

# Restrict the analysis to supported phenotypic distributions: 
# 'binomial' (e.g., case/control disease status) or 'gaussian' (continuous traits)
 if (!{trait %in% c('binomial', 'gaussian')}) {
	stop("trait must be one of {'binomial', 'gaussian'}")
 }

# -------------------------------------------------------------------------
# 3. Method and Kernel Object Initialization
# -------------------------------------------------------------------------
# Map the temporary variables (defined in the parent script) to the working variables
kernel = kernel.tmp
pvMethod = pvMethod.tmp

# Instantiate the internal configuration objects for the p-value calculation 
# method and the kernel. 
pvMethod <- .newPvMethodObj(pvMethod = pvMethod, 
                            verbose = verbose)

kernel <- .newKernelObj(kernel = kernel, 
                        weight = weight,
                        verbose = verbose)   

# -------------------------------------------------------------------------
# 4. Hi-C Interaction Weight Vector (rVec) Construction
# -------------------------------------------------------------------------
# The 'cv' variable acts as a smoothing or bandwidth parameter.
if (cv < 1e-6) {
  # If cv is effectively zero, the script ignores the 3D chromosomal structure.
  # It creates a vector of zeros, placing a weight of 1.0 ONLY on the focal locus.
  rVec <- numeric(length = nmarker)
  rVec[msnplocation] <- 1.0
} else {
  # If cv > 0, it calculates spatial weights for the surrounding loci using 
  # the 'pga' function (the Gamma CDF defined earlier in your pipeline). 
  # This scales the Hi-C interaction z-scores (Zvec) into functional weights.
  rVec <- pga(Zvec,1/cv,vn)
  
  # Ensure the focal locus itself always retains a maximum interaction weight of 1.0
  rVec[msnplocation] <- 1.0
}         

# -------------------------------------------------------------------------
# 5. Local Kernel Calculation
# -------------------------------------------------------------------------
# Compute the factorized local kernel matrix (ZZ1) combining the genetic 
# data (snp), the base kernel configuration, and the Hi-C physical 
# interaction weights (rVec).
ZZ1 <- calcLocalKernel(Rmc = rVec,
                     snp = snp,
                     kernel = kernel)