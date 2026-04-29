# =========================================================================
# Script Name: pvObs.R
# Description: This script ensures the covariate matrix has an intercept, 
#              fits the null model (phenotype vs. covariates without genetics), 
#              and calculates the observed p-value for the genetic association test.
# =========================================================================

# -------------------------------------------------------------------------
# 1. Non-Genetic Covariate Matrix (X) Formatting
# -------------------------------------------------------------------------
# Ensure that an intercept column (a column of 1s) is present in the 
# non-genetic covariate matrix (X). This is required for proper regression modeling.

if (is.null(x = X)) {
  # If no covariates were provided at all, create an X matrix that is 
  # simply an intercept vector with a length matching the number of individuals (nrow(snp)).
  X <- rep(x = 1.0, times = nrow(snp))
} else {
  # If X exists, check if it already contains an intercept column.
  ones <- rep(x = 1.0, times = nrow(snp))
  hasIntercept <- FALSE
  
  # Loop through each column of X to see if any column is entirely composed of 1s.
  for( i in 1L:ncol(X)) {
    hasIntercept <- isTRUE(all.equal(target = X[,i], current = ones ))
    if (hasIntercept) break # Stop searching once an intercept is found
  }
  
  # If no intercept column was found in the loop above, append one to the 
  # beginning of the matrix using cbind.
  if (!hasIntercept) {
    if (verbose) cat("\nintercept added to non-genetic covariate matrix\n")
    X <- cbind(1.0, X)
  }
}

# -------------------------------------------------------------------------
# 2. Fit the Null Model
# -------------------------------------------------------------------------
# Fit the null model, which regresses the phenotype (yy) on the non-genetic 
# covariates (X) without considering the genetic markers. This sets the baseline 
# for the score test.   
nullResult <- .newNullResult(trait = trait, 
                             X = X, 
                             yy = y, 
                             verbose = verbose)

##############################################
# 3. Calculate the Observed P-Value
##############################################  
# Call the main computational function to perform the variance component 
# score test. It tests the null model residuals against the genetic kernel 
# matrix (ZZ1), which has been weighted by the Hi-C interaction vector (rVec).
result <- mainCode(Rmc = rVec, 
                   snp = snp, 
                   kernel = kernel, 
                   nullResult = nullResult, 
                   X = X, 
                   pvMethod = pvMethod,
                   ZZ1 = ZZ1)

# If verbose logging is enabled, print the final observed p-value to the console.
if (verbose) cat(result$pvObs, " ")