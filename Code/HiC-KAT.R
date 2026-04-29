# =========================================================================
# Script Name: HiC-KAT.R
# Description: Main execution script for the HICKAT (Hi-C Kernel Association Test) 
#              pipeline. It loads initial datasets, defines the target region, 
#              and runs a function to calculate and combine p-values across 
#              multiple spatial bandwidths.
# =========================================================================

# -------------------------------------------------------------------------
# HICKAT Function Definition
# -------------------------------------------------------------------------
HICKAT <- function(sample, q.therashold, HiC.data, which.chr, rslocus, input.genotype, loci.vec, cv.vec){
  
  # Load essential functions globally
  source("../Code/Essential functions and data_cluster.R")
  
  # Source data processing scripts locally within the function's environment
  source("../Code/mainCode.R", local = TRUE)
  source("../Code/snp_MAF.R", local = TRUE)
  source("../Code/Zvec.R", local = TRUE)
  
  # Set up kernel configuration parameters
  weighted = TRUE
  weight = NULL
  trait = "gaussian"       # Indicates a continuous phenotype
  kernel.tmp = 'linear'    # Linear kernel
  pvMethod.tmp = 'Liu'     # Davies/Liu method for p-value calculation
  verbose = TRUE
  vn <- 10                 # Variance parameter for the Gamma distribution used in Zvec weighting
  
  
  p <- NULL # Initialize empty vector to store p-values
  
  # Loop through each bandwidth parameter (cv)
  for (i in 1:length(cv.vec)){
    cv <- cv.vec[i]
    
    # Calculate the kernel weights and observed p-value for the current cv parameter
    source("../Code/kernel_trait2.R", local = TRUE)
    source("../Code/pvObs.R", local = TRUE)
    
    # Append the calculated observed p-value to the vector
    p <- c(p, result$pvObs)
  }
  
  # -------------------------------------------------------------------------
  # Aggregated Cauchy Association Test (ACAT)
  # -------------------------------------------------------------------------
  # Combine the p-values from the different bandwidths into a single robust p-value.
  # This accounts for the correlation between tests performed at different cv values.
  cct <- tan((0.5-p)*pi)
  cct.mn <- mean(cct, na.rm = TRUE)
  p_acat <- 0.5 - atan(cct.mn) / pi
  
  return(p_acat)
}