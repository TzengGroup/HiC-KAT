# =========================================================================
# Script Name: Essential functions and data_cluster.R
# Description: Initialization script for HiC-KAT analysis. 
#             This script prepares the environment by loading. 
#              required libraries, custom analytical functions, and defining 
#              basic statistical helper functions.
# =========================================================================

# -------------------------------------------------------------------------
# 1. Load Required R Packages
# -------------------------------------------------------------------------
library(rARPACK)
library(CompQuadForm)
library(methods)

# -------------------------------------------------------------------------
# 2. Load Custom Hi-C Package Functions
# -------------------------------------------------------------------------
# Loading functions 
# These represent the core mathematical/statistical models used in the analysis.
load("../Code/f_kernel_v2.RData")     # Functions related to kernel matrix construction
load("../Code/f_pvmethod_v2.RData")   # Functions for p-value calculations
load("../Code/f_null_v2.RData")       # Functions for null distribution modeling
load("../Code/f_main_calc_v2.RData")  # Main calculator/wrapper functions

#-------------------------------------------------------------------------
# 4. Helper Functions
# -------------------------------------------------------------------------

#' Given variance, calcualte gamma cdf 
#' 
#' This function re-parameterizes the standard gamma distribution. Instead of 
#' shape and rate, it takes the expected mean and variance to calculate the CDF.
#' 
#' @param xx The numeric value(s) at which to evaluate the CDF
#' @param mn The mean of the desired gamma distribution
#' @param vn The variance of the desired gamma distribution
#' @return The probability P(X <= xx) for the corresponding Gamma distribution
pga <- function(xx, mn, vn){ 
  b <- mn/vn      # Calculate rate parameter (beta)
  a <- mn*b       # Calculate shape parameter (alpha)
  return(pgamma(xx, a, b))
}

#' Check if a vector is empty
#' 
#' @param x A vector
#' @return TRUE if the vector has a length of 0, FALSE otherwise
vector.is.empty <- function(x) {
  return(length(x) == 0)
}

# Clean up any leftover 'tmp' variable from the global environment to free memory
rm(tmp)
