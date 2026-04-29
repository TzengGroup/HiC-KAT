#' Calculate the local kernel matrix and its rank factorization
#'
#' @param Rmc vector subset of variant similarity matrix R for the mth marker 
#'   and cvth c value (R[,m,cv])
#' @param snp nxM genotype snp matrix
#' @param kernel object of class Kernel
#' @param nullResult fit results of null model
#' @param X matrix of non-genetic covariates
#' @param pvMethod object of class pvMethod
#' @param ZZ1 pre-calculated local kernel matrix factorization 
#' @param ... ignored
#'
#' @return a list containing
#'   pvObs : p-value of observed test statistic
#'
#' @keywords internal
mainCode <- function(Rmc, 
                     snp, 
                     kernel, 
                     nullResult, 
                     X, 
                     pvMethod, 
                     ZZ1, ...) {

  # Calculate the influence function (psi). This typically represents the 
  # individual-level contributions to the score test statistic under the null model.
  psi <- .calcPsi(nullResult = nullResult, X = X, ZZ1 = ZZ1)

  # Perform eigen decomposition on the empirical covariance matrix of the 
  # influence function (crossprod(psi) / sample_size). 
  # We use `eigs_sym` (likely from the rARPACK package) to efficiently compute 
  # the top 'k' eigenvalues of a symmetric matrix without computing the full decomposition.
  ee <- tryCatch(expr = suppressWarnings(expr = eigs_sym(A = crossprod(x=psi)/nrow(x = snp),
                                                         k = ncol(x = psi),
                                                         which = "LM", # "LM" requests eigenvalues with the Largest Magnitude
                                                         opts = list("retvec" = FALSE))), # Set to FALSE because we only need values, not eigenvectors
                 error = function(e){
                   print(x = e$message)
                   stop("unable to obtain eigen decomposition")
                 })


  # Filter out eigenvalues that are zero or effectively zero due to floating-point
  # numerical precision limits (threshold set at 1e-10).
  nonzeroev <- ee$value[abs(x = ee$value) > 1e-10]

  # Calculate the asymptotic p-value of the observed test statistic using 
  # the specified p-value calculation method (e.g., Davies method, CompQuadForm) 
  # and the non-zero eigenvalues of the asymptotic covariance matrix.
  pvObs <- .calcPV(method = pvMethod, 
                   psi = psi, 
                   ev = nonzeroev)
  
  return( list("pvObs" = pvObs) )
       
}