# Hi-C Informed Kernel Association Test (HICKAT)

## Overview

The Hi-C Informed Kernel Association Test (HICKAT) is an advanced statistical framework designed to integrate three-dimensional (3D) genome architecture into variant-set association analyses for whole-genome sequencing (WGS) data.  

While traditional variant-set analyses (like SKAT) aggregate signals to improve statistical power for rare variants, they typically ignore spatial chromatin interactions. HICKAT extends this by offering a gene-agnostic, whole-genome testing strategy. It adaptively incorporates information from loci that physically interact with a target variant set in 3D space, ensuring that interacting loci with higher contact confidence contribute more to the association test. 

## Core Methodology

The HICKAT workflow consists of four primary steps:  

1. Obtain Hi-C Contact Confidence: The genome is partitioned into 10-kb regions (loci). Hi-C contact data is evaluated using q-values (via tools like Fit-Hi-C) to quantify the statistical significance of spatial proximity between locus pairs.
2. Convert q-values to Borrowing Weights: Contact q-values are first transformed into Z-scores. For significant interactions ($q < 0.05$), these Z-scores are mapped to a weight value ($w$) between 0 and 1 using a Gamma cumulative distribution function (CDF). A controlling parameter, $c$, regulates the degree of information borrowing.
3. Construct the Hi-C Informed Kernel: Using a Kernel Machine (KM) regression framework, genetic similarity is calculated. For a target locus, neighboring loci within a $\pm 1$ Mb window are incorporated into the kernel, weighted by their spatial interaction confidence ($w$).
4. Perform Association Testing: The test evaluates the null hypothesis across a grid of borrowing scales (e.g., $c \in \{1/100, 1/7, 1/6, 1/5\}$), where $c=1/100$ effectively reduces the model to the baseline SKAT test without spatial borrowing. The p-values from these different bandwidths are then combined into a single robust p-value using the Aggregated Cauchy Association Test (ACAT). 

## Project Structure and Scripts

+ HiC-KAT.R (Main Execution): The primary entry point. It defines the HICKAT wrapper function, iterates over a vector of spatial bandwidths (cv.vec), and calculates the final combined ACAT p-value.  
+ Essential functions and data_cluster.R: Sets up the global environment. It loads necessary R packages (rARPACK, CompQuadForm) and helper functions (like the Gamma CDF pga used for weight conversion).
+ mainCode.R: Performs the core mathematical operations, deriving the asymptotic p-value using eigenvalue decomposition on the kernel matrices.
+ snp_MAF.R: Identifies genomic bins interacting with the focal locus within the $\pm 1$ Mb window , extracts the sparse genotype matrices, and filters SNPs based on a Minor Allele Frequency (MAF) threshold.
+ Zvec.R: Maps the locus-level Hi-C interaction weights (Z-scores) down to the individual SNP level, creating the interaction vector Zvec.
+ kernel_trait2.R: Initializes the specific kernel matrices (linear or burden) and configures marker weighting. It applies the spatial weights to compute the factorized local kernel matrix (ZZ1).
+ pvObs.R: Prepares the non-genetic covariate matrix (ensuring an intercept is present), fits the null model, and triggers the score-like test to calculate the observed p-value.

## Example Execution

The `Example/Example.R` script provides a clear template for how to initialize the environment, load the necessary datasets, and execute the HICKAT pipeline. 

Here is a breakdown of the execution steps:

1. **Environment & Data Setup:**
   First, the script sets the working directory and loads the simulated genotype matrices (`Dummy_Genotype_201_2000.RData`) along with the phenotype and covariate data (`yx.RData`).
2. **Target Region Definition:**
   The focal locus (e.g., `95055000`) is defined. The script then generates a spatial window (`loci`) spanning $\pm 1$ Megabase (100 bins of 10,000 bp each) around this focal point to determine which neighboring loci will be evaluated for 3D interactions.
3. **Function Sourcing:**
   The main pipeline function is loaded into the environment by sourcing `HiC-KAT.R`.
4. **Execution:**
   The `HICKAT` function is called with the following key parameters:
   * `sample = 2000`: Number of individuals in the dataset.
   * `q.therashold = 0.05`: The False Discovery Rate (FDR) threshold for filtering rare variants.
   * `HiC.data = "LCL"`: The cell line or tissue type used to fetch the appropriate Hi-C contact map.
   * `cv.vec = c(1/4, 1/5, 1/6, 1/7, 1/8, 1/100)`: The grid of spatial bandwidth parameters used to scale the interaction weights. (Note: `1/100` serves as the baseline SKAT model without spatial borrowing).

```R
# Example execution call
HICKAT(sample = 2000, 
       q.therashold = 0.05, 
       HiC.data = "LCL", 
       which.chr = "chr5", 
       rslocus = 95055000, 
       input.genotype = simulated_data, 
       loci.vec = loci, 
       cv.vec = c(1/4, 1/5, 1/6, 1/7, 1/8, 1/100))
