## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  # install development version of package
#  install.packages("devtools")
#  library("devtools")
#  devtools::install_github("ajbass/sffdr")

## -----------------------------------------------------------------------------
library(sffdr)
set.seed(123)
data(bmi)
head(sumstats)
p <- sumstats$bmi
z <- as.matrix(sumstats[,-1])

## -----------------------------------------------------------------------------
# Create model: choose knots at small quantiles
mpi0 <- pi0_model(z = z,
                  knots = c(0.01, 0.025, 0.05, 0.1))

# Estimation of functional pi0 using the design matrix in mpi0
fpi0 <- fpi0est(p = p,
                z = mpi0$zt,
                pi0_model = mpi0$fmod)

## ----eval = FALSE-------------------------------------------------------------
#  # Create design matrix (can include other variables (e.g., MAF) or specify more complicated models)
#  fmod <- "~ns(bfp, knots = c(0.01, 0.025, 0.05, 0.1))"
#  fpi0_mod <- fpi0est(p = p,
#                      z = mpi0$zt,
#                      pi0_model = fmod)

## -----------------------------------------------------------------------------
# apply sfFDR
sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0)   

# plot significance results
plot(sffdr_out, rng = c(0, 1e-6))

# Functional P-values, Q-values, and local FDR
fp <- sffdr_out$fpvalues
fq <- sffdr_out$fqvalues
flfdr <- sffdr_out$flfdr

## -----------------------------------------------------------------------------
# Boolean to specify which SNPs are independent (e.g., pruning)
# All SNPs are LD-independent in this example data set 
indep_snps <- rep(TRUE, length(p))

# Create model 
mpi0 <- pi0_model(z = z,
                  indep_snps = indep_snps,
                  knots = c(0.01, 0.025, 0.05, 0.1))

# Estimation fpi0 using design matrix from mpi0
fpi0 <- fpi0est(p = p,
                z = mpi0$zt,
                indep_snps = indep_snps,
                pi0_model = mpi0$fmod)

# Estimate FDR quantities and functional p-value
sffdr_out <- sffdr(p,
                   fpi0 = fpi0$fpi0,
                   indep_snps = indep_snps)

