## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----eval = FALSE-------------------------------------------------------------
# # install development version of package
# install.packages("devtools")
# devtools::install_github("ajbass/sffdr")

## -----------------------------------------------------------------------------
library(sffdr)
data(bmi)

# Define primary p-values and informative p-values
p <- sumstats$bmi
z <- as.matrix(sumstats[,-1])

head(sumstats)

## -----------------------------------------------------------------------------
# Create model
mpi0 <- pi0_model(z)

# Estimation of functional pi0
fpi0 <- fpi0est(p, mpi0)

## -----------------------------------------------------------------------------
# Apply sfFDR
sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0)   

# Extract results
results <- data.frame(
  p_value = p,
  fp_value = sffdr_out$fpvalues,
  fq_value = sffdr_out$fqvalues,
  flfdr = sffdr_out$flfdr
)

# Plot significance results (focusing on high-significance SNPs)
plot(sffdr_out, rng = c(0, 1e-6))

## ----eval = FALSE-------------------------------------------------------------
# # Fit model using LD-independent subset (indep_snps: TRUE if SNP is LD-independent)
# mpi0 <- pi0_model(z = z, indep_snps = indep_snps)
# fpi0 <- fpi0est(p, mpi0, indep_snps = indep_snps)
# 
# # Recommended: cross-fit by chromosome (chr: the chromosome of each SNP)
# sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0, folds = chr, nfolds = 5)
# 
# # Without LD-independent SNPs (slower): pi0 is estimated in the same folds
# fpi0 <- fpi0est(p, pi0_model(z = z), folds = chr, nfolds = 5)
# sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0, folds = chr, nfolds = 5)

## ----eval = FALSE-------------------------------------------------------------
# # LD score (log scale) added to the pi0 model built by pi0_model()
# zt <- cbind(mpi0$zt, ld = log(ld_score))
# fmod <- update(mpi0$fmod, ~ . + ns(ld, df = 3))
# fpi0 <- fpi0est(p, z = zt, pi0_model = fmod, indep_snps = indep_snps)
# sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0, folds = chr, nfolds = 5)

## ----eval = FALSE-------------------------------------------------------------
# w <- 1 / ld_score
# mpi0 <- pi0_model(z, weights = w)
# fpi0 <- fpi0est(p, mpi0, weights = w)
# sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0, folds = chr, nfolds = 5)

## ----eval = FALSE-------------------------------------------------------------
# mpi0 <- pi0_model_echisq(echisq, indep_snps = indep_snps)
# fpi0 <- fpi0est(p, mpi0, indep_snps = indep_snps)
# sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0, surrogate = echisq, folds = chr, nfolds = 5)

## ----eval = FALSE-------------------------------------------------------------
# mpi0 <- pi0_model_echisq(cbind(echisq = echisq, ld = ld_score), indep_snps = indep_snps)
# fpi0 <- fpi0est(p, mpi0, indep_snps = indep_snps)
# sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0, folds = chr, nfolds = 5)

## ----eval = FALSE-------------------------------------------------------------
# # Create model (can include other variables (e.g., MAF) or specify more complicated models)
# fmod <- "~ns(bfp, knots = c(0.01, 0.025, 0.05, 0.1))"
# fpi0_mod <- fpi0est(p = p,
#                     z = mpi0$zt,
#                     pi0_model = fmod)

