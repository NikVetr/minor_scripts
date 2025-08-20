################################################################################
##  UniLasso vs Lasso  – home-court design only  (fast, with progress bar)
################################################################################
library(glmnet)
library(uniLasso)         # devtools::install_github("trevorhastie/uniLasso")
library(MASS)
library(pbapply)          # text progress bar

## --- simulate one home-court dataset ---------------------------------------
sim_home <- function(n = 300, p = 1000, rho = 0.8, k = 100) {
  Sigma <- rho ^ abs(outer(1:p, 1:p, "-"))
  X     <- mvrnorm(n, rep(0, p), Sigma)
  
  beta <- rep(0, p)
  idx  <- seq(1, 2*k, by = 2)          # every other column
  beta[idx] <- runif(k, .5, 2)
  
  sigma <- sd(X %*% beta) / 3.4        # target SNR ≈ 3.4
  y     <- X %*% beta + rnorm(n, 0, sigma)
  list(X = X, y = y, beta = beta, idx = idx)
}

## --- one replicate ----------------------------------------------------------
one_rep <- function(seed) {
  set.seed(seed)
  tr <- sim_home(); te <- sim_home()
  
  las <- cv.glmnet(tr$X, tr$y, nfolds = 10, intercept = TRUE, standardize = FALSE)
  uni <- cv.uniLasso(tr$X, tr$y, nfolds = 10)
  
  bl  <- as.numeric(coef(las, s = "lambda.min"))[-1]
  bu  <- as.numeric(coef(uni, s = "lambda.min"))[-1]
  
  mse <- function(b) mean((te$y - te$X %*% b)^2)
  sup <- function(b, t = 1e-6) which(abs(b) > t)
  
  c(mse_las = mse(bl), mse_uni = mse(bu),
    supp_las = length(sup(bl)), supp_uni = length(sup(bu)))
}

## --- run B replicates with a progress bar -----------------------------------
B <- 1
res <- pbapply::pbreplicate(B, one_rep(sample.int(1e6, 1)), simplify = "array")
out <- as.data.frame(t(res))

cat("\nHome-court summary (B = ", B, "):\n", sep = "")
summary_tbl <- cbind(
  metric = c("Test-MSE", "Support"),
  Lasso  = c(sprintf("%.2f ± %.2f",
                     mean(out$mse_las), sd(out$mse_las)/sqrt(B)),
             sprintf("%.1f", mean(out$supp_las))),
  UniLasso = c(sprintf("%.2f ± %.2f",
                       mean(out$mse_uni), sd(out$mse_uni)/sqrt(B)),
               sprintf("%.1f", mean(out$supp_uni)))
)
print(summary_tbl, row.names = FALSE, right = FALSE)
