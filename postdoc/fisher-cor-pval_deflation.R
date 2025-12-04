
cor_test_hac <- function(x, y, bw = NULL, prewhite = TRUE){
  stopifnot(length(x) == length(y))
  n <- length(x)
  df <- data.frame(y = as.numeric(y), x = as.numeric(x), t = 1:n)
  
  m  <- lm(y ~ x, data = df)
  
  # simple bandwidth rule-of-thumb if none supplied
  if (is.null(bw)) {
    # newey-west plug-in style; conservative if n is small
    bw <- floor(4 * (n / 100)^(2/9))
  }
  
  # HAC vcov (Bartlett kernel), optional prewhitening
  Vhac <- sandwich::NeweyWest(m, lag = bw, prewhite = prewhite, adjust = TRUE)
  tt   <- lmtest::coeftest(m, vcov. = Vhac)["x", ]
  
  # return t and p for H0: beta_x = 0  (equivalent to rho = 0)
  list(t = unname(tt[ "t value" ]), p = unname(tt[ "Pr(>|t|)" ]), bw = bw)
}

cor_test_effdf <- function(x, y, K = NULL){
  stopifnot(length(x) == length(y))
  n  <- length(x)
  r  <- cor(x, y)
  
  # choose truncation K if not supplied
  if (is.null(K)) K <- floor(n^(1/3))
  acx <- acf(x, lag.max = K, plot = FALSE)$acf[-1]
  acy <- acf(y, lag.max = K, plot = FALSE)$acf[-1]
  
  # plug-in EDF
  denom <- 1 + 2 * sum(acx * acy)
  neff  <- n / max(denom, 1e-8)  # guard against division by ~0
  
  # t-test using neff
  tstat <- r * sqrt((neff - 2) / max(1 - r^2, 1e-12))
  pval  <- 2 * pt(abs(tstat), df = max(neff - 2, 1), lower.tail = FALSE)
  
  list(r = r, t = unname(tstat), p = unname(pval), neff = neff, K = K)
}

cor_test_tts <- function(x, y, B = 2000, q = 0.1){
  stopifnot(length(x) == length(y))
  n   <- length(x)
  m   <- floor(q * n)        # truncation
  idx <- (m + 1):(n - m)     # usable interior
  r0  <- cor(x[idx], y[idx])
  
  # generate circular shifts; if B >= n use all unique shifts
  shifts <- if (B >= (n - 1)) 1:(n - 1) else sample(1:(n - 1), B, replace = FALSE)
  
  rperm <- numeric(length(shifts))
  for (i in seq_along(shifts)) {
    s    <- shifts[i]
    ys   <- c(y[(s + 1):n], y[1:s])  # circular shift
    rperm[i] <- cor(x[idx], ys[idx])
  }
  
  p <- (1 + sum(abs(rperm) >= abs(r0))) / (length(rperm) + 1)
  list(r = r0, p = p, shifts = length(shifts), q = q)
}

cor.sim <- function(r_betw = 0.0, 
                    r_within = 0.9,
                    p = 100,
                    n = 2,
                    permute_R = F,
                    corr_struct = c("block", "AR(1)")[2],
                    corr_test = c("sample", "effdf", "tts", "hac")[1]){
  
  #induce correlation within features
  if(corr_struct == "block"){
    #basic two block matrix
    Rsub <- diag(p/2) + r_within - diag(p/2) * r_within
    R <- diag(p)
    R[1:(p/2), 1:(p/2)] <- R[p/2+1:(p/2), p/2+1:(p/2)] <- Rsub  
  } else if(corr_struct == "AR(1)"){
    #AR(1) correlation structutre
    R <- r_within ^ abs(outer(1:p, 1:p, "-"))
  }
  
  if(permute_R){
    R_reord <- sample(1:p)
    R <- R[R_reord, R_reord]
  }
  
  #transform data within features
  L <- t(chol(R))
  x <- matrix(rnorm(p*2), p, n)
  x <- L %*% x
  
  #induce true correlation between features
  R2  <- diag(2) + r_betw - diag(2) * r_betw
  L2 <- t(chol(R2))
  x <- t(L2 %*% t(x))
  
  #test hypothesis
  if(corr_test == "sample"){
    out <- cor.test(x[,1], x[,2])
    return(c(sample_cor = out$estimate, pval = out$p.value))
  } else if(corr_test == "effdf"){
    out <- cor_test_effdf(x[,1], x[,2])
    return(c(sample_cor = out$r, pval = out$p))
  } else if(corr_test == "tts"){
    out <- cor_test_tts(x[,1], x[,2])
    return(c(sample_cor = out$r, pval = out$p))
  } else if(corr_test == "hac"){
    out <- cor_test_hac(x[,1], x[,2])
    return(c(sample_cor = out$t, pval = out$p))
  }
}

cor.sim <- function(r_betw = 0.0, 
                    r_within = 0.9,
                    p = 100,
                    n = 2,
                    permute_R = F,
                    corr_struct = c("block", "AR(1)")[2],
                    corr_test = c("sample")[1]){
  
  #induce correlation within features
  if(corr_struct == "block"){
    #basic two block matrix
    Rsub <- diag(p/2) + r_within - diag(p/2) * r_within
    R <- diag(p)
    R[1:(p/2), 1:(p/2)] <- R[p/2+1:(p/2), p/2+1:(p/2)] <- Rsub  
  } else if(corr_struct == "AR(1)"){
    #AR(1) correlation structutre
    R <- r_within ^ abs(outer(1:p, 1:p, "-"))
  }
  
  if(permute_R){
    R_reord <- sample(1:p)
    R <- R[R_reord, R_reord]
  }
  
  #transform data within features
  L <- t(chol(R))
  x <- matrix(rnorm(p*2), p, n)
  x <- L %*% x
  
  #induce true correlation between features
  R2  <- diag(2) + r_betw - diag(2) * r_betw
  L2 <- t(chol(R2))
  x <- t(L2 %*% t(x))
  
  #test hypothesis
  if(corr_test == "sample"){
    out <- cor.test(x[,1], x[,2])
    return(c(sample_cor = out$estimate, pval = out$p.value))
  }
}

res <- data.frame(t(replicate(n = 1E3, 
                              cor.sim(r_within = 0.7,
                                      r_betw = 0.0,
                                      p = 200,
                                      corr_struct = c("block", "AR(1)")[2],
                                      permute_R = F)
)))
par(mfrow = c(2,1), mar = c(4,3,1,1))
hist(res$sample_cor, breaks = 100, main = "sample correlation under null is technically unbiased")
abline(v = mean(res$sample_cor), col = 2, lwd = 3)
hist(res$pval, breaks = 100, main = "but p-values under the null are super deflated")

