
panel.hist.ks <- function(x, ...) {
  x <- x[!is.na(x)]
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.8))
  h <- hist(x, plot = FALSE)
  if (length(h$counts) == 0) {
    return()
  }
  
  y <- h$counts
  y <- y / max(y)
  
  breaks <- h$breaks
  nB <- length(breaks)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", border = "black")
  if(any(x > 1 | x < 0)){
    x <- 10^x
  }
  ks <- suppressWarnings(ks.test(x, "punif"))
  txt <- paste0("ks-pval =\n", formatC(ks$p.value, format = "e", digits = 2))
  text(x = 0, y = 1.4, labels = txt, cex = 0.8, pos = 4)
}
panel.scatter.1to1 <- function(x, y, ...) {
  points(x, y, pch = 19, ...)
  abline(a = 0, b = 1, lty = 2, lwd = 2, col = 2)
}

softmax <- function(x) exp(x) / sum(exp(x))
simulate_null <- function(n, nrows, ncols, 
                          row_sd = NULL, col_sd = NULL,
                          guard_chi2 = T){
  if(is.null(row_sd)){
    row_sd <- rep(1, nrows)
  }
  if(is.null(col_sd)){
    col_sd <- rep(1, ncols)
  }
  row_probs <- softmax(rnorm(nrows) * row_sd)
  col_probs <- softmax(rnorm(ncols) * col_sd)
  cell_probs <- outer(row_probs, col_probs, FUN = "*")
  ct <- matrix(rmultinom(n = 1, size = n, prob = cell_probs), nrows, ncols)
  out <- c(fisher = fisher.test(ct, workspace = 2e8)$p.value,
           chi2 = chisq.test(ct, simulate.p.value = guard_chi2 & any(ct <= 5), 
                             B = 1E4)$p.value,
           min_cell = min(ct))
  return(out)
}


n <- 1E3
nrep <- 2E3
nrows <- 2
ncols <- 3
guard_chi2 <- T
row_sd <- rep(1, nrows)
col_sd <- rep(1, ncols)
row_sd <- c(1, rep(3, nrows-1))
col_sd <- c(1,rep(3, ncols-1))
pvals <- data.frame(do.call(rbind, replicate(n = nrep, 
                                             simulate_null(n, nrows, ncols, row_sd, col_sd, 
                                                           guard_chi2 = guard_chi2), 
                                             simplify = F)))
cols <- adjustcolor(c(4,1), 0.4)[(pvals$min_cell > 5) + 1]
pairs(log10(pvals[,c("fisher", "chi2")]), 
      diag.panel = panel.hist.ks, 
      lower.panel = panel.scatter.1to1, 
      upper.panel = panel.scatter.1to1, 
      col = cols)
legend("topleft", col = c(adjustcolor(c(1,4), 0.4), 2), pch = c(19, 19, NA), lty = c(NA, NA, 2),
       legend = c("n > 5", "n <= 5", "1-to-1 line"), xpd = NA, bty = "n" , horiz = T,
       cex = 0.8, x.intersp = 0.5)


#one off test for sohaib
ct <- matrix(c(46911, 5, 13,  1), 2, 2)
fisher.test(ct)$p.value
chisq.test(ct, simulate.p.value = T, B = 1E7)
Exact::exact.test(ct, model = "Binomial")
Exact::exact.test(ct, model = "Binomial", method = "csm")
ct <- matrix(c(46911, 5, 13,  1), 2, 2)
exact2x2::uncondExact2x2(x1 = ct[1,1],
                         n1 = sum(ct[1,]),
                         x2 = ct[2,1],
                         n2 = sum(ct[2,]), 
                         method = "score")
exact2x2::uncondExact2x2(x1 = ct[1,1],
                         n1 = sum(ct[1,]),
                         x2 = ct[2,1],
                         n2 = sum(ct[2,]))
exact2x2::uncondExact2x2(x1 = ct[1,1],
                         n1 = sum(ct[1,]),
                         x2 = ct[2,1],
                         n2 = sum(ct[2,]),
                         parmtype   = "difference",
                         nullparm   = 0,
                         alternative = "two.sided",
                         method     = "wald-pooled",
                         tsmethod   = "square")

#binomial Bayes
nMC <- 1E7
mean((rbeta(nMC, ct[1,1] + 1, ct[1,2] + 1) - 
  rbeta(nMC, ct[2,1] + 1, ct[2,2] + 1)) < 0)

#multinomial freq
row_prop <- rowSums(ct) / sum(ct)
col_prop <- colSums(ct) / sum(ct)
p_null <- as.vector(outer(row_prop, col_prop))
chisq_2x2 <- function(tab) {
  n <- sum(tab)
  r <- rowSums(tab)
  c <- colSums(tab)
  expected <- outer(r, c) / n
  nz <- expected > 0
  x2 <- sum((tab[nz] - expected[nz])^2 / expected[nz])
  return(x2)
}

TS_obs <- chisq_2x2(ct)
nMC <- 1E6
TS_sim <- replicate(nMC, chisq_2x2(matrix(rmultinom(1, size = N, prob = p_null), nrow = 2)))
mean(TS_sim >= TS_obs)
