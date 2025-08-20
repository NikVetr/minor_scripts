set.seed(1)

#specify simulation parameters
n_g1 <- 6
n_g2 <- 4
n <- n_g1 + n_g2
a <- c(3,5)
b <- c(2,4)
s <- c(1,1)
g <- c(rep(1, n_g1), rep(2, n_g2))
gf <- as.factor(g)
# g_s2z <- c(-1/2, 1/2)[g]
g_s2z <- c(-n_g1/n, n_g2/n)[g] #can also adjust these under imbalanced design if necessary
g_s2z <- c(-1/2, 1/2)[g] #can also adjust these under imbalanced design if necessary

#simulate data
x1 <- rnorm(n_g1)
x2 <- rnorm(n_g2)
x <- c(x1, x2)
xg_s2z <- x * g_s2z
e <- rnorm(n) * s[g]
y <- a[g] + b[g] * x + e

#fit models
fit_main <- lm(y ~ x)
coef_main <- summary(fit_main)$coefficients

fit_sep1 <- lm(y[g==1] ~ x[g==1])
coef_sep1 <- summary(fit_sep1)$coefficients

fit_sep2 <- lm(y[g==2] ~ x[g==2])
coef_sep2 <- summary(fit_sep2)$coefficients

fit_inter <- lm(y ~ x + gf + x*gf)
coef_inter <- summary(fit_inter)$coefficients

fit_inter_s2z <- lm(y ~ x + g_s2z + xg_s2z)
coef_inter_s2z <- summary(fit_inter_s2z)$coefficients

#recover models from each other
coef_sep_recon <- list(
  #these get group specific estimates of error_variance (sigma)
  axg = cbind(
    a_x_g_est = coef_sep1[1,1] - coef_sep2[1,1],
    a_x_g_se = sqrt(coef_sep1[1,2]^2 + coef_sep2[1,2]^2)
  ),
  bxg = cbind(
    b_x_g_est = coef_sep1[2,1] - coef_sep2[2,1],
    b_x_g_se = sqrt(coef_sep1[2,2]^2 + coef_sep2[2,2]^2)
  ), 
  a = cbind(
    a_est = (coef_sep1[1,1] + coef_sep2[1,1]) / 2,
    a_se = sqrt((coef_sep1[1,2]^2 + coef_sep2[1,2]^2) / 4)
  ), 
  b = cbind(
    b_est = (coef_sep1[2,1] + coef_sep2[2,1]) / 2,
    b_se = sqrt((coef_sep1[2,2]^2 + coef_sep2[2,2]^2) / 4)
  )
)

#inspect results
coef_main
coef_inter_s2z
coef_inter
coef_sep_recon

#get main effects from no-interaction model from dummy interaction fit
coef_main
avg_main <- (function(fit) {
  b <- coef(fit);          V <- vcov(fit);         df <- df.residual(fit)
  
  ## L-vectors that take the simple (½,½) average of the two groups
  L_int   <- setNames(numeric(length(b)), names(b))
  L_slope <- L_int
  L_int[c("(Intercept)",  "gf2")]  <- c(1,  0.5)
  L_slope[c("x",          "x:gf2")]<- c(1,  0.5)
  
  pick <- function(L) {
    est <- sum(L * b)
    se  <- sqrt(drop(t(L) %*% V %*% L))
    t   <- est / se
    p   <- 2 * pt(abs(t), df, lower.tail = FALSE)
    c(est = est, se = se, t = t, p = p)
  }
  
  rbind(intercept = pick(L_int),
        slope     = pick(L_slope))
})(fit_inter) 
avg_main
coef_inter_s2z