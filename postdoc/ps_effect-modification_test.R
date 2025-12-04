# ---------- simulate ----------
nA <- 1E5; nB <- 1E5
ZA <- data.frame(age = rnorm(nA, 40, 8), sex = rbinom(nA, 1, 0.4))
ZB <- data.frame(age = rnorm(nB, 60, 8), sex = rbinom(nB, 1, 0.6))

XA <- rnorm(nA, 0.1*ZA$age + 0.5*ZA$sex, 1)
XB <- rnorm(nB, 0.1*ZB$age + 0.5*ZB$sex, 1)

betaX <- 0.5
YA <- betaX*XA + 0.02*ZA$age + 0.3*ZA$sex + 0.01*XA*ZA$age + rnorm(nA, 0, 1)
YB <- betaX*XB + 0.02*ZB$age + 0.3*ZB$sex + 0.01*XB*ZB$age + rnorm(nB, 0, 1)

A <- data.frame(S = 0, Y = YA, X = XA, age = ZA$age, sex = ZA$sex)
B <- data.frame(S = 1, Y = YB, X = XB, age = ZB$age, sex = ZB$sex)
dat <- rbind(A, B)

# ---------- helper: standardized mean difference ----------
smd <- function(x, w, S){
  m1 <- weighted.mean(x[S==1], w[S==1])
  m0 <- weighted.mean(x[S==0], w[S==0])
  v1 <- sum(w[S==1]*(x[S==1]-m1)^2)/sum(w[S==1])
  v0 <- sum(w[S==0]*(x[S==0]-m0)^2)/sum(w[S==0])
  (m1 - m0) / sqrt((v1+v0)/2)
}

# ---------- selection model & weights to map A → B on Z ----------
sel_fit <- glm(S ~ age + sex + I(age^2) + age:sex, data = dat, family = binomial())
pS1 <- predict(sel_fit, type = "response")
w_all <- rep(1, nrow(dat))
w_all[dat$S==0] <- pS1[dat$S==0] / (1 - pS1[dat$S==0])  # inverse odds for A; B has weight 1

cat("# debug balance: SMD(age)  before:", round(smd(dat$age, rep(1,nrow(dat)), dat$S), 3), "\n")
cat("# debug balance: SMD(age)  after :", round(smd(dat$age, w_all, dat$S), 3), "\n")
cat("# debug balance: SMD(sex)  before:", round(smd(dat$sex, rep(1,nrow(dat)), dat$S), 3), "\n")
cat("# debug balance: SMD(sex)  after :", round(smd(dat$sex, w_all, dat$S), 3), "\n")

# ---------- function: fit one spec and extract comparable 'slope of X' ----------
# spec_id:
#   1: no weights, model Y ~ X
#   2: weights (A→B), model Y ~ X
#   3: weights (A→B), model Y ~ X + age + sex
#   4: weights (A→B), model Y ~ X*age + X*sex + age + sex  (compute average derivative under B's Z)
get_slope_pair <- function(dat, w_all, spec_id){
  stopifnot(spec_id %in% 1:4)
  # choose weights
  w_use <- switch(as.character(spec_id),
                  "1" = rep(1, nrow(dat)),
                  "2" = w_all,
                  "3" = w_all,
                  "4" = w_all)
  # subset A and B
  idxA <- which(dat$S==0); idxB <- which(dat$S==1)
  dA <- dat[idxA, ]; dB <- dat[idxB, ]
  wA <- w_use[idxA]; wB <- w_use[idxB]
  
  # fit models
  fmla <- switch(as.character(spec_id),
                 "1" = Y ~ X,
                 "2" = Y ~ X,
                 "3" = Y ~ X + age + sex,
                 "4" = Y ~ X*age + X*sex + age + sex)
  
  fitA <- lm(fmla, data = dA, weights = wA)
  fitB <- lm(fmla, data = dB, weights = wB)
  
  # extract comparable slope
  if (spec_id %in% c(1,2,3)){
    # coefficient on X is the estimand
    bA <- unname(coef(fitA)["X"])
    bB <- unname(coef(fitB)["X"])
  } else {
    # spec 4: average derivative wrt X under B's Z
    bA <- {
      cf <- coef(fitA)
      bx  <- if ("X" %in% names(cf)) cf["X"] else 0
      bxA <- if ("X:age" %in% names(cf)) cf["X:age"] else 0
      bxS <- if ("X:sex" %in% names(cf)) cf["X:sex"] else 0
      Zage <- dB$age
      Zsex <- dB$sex
      mean(bx + bxA*Zage + bxS*Zsex)
    }
    bB <- {
      cf <- coef(fitB)
      bx  <- if ("X" %in% names(cf)) cf["X"] else 0
      bxA <- if ("X:age" %in% names(cf)) cf["X:age"] else 0
      bxS <- if ("X:sex" %in% names(cf)) cf["X:sex"] else 0
      Zage <- dB$age
      Zsex <- dB$sex
      mean(bx + bxA*Zage + bxS*Zsex)
    }
  }
  
  # debug prints
  cat(sprintf("# debug spec %d: slope_A = %.3f, slope_B = %.3f, diff (A-B) = %.3f\n",
              spec_id, bA, bB, bA - bB))
  list(A = bA, B = bB, diff = bA - bB,
       fitA = fitA, fitB = fitB)
}

# ---------- run all four specs ----------
res1 <- get_slope_pair(dat, w_all, 1)  # no weights, Y ~ X
res2 <- get_slope_pair(dat, w_all, 2)  # weights,   Y ~ X
res3 <- get_slope_pair(dat, w_all, 3)  # weights,   Y ~ X + age + sex
res4 <- get_slope_pair(dat, w_all, 4)  # weights,   Y ~ X*age + X*sex + age + sex (avg derivative under B's Z)

# ---------- concise results table ----------
tab <- rbind(
  "1) no wts; Y~X"                           = c(res1$A, res1$B, res1$diff),
  "2) wts;   Y~X"                            = c(res2$A, res2$B, res2$diff),
  "3) wts;   Y~X+age+sex"                    = c(res3$A, res3$B, res3$diff),
  "4) wts;   Y~X*age+X*sex+age+sex (avg dX)" = c(res4$A, res4$B, res4$diff)
)
colnames(tab) <- c("slope_A", "slope_B", "A_minus_B")
print(round(tab, 3))
