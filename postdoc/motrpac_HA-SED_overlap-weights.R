library(data.table)

#read in and integrate data
d <- fread("~/data/motrpac_precovid_clinical_obs.txt")
d_ha <- fread("~/data/motrpac_ha_clinical_obs.txt")
d_ha$modality <- c("ATHResist" = "RE", "ADUEndur" = "EE",
                   "ATHEndur" = "EE", "ADUControl" = "Control",
                   "ADUResist" = "RE")[d_ha$modality]
pid2modality <- setNames(d_ha$modality, d_ha$pid)
d$modality <- pid2modality[as.character(d$pid)]

#read in plotting functions
source("/Users/nikgvetr/repos/flexible-sum-to-zero/Nik/R/functions.R")
source("/Users/nikgvetr/repos/motrpac_ha/R/functions.R")

#####

#split data of interest up
var_index <- 5
target_variable <- c("age_years", 
                     "sex", 
                     "body_weight",
                     "bmi", 
                     "body_height_cm")[var_index]
label <- c("Age", "Sex", "Weight (kg)", "BMI", "Height (cm)")[var_index]
target_data <- d[[target_variable]]
target_split <- split(target_data, d$activity_level)
names(target_split) <- c("HA", "SED")
if(class(target_data) %in% c("numeric", "integer")){
  par(mfrow = c(1,2), mar = c(5,4,2,1))
  multihist(target_split, nbreaks = 20, 
            legend_location = "topright", 
            leg_loc_y_extra = 0.05, leg_loc_x_extra = -0.35, 
            main = label, xlab = label,
            no_border = T, top_border = T, color_border = T)
  qq_compare(xlist = target_split[2:1], qs = 1:99/100, 
             lwd = 3, main = paste0(label, " Quantiles"), cex.main = 1)
} else {
  par(mfrow = c(1,2), mar = c(3,4,2,1))
  stacked_barplot_list(target_split, main = paste0(paste0(rep(" ", 50), collapse = ""), label), use_proportion = T, leg_loc_x_extra = 100, xpd = NA)
  stacked_barplot_list(target_split, main = "", use_proportion = F, leg_loc_x_extra = -1, run_chisq = F)
}

#### calc overlap weights ####

# inputs: data.frame or data.table `d` with columns:
# activity_level (chr), sex (chr), age (num), bmi (num)

# - helpers -

kish_ess <- function(w) {
  # kish (1965/1992): (sum w)^2 / sum(w^2)
  s <- sum(w)
  if (any(w < 0) || s <= 0) stop("weights must be nonnegative with positive sum")
  (s * s) / sum(w * w)
}

normalize_weights <- function(w, target_sum) {
  # scale so sum(w) = target_sum (use n so ess is interpretable as 'out of n')
  if (sum(w) == 0) stop("sum of weights is zero; cannot normalize")
  w * (target_sum / sum(w))
}

# standardized mean difference (optional balance diagnostic; base-R only)
smd <- function(x, a, w = NULL) {
  # x numeric; a in {0,1}; w optional weights
  if (is.null(w)) w <- rep(1, length(a))
  m1 <- sum(w[a == 1] * x[a == 1]) / sum(w[a == 1])
  m0 <- sum(w[a == 0] * x[a == 0]) / sum(w[a == 0])
  v1 <- sum(w[a == 1] * (x[a == 1] - m1)^2) / sum(w[a == 1])
  v0 <- sum(w[a == 0] * (x[a == 0] - m0)^2) / sum(w[a == 0])
  (m1 - m0) / sqrt(0.5 * (v1 + v0))
}

# - prepare data 

df <- as.data.frame(d)
df$age <- df$age_years
keep <- c("activity_level", "sex", "age", "bmi")
df <- df[complete.cases(df[, keep]), keep]

# binary treatment: 1 = Adult Highly Active, 0 = Adult Sedentary
df$A <- as.integer(df$activity_level == "Adult Highly Active")
df <- df[df$activity_level %in% c("Adult Highly Active", "Adult Sedentary"), ]

# sex as factor (from your summary: levels Female/Male)
df$sex <- factor(df$sex)

message("n (complete cases in the two groups): ", nrow(df))
message("group counts:\n", paste(capture.output(print(table(df$activity_level))), collapse = "\n"))

# - propensity model (flexible but base-R) 
# natural splines for age/bmi to avoid strong linearity assumptions
fm <- A ~ sex + splines::ns(age, df = 4) + splines::ns(bmi, df = 4)

ps_fit <- glm(fm, data = df, family = binomial())
ps_hat <- predict(ps_fit, type = "response")

# clip tiny numerical extremes for stability (logistic won't give exact 0/1)
eps <- 1e-6
ps_hat <- pmin(pmax(ps_hat, eps), 1 - eps)

df$ps <- ps_hat

message(sprintf("propensity summary: min=%.3f, 1stQ=%.3f, med=%.3f, 3rdQ=%.3f, max=%.3f",
                min(df$ps), quantile(df$ps, 0.25), median(df$ps),
                quantile(df$ps, 0.75), max(df$ps)))
message(sprintf("share with ps<0.05: %.1f%%; ps>0.95: %.1f%%",
                100 * mean(df$ps < 0.05), 100 * mean(df$ps > 0.95)))

# - overlap weights -
# li–morgan–zaslavsky ow: w = 1-ps for treated, w = ps for control
w_ow <- ifelse(df$A == 1, 1 - df$ps, df$ps)

# normalize so sum(w) = n (useful for ESS interpretability and stable reporting)
w_ow_norm <- normalize_weights(w_ow, target_sum = nrow(df))

df$w_ow <- w_ow_norm

message(sprintf("overlap weights: min=%.4f, mean=%.4f, max=%.4f",
                min(df$w_ow), mean(df$w_ow), max(df$w_ow)))

# - kish ESS 
ess_all <- kish_ess(df$w_ow)
ess_t1  <- kish_ess(df$w_ow[df$A == 1])
ess_t0  <- kish_ess(df$w_ow[df$A == 0])

message(sprintf("kish ESS (overall): %.1f out of %d", ess_all, nrow(df)))
message(sprintf("kish ESS (Adult Highly Active): %.1f out of %d",
                ess_t1, sum(df$A == 1)))
message(sprintf("kish ESS (Adult Sedentary): %.1f out of %d",
                ess_t0, sum(df$A == 0)))

# - optional: quick balance check before/after 
message(sprintf("SMD(age) unweighted=%.3f, OW=%.3f",
                smd(df$age, df$A),
                smd(df$age, df$A, w = df$w_ow)))
message(sprintf("SMD(bmi) unweighted=%.3f, OW=%.3f",
                smd(df$bmi, df$A),
                smd(df$bmi, df$A, w = df$w_ow)))
# for sex (binary), code to 0/1 then use smd
sex01 <- as.integer(df$sex == levels(df$sex)[1])
message(sprintf("SMD(sex==%s) unweighted=%.3f, OW=%.3f",
                levels(df$sex)[1],
                smd(sex01, df$A),
                smd(sex01, df$A, w = df$w_ow)))

# - minimal unit tests 
# 1) ps in (0,1), weights bounded and nonnegative
stopifnot(all(df$ps > 0 & df$ps < 1))
stopifnot(all(df$w_ow >= 0))

# 2) normalization sanity
stopifnot(abs(sum(df$w_ow) - nrow(df)) < 1e-8)

# 3) kish ESS: equal weights -> ESS = n
stopifnot(all.equal(kish_ess(rep(1, 10)), 10))

# 4) OW bounded by 1 (since ps in (0,1))
stopifnot(max(w_ow) <= 1 && min(w_ow) >= 0)

# - results you likely want to keep 
# df$w_ow   # normalized overlap weights (sum = n)
# ess_all; ess_t1; ess_t0


#### interpret ####
# assumes you already built:
# df with columns: A (0/1), sex (factor), age (numeric), bmi (numeric), ps (propensity), w_ow (normalized overlap weights)

#  plotting helpers (base R only) 
plot_weight_hist <- function(w, main = "overlap weights (normalized)", bins = 30) {
  hist(w, breaks = bins, main = main, xlab = "weight", col = "gray", border = "white")
  abline(v = median(w), lty = 2)
}

plot_ps_densities <- function(ps, A, main = "propensity score densities by group") {
  d1 <- density(ps[A == 1])
  d0 <- density(ps[A == 0])
  rng <- range(c(d1$y, d0$y))
  plot(d1, main = main, xlab = "propensity score", ylim = rng)
  lines(d0, lty = 2)
  legend("topright", c("A=1", "A=0"), lty = c(1,2), bty = "n")
}

plot_ps_qq <- function(ps, A, main = "propensity score Q–Q (A=1 vs A=0)") {
  q <- seq(0.01, 0.99, by = 0.01)
  q1 <- quantile(ps[A == 1], q); q0 <- quantile(ps[A == 0], q)
  plot(q0, q1, xlab = "A=0 quantiles", ylab = "A=1 quantiles", main = main)
  abline(0,1,lty=2)
}

#  common support + share outside 
common_support <- function(ps, A) {
  min1 <- min(ps[A==1]); max1 <- max(ps[A==1])
  min0 <- min(ps[A==0]); max0 <- max(ps[A==0])
  lower <- max(min1, min0); upper <- min(max1, max0)
  outside1 <- mean(ps[A==1] < lower | ps[A==1] > upper)
  outside0 <- mean(ps[A==0] < lower | ps[A==0] > upper)
  list(bounds = c(lower = lower, upper = upper),
       frac_outside = c(A1 = outside1, A0 = outside0))
}

#  SMD table, including interactions / nonlinearities 
smd <- function(x, a, w = NULL) {
  if (is.null(w)) w <- rep(1, length(a))
  m1 <- sum(w[a == 1] * x[a == 1]) / sum(w[a == 1])
  m0 <- sum(w[a == 0] * x[a == 0]) / sum(w[a == 0])
  v1 <- sum(w[a == 1] * (x[a == 1] - m1)^2) / sum(w[a == 1])
  v0 <- sum(w[a == 0] * (x[a == 0] - m0)^2) / sum(w[a == 0])
  (m1 - m0) / sqrt(0.5 * (v1 + v0))
}

smd_table <- function(df, w = NULL) {
  # base set + some interactions/nonlinearities to probe multivariate structure
  sex01 <- as.integer(df$sex == levels(df$sex)[1])
  X <- data.frame(
    age = df$age,
    bmi = df$bmi,
    sex01 = sex01,
    age2 = df$age^2,
    bmi2 = df$bmi^2,
    age_bmi = df$age * df$bmi,
    age_sex = df$age * sex01,
    bmi_sex = df$bmi * sex01
  )
  nm <- names(X)
  out <- data.frame(var = nm,
                    SMD_unw = sapply(X, function(z) smd(z, df$A)),
                    SMD_OW  = sapply(X, function(z) smd(z, df$A, w = w)),
                    row.names = NULL)
  out
}

#  multivariate overlap via nearest-opposite Mahalanobis 
nearest_opp_mahal <- function(df, vars = c("age","bmi","sex")) {
  # code sex to 0/1 for distance
  sex01 <- as.integer(df$sex == levels(df$sex)[1])
  X <- cbind(df$age, df$bmi, sex01)
  colnames(X) <- c("age","bmi","sex01")
  S <- cov(X)  # pooled covariance
  invS <- solve(S)
  
  # function to get min Mahalanobis distance from x to rows of Y
  min_mahal <- function(x, Y) {
    # vectorized Mahalanobis to all rows in Y
    dif <- t(t(Y) - x)
    d2 <- rowSums((dif %*% invS) * dif)
    sqrt(min(d2))
  }
  
  X1 <- X[df$A==1,,drop=FALSE]
  X0 <- X[df$A==0,,drop=FALSE]
  
  d1 <- apply(X1, 1, min_mahal, Y = X0)
  d0 <- apply(X0, 1, min_mahal, Y = X1)
  
  list(d_to_opp_A1 = d1, d_to_opp_A0 = d0)
}

#  compact “is this good or bad?” summary 
diagnostic_summary <- function(df) {
  cs <- common_support(df$ps, df$A)
  cat(sprintf("common support bounds: [%.3f, %.3f]\n", cs$bounds[1], cs$bounds[2]))
  cat(sprintf("fraction outside bounds: A=1: %.1f%%, A=0: %.1f%%\n",
              100*cs$frac_outside[1], 100*cs$frac_outside[2]))
  
  smd_tab <- smd_table(df, w = df$w_ow)
  flag <- abs(smd_tab$SMD_OW) > 0.1
  cat(sprintf("SMD (|.|>0.1 flagged): %d/%d flagged after OW\n",
              sum(flag), nrow(smd_tab)))
  
  nn <- nearest_opp_mahal(df)
  q1 <- quantile(nn$d_to_opp_A1, c(.5,.9,.99))
  q0 <- quantile(nn$d_to_opp_A0, c(.5,.9,.99))
  cat(sprintf("nearest opposite-group Mahalanobis distance (A=1) med/90th/99th: %.2f / %.2f / %.2f\n",
              q1[1], q1[2], q1[3]))
  cat(sprintf("nearest opposite-group Mahalanobis distance (A=0) med/90th/99th: %.2f / %.2f / %.2f\n",
              q0[1], q0[2], q0[3]))
  
  invisible(list(common_support = cs, smd = smd_tab, nn = nn))
}

#  run diagnostics 
par(mfrow = c(2,2))
plot_weight_hist(df$w_ow)
plot_ps_densities(df$ps, df$A)
plot_ps_qq(df$ps, df$A)
plot(df$age, df$bmi, col = ifelse(df$A==1, "black", "gray"),
     xlab = "age", ylab = "bmi",
     main = "age × bmi by group (A=1 black, A=0 gray)")
par(mfrow = c(1,1))

res_diag <- diagnostic_summary(df)

#  unit tests (concise) 
stopifnot(is.numeric(res_diag$common_support$bounds),
          res_diag$common_support$bounds[1] <= res_diag$common_support$bounds[2])
# SMDs should not get worse after OW for modeled terms (heuristic)
stopifnot(mean(abs(res_diag$smd$SMD_OW)) <= mean(abs(res_diag$smd$SMD_unw)) + 1e-6)
# Mahalanobis distances must be nonnegative
stopifnot(all(res_diag$nn$d_to_opp_A1 >= 0),
          all(res_diag$nn$d_to_opp_A0 >= 0))


#### more plots ####

# assumes you already have:
# df with A (0/1), sex (factor), age, bmi, ps in (0,1), w_ow (normalized)
# ess_all, ess_t1, ess_t0 already computed (Kish ESS overall, HA, SED)
# nearest_opp_mahal() from earlier message (returns d_to_opp_A1, d_to_opp_A0)

# - colors & labels -
cols <- adjustcolor(rainbow(2), 0.5)  # semi-transparent
col_HA  <- cols[1]
col_SED <- cols[2]
lab_HA  <- "HA"
lab_SED <- "SED"

# - 1) Kish ESS plot -
N_overall <- nrow(df)
N_HA  <- sum(df$A == 1)
N_SED <- sum(df$A == 0)

ESS <- c(Overall = ess_all, HA = ess_t1, SED = ess_t0)
N   <- c(Overall = N_overall, HA = N_HA, SED = N_SED)
pct <- round(100 * ESS / N, 1)

par(mfrow = c(1,2))  # will also use the second slot for the Mahalanobis plot later

bp <- barplot(ESS,
              col = c("gray70", col_HA, col_SED),
              border = "gray30",
              ylim = c(0, max(N) * 1.15),
              main = "Kish ESS (overlap weights)",
              ylab = "Effective sample size",
              las = 1)
# annotate with "ESS (pct% of N=..)"
text(x = bp, y = ESS,
     labels = paste0(sprintf("%.1f", ESS), " (", pct, "% of N=", N, ")"),
     pos = 3, cex = 0.8)
# show the raw N as faint reference points
points(bp, N, pch = 4, col = "gray20")
legend("topright", c("ESS", "Raw N"),
       fill = c("gray70", NA), border = c("gray30", NA), pch = c(NA,4),
       bty = "n", inset = 0.02)

# - 2) Mahalanobis distance histogram (nearest opposite-group) -
# (recompute if needed)
# builds covariance matrices to use inside Mahalanobis distance
# cov_type: "overall", "pooled_within", "overall_OW", "pooled_within_OW"
build_cov <- function(df, cov_type = "overall") {
  sex01 <- as.integer(df$sex == levels(df$sex)[1])
  X <- cbind(df$age, df$bmi, sex01)
  colnames(X) <- c("age","bmi","sex01")
  
  wcov <- function(X, w = NULL) {
    if (is.null(w)) return(cov(X))
    w <- as.numeric(w); w <- w / sum(w)
    mu <- colSums(X * w)
    xc <- sweep(X, 2, mu, "-")
    (t(xc) %*% (xc * w))  # probability-weighted covariance
  }
  
  if (cov_type == "overall") {
    S <- cov(X)
  } else if (cov_type == "overall_OW") {
    S <- wcov(X, df$w_ow)
  } else if (cov_type == "pooled_within") {
    X1 <- X[df$A==1,,drop=FALSE]; X0 <- X[df$A==0,,drop=FALSE]
    n1 <- nrow(X1); n0 <- nrow(X0)
    S1 <- cov(X1); S0 <- cov(X0)
    S <- ((n1-1)*S1 + (n0-1)*S0) / (n1 + n0 - 2)
  } else if (cov_type == "pooled_within_OW") {
    X1 <- X[df$A==1,,drop=FALSE]; X0 <- X[df$A==0,,drop=FALSE]
    w1 <- df$w_ow[df$A==1]; w0 <- df$w_ow[df$A==0]
    S1 <- wcov(X1, w1); S0 <- wcov(X0, w0)
    S <- (sum(w1) * S1 + sum(w0) * S0) / (sum(w1) + sum(w0))
  } else stop("unknown cov_type")
  
  # tiny ridge in case S is near-singular
  if (det(S) <= .Machine$double.eps) {
    lam <- 1e-6
    S <- S + diag(lam, ncol(S))
  }
  S
}

nearest_opp_mahal <- function(df, cov_type = "pooled_within") {
  S <- build_cov(df, cov_type)
  invS <- solve(S)
  sex01 <- as.integer(df$sex == levels(df$sex)[1])
  X <- cbind(df$age, df$bmi, sex01)
  
  min_mahal <- function(x, Y) {
    dif <- t(t(Y) - x)
    d2 <- rowSums((dif %*% invS) * dif)
    sqrt(min(d2))
  }
  
  X1 <- X[df$A==1,,drop=FALSE]; X0 <- X[df$A==0,,drop=FALSE]
  d1 <- apply(X1, 1, min_mahal, Y = X0)  # HA to nearest SED
  d0 <- apply(X0, 1, min_mahal, Y = X1)  # SED to nearest HA
  list(d_to_opp_HA = d1, d_to_opp_SED = d0, cov_type = cov_type, S = S)
}

nn <- nearest_opp_mahal(df)

# overlayed histograms with matched breaks
brks <- pretty(c(nn$d_to_opp_HA, nn$d_to_opp_SED), n = 30)
h0 <- hist(nn$d_to_opp_SED, breaks = brks, plot = FALSE)  # SED->HA
h1 <- hist(nn$d_to_opp_HA, breaks = brks, plot = FALSE)  # HA->SED

plot(h0, col = col_SED, border = NA, main = "Nearest-opposite Mahalanobis",
     xlab = "distance")
plot(h1, col = col_HA,  border = NA, add = TRUE)
abline(v = median(nn$d_to_opp_HA), lty = 2, 
       col = adjustcolor(col_HA, 1), lwd = 3)
abline(v = median(nn$d_to_opp_SED), lty = 2, 
       col = adjustcolor(col_SED, 1), lwd = 3)
legend("topright",
       legend = c(paste0(lab_SED, " → ", lab_HA),
                  paste0(lab_HA,  " → ", lab_SED),
                  "medians"),
       fill = c(col_SED, col_HA, NA), border = c(NA, NA, NA),
       lty = c(NA, NA, 2), bty = "n")
#####
# assumes df has A (0/1), ps, w_ow, age, bmi, sex
# by construction: A==1 => HA (highly active), A==0 => SED (sedentary)


par(mfrow = c(1,2))

# propensity score densities by group (HA vs SED)
d_HA  <- density(df$ps[df$A==1], from = 0)
d_SED <- density(df$ps[df$A==0], from = 0)
d_HA$y[1] <- d_SED$y[1] <- 0
ylim <- range(c(d_HA$y, d_SED$y))

plot(d_SED, main = "propensity score densities by group", xlab = "propensity score",
     ylim = ylim, type = "n")  # empty plot, we'll add polygons manually
# add filled polygons for each group
polygon(d_SED$x, d_SED$y, col = adjustcolor(col_SED, alpha.f = 0.5), border = col_SED)
polygon(d_HA$x,  d_HA$y,  col = adjustcolor(col_HA, alpha.f = 0.5), border = col_HA)
legend("topright", legend = c(lab_SED, lab_HA), fill = c(col_SED, col_HA), bty = "n")


# 1) overlap weights (normalized)
hist(df$w_ow,
     breaks = 30, col = "gray80", border = "white",
     main = "overlap weights (normalized)", xlab = "weight")
abline(v = median(df$w_ow), lty = 2)

par(mfrow = c(1,1))

#### visualize overlap dist ####

# --- target weights for the overlap population (no group label) ------------
w_target <- df$ps * (1 - df$ps)   # proportional to e(x){1-e(x)}
stopifnot(all(w_target >= 0), sum(w_target) > 0)

# --- sample from the overlap population ------------------------------------
N_samp <- 5000
idx <- sample.int(nrow(df), size = N_samp, replace = TRUE, prob = w_target)
ovl <- df[idx, c("age","bmi","sex","A")]   # A retained only for color in plots; not part of target

# --- sanity check: weighted vs sampled moments -----------------------------
w_target_n <- w_target / sum(w_target)
w_mean_age <- sum(w_target_n * df$age)
w_mean_bmi <- sum(w_target_n * df$bmi)
message(sprintf("weighted means (age,bmi) under e(1-e): (%.2f, %.2f)", w_mean_age, w_mean_bmi))
message(sprintf("sampled means (age,bmi) from overlap sample: (%.2f, %.2f)",
                mean(ovl$age), mean(ovl$bmi)))

# --- simple visualizations of the overlap population -----------------------
cols <- adjustcolor(rainbow(2), 0.5); col_HA <- cols[1]; col_SED <- cols[2]
lab_HA <- "HA"; lab_SED <- "SED"

target_variables <- c("age_years", 
                      "sex", 
                      "body_weight",
                      "bmi", 
                      "body_height_cm")
labels <- c("Age", "Sex", "Weight (kg)", "BMI", "Height (cm)")
par(mfcol = c(3,3))

# 1) overlap-population age and BMI
for(var_index in c(1,4)){
  
  target_variable <- target_variables[var_index]
  target_data <- d[[target_variable]]
  label <- labels[var_index]
  target_split <- split(target_data, d$activity_level)
  names(target_split) <- c("HA", "SED")
  par(mar = c(5,4,2,1))
  
  #split pops
  multihist(target_split, nbreaks = 20, 
            legend_location = "topright", 
            leg_loc_y_extra = 0.05, leg_loc_x_extra = -0.275, 
            main = label, xlab = label,
            no_border = T, top_border = T, color_border = T)
  
  #pooled pop
  multihist(list(pooled = target_data), nbreaks = 20, 
            legend_location = "topright", 
            leg_loc_y_extra = 0.05, leg_loc_x_extra = -0.275, 
            main = "", xlab = label, run_ks.test = F, 
            plot_cols = adjustcolor("purple", 0.4),
            no_border = T, top_border = T, color_border = T)
  
  #overlap pop
  ovl_data <- ovl[,tolower(label)]
  multihist(list(overlap = ovl_data), nbreaks = 20, 
            legend_location = "topright", 
            leg_loc_y_extra = 0.05, leg_loc_x_extra = -0.275, 
            main = "", xlab = label, run_ks.test = F, 
            plot_cols = adjustcolor("violetred", 0.4),
            no_border = T, top_border = T, color_border = T)
}

# sex
var_index <- 2
target_variable <- target_variables[var_index]
target_data <- d[[target_variable]]
label <- labels[var_index]
target_split <- split(target_data, d$activity_level)
names(target_split) <- c("HA", "SED")

#split
stacked_barplot_list(target_split, main = paste0(paste0(rep(" ", 50), collapse = ""), label), 
                     use_proportion = T, leg_loc_x_extra = -1, xpd = NA, 
                     run_chisq = F, names_cex = 1, main_off.x = -0.05, 
                     leg_loc_y_extra = 0.2, legend_args = list(horiz = T))
abline(h = 0.5, col = adjustcolor(1, 0.4), lty = 2)

#pooled
stacked_barplot_list(list(pooled = target_data), main = paste0(paste0(rep(" ", 50), collapse = ""), label), 
                     use_proportion = T, leg_loc_x_extra = -1, xpd = NA, 
                     run_chisq = F, names_cex = 1, main_off.x = -0.05, 
                     leg_loc_y_extra = 0.2, legend_args = list(horiz = T))
abline(h = 0.5, col = adjustcolor(1, 0.4), lty = 2)

#overlap
ovl_data <- ovl[,tolower(label)]
stacked_barplot_list(list(overlap = ovl_data), main = paste0(paste0(rep(" ", 50), collapse = ""), label), 
                     use_proportion = T, leg_loc_x_extra = -1, xpd = NA, 
                     run_chisq = F, names_cex = 1, main_off.x = -0.05, 
                     leg_loc_y_extra = 0.2, legend_args = list(horiz = T))
abline(h = 0.5, col = adjustcolor(1, 0.4), lty = 2)

#### end ####
