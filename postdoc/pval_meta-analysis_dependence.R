source("~/scripts/minor_scripts/postdoc/scatter_hist.R")
source("~/scripts/minor_scripts/postdoc/pval_meta-analysis_dependence_functions.R")

#### basic analysis ####
n_rep <- 1E2
comb_methods <- list(
  fisher   = combine_fisher,
  harmonic = combine_harmonic,
  cauchy   = combine_cauchy,
  tippett     = combine_tippett,
  simes = combine_simes,
  stouffer = combine_stouffer,
  brown_cs = combine_brownCS
)

# example: actually run it
base_vals <- list(
  n_pval = 2^9,
  n_rep = n_rep,
  mu_alt = 1,
  prop_alt = 0.1,
  rho_cs = 0.75,
  rho_ar1 = 0.95,
  alpha = 0.1
)
results <- run_full_study(
  B = base_vals$n_rep,   # number of simulation replicates
  m = base_vals$n_pval,      # dimension
  mu_alt = base_vals$mu_alt,
  prop_alt = base_vals$prop_alt,
  rho_cs = base_vals$rho_cs,
  rho_ar1 = base_vals$rho_ar1,
  alpha = base_vals$alpha,
  comb_methods = comb_methods
)

#### plotting scatter-hists ####
source("~/scripts/minor_scripts/postdoc/scatter_hist.R")
plot_as_heatmap <- F
figures_dir <- "~/pval_meta-anal/"
subfigures_dir <- paste0(figures_dir, "subfigures/")
figure_paths <- c()

for(use_log10 in c(F,T)){
  
  
  depstructs <- c("indep", "cs", "ar1")
  conds <- c("null", "alt")
  
  figure_path <- paste0(subfigures_dir, "/meta-analysis_scatter-hist... ", 
                        paste0("n_pval = ", base_vals$n_pval, ", "),
                        paste0("mu_alt = ", base_vals$mu_alt, ", "),
                        paste0("prop_alt = ", base_vals$prop_alt, ", "),
                        paste0("rho_cs = ", base_vals$rho_cs, ", "),
                        paste0("rho_ar1 = ", base_vals$rho_ar1, ", "),
                        paste0("alpha = ", base_vals$alpha, ", "),
                        ifelse(use_log10, "log10transformed", "natural scale"),
                        ".pdf")
  figure_paths <- c(figure_paths, figure_path)
  cairo_pdf(filename = figure_path, width = 12, height = 8)
  par(mfcol = c(2,3), mar = c(3,3,0,0), oma = c(0,6,10,0))
  
  for(depstruct in depstructs){
    for(cond in conds){
      harmonic <- results[[depstruct]][[cond]]$combined_p$harmonic
      cauchy <-  results[[depstruct]][[cond]]$combined_p$cauchy
      alpha_thresh <- c(base_vals$alpha, base_vals$alpha)
      xbreaks <- ybreaks <- 0:10/10
      xlab <- "harmonic mean p-value"
      ylab <- "cauchy method p-value"
      if(use_log10){
        harmonic <- log10(harmonic)
        cauchy <-  log10(cauchy)
        alpha_thresh <- log10(alpha_thresh)
        xbreaks <- ybreaks <- pretty(c(harmonic, cauchy), n = 10)
        xlab <- "harmonic mean log₁₀(p-value)"
        ylab <- "cauchy method log₁₀(p-value)"
      }
      
      scatter_hist(x = harmonic, y = cauchy, col.hpt = adjustcolor(2,0.9),
                   shade_right_tail = F,
                   highlight_pt = alpha_thresh, shade_hist_tails = T, plot_pt = F,
                   target.ticks = 10, xbreaks = xbreaks, ybreaks = ybreaks,
                   xlab = xlab, ylab = ylab, asp = 1, 
                   add_one_to_one_line = T, cex.lab = 1, equal_axes = T,
                   plot_as_heatmap = plot_as_heatmap, col.pts = NULL, kde_heatmap = F)
      
      #check tl corner entries -- not meant to be run inside function
      check_tl <- F
      if(check_tl){
        tl_hits <- harmonic < 0.1 & cauchy > 0.9
        ntl_hits <- harmonic < 0.1 & cauchy < 0.1
        # ntl_hits <- !tl_hits
        tl_corner_ps <- results$indep$alt$settings$p_mat[tl_hits,]
        ntl_corner_ps <- results$indep$alt$settings$p_mat[ntl_hits,]
        
        mean(apply(tl_corner_ps, 1, var))
        mean(apply(ntl_corner_ps, 1, var))
        
        mean(apply(tl_corner_ps, 1, mean))
        mean(apply(ntl_corner_ps, 1, mean))
        
        #try minimum p value
        num_mins <- c(1, 2, 4, 8, 16)
        for(num_min in num_mins){
          stat_tlps <- apply(tl_corner_ps, 1, function(tmp) mean(sort(tmp)[1:num_min]))
          stat_ntlps <- apply(ntl_corner_ps, 1, function(tmp) mean(sort(tmp)[1:num_min]))
          
          breaks <- seq(min(c(stat_tlps, stat_ntlps)),
                        max(c(stat_tlps, stat_ntlps)),
                        length.out = 51)
          h_tlps <- hist(stat_tlps, breaks = breaks, plot = F)
          h_ntlps <- hist(stat_ntlps, breaks = breaks, plot = F)
          hist(stat_tlps, freq = F, breaks = breaks, 
               col = adjustcolor("orange", 0.5),
               xlim = range(c(stat_tlps, stat_ntlps)),
               ylim = range(c(h_tlps$density, h_ntlps$density)),
               xlab = ifelse(num_min == 1, "smallest p-value", 
                             paste0("mean of smallest ", num_min, " p-values")),
               main = "")
          hist(stat_ntlps, breaks = breaks,
               freq = F, add = T, col = adjustcolor("blue", 0.5))
          legend("topright",
                 legend = c("harmonic < 0.1 & cauchy > 0.9",
                            "harmonic < 0.1 & cauchy < 0.1"),
                 fill   = c(adjustcolor("orange", 0.5),
                            adjustcolor("blue",   0.5)),
                 border = NA,
                 bty    = "n",
                 cex    = 0.9)
          
        }
        
        
      }
      
      if(depstruct == depstructs[1]){
        horiz_labs <- c("null" = paste0("Null Hypothesis\n(all μ = 0)"),
                        "alt" = paste0("Alternative Hypothesis\n(μ* = ", base_vals$mu_alt, ", Pr(μ*) = ", base_vals$prop_alt, ")"))
        mtext(paste0(horiz_labs[cond]), outer = F, xpd = NA, side = 2, line = 4, cex = 1.5, font = 2)
      }
      if(cond == conds[1]){
        vert_labs <- c("indep" = paste0("Independent\nP-Values"),
                       "cs" = paste0("Compound Symmetry\n(rᵢⱼ = ", base_vals$rho_cs, ")"),
                       "ar1" = paste0("AR(1)\n(rᵢⱼ = ", base_vals$rho_ar1, ")"))
        mtext(paste0(vert_labs[depstruct]), outer = F, xpd = NA, side = 3, line = 0, cex = 1.5, font = 2)
      }
    }
  }
  mtext(paste0("Comparing Aggregation Methods for ", base_vals$n_pval," P-Values Generated Using Φ⁻¹(MVN) at α = ", base_vals$alpha), 
        outer = T, xpd = NA, side = 3, line = 6, cex = 1.65)
  
  
  dev.off()
  
  magick::image_write(magick::image_read_pdf(figure_path, density = 300), 
                      path = gsub("\\.pdf$", ".png", figure_path), format = "png")
  
}

all_figures_path <- paste0(figures_dir, "meta-analysis_scatter-hist... ", 
                           paste0("n_pval = ", base_vals$n_pval, ", "),
                           paste0("mu_alt = ", base_vals$mu_alt, ", "),
                           paste0("prop_alt = ", base_vals$prop_alt, ", "),
                           paste0("rho_cs = ", base_vals$rho_cs, ", "),
                           paste0("rho_ar1 = ", base_vals$rho_ar1, ", "),
                           paste0("alpha = ", base_vals$alpha),
                           ".pdf")
pdftools::pdf_combine(input = figure_paths, output = all_figures_path)

#### sweep analysis ####

#base values
base_vals <- list(
  n_pval = 2^9,
  n_rep = n_rep,
  mu_alt = 1,
  prop_alt = 0.1,
  rho_cs = 0.75,
  rho_ar1 = 0.95,
  alpha = 0.1
)

#ranges of values to check
range_vals <- list(
  n_pval = 2^(1:10),
  mu_alt = (0:10)/4,
  prop_alt = c((0.5)^(0:10), 0),
  rho = c(0, 1 - 0.5^(1:8))
)

sweep_results <- run_all_sweeps(
  base_vals  = base_vals,
  range_vals = range_vals,
  conf = 0.90,
  comb_methods = comb_methods,
  verbose = TRUE
)

#### plot sweep results ####
names(comb_methods)
methods_inds <- c(2,3)
all_methods <- names(comb_methods)
all_cols <- RColorBrewer::brewer.pal(8, "Dark2")
methods <- all_methods[methods_inds]
cols <- all_cols[methods_inds]
figure_paths <- c()
for(sp in names(range_vals)){
  
  figure_path <- paste0(subfigures_dir, "/meta-analysis_coverage-at-alpha-=",
                        base_vals$alpha, "_", paste0(methods, collapse = "-"), "_", sp, ".pdf")
  figure_paths <- c(figure_paths, figure_path)
  cairo_pdf(filename = figure_path, width = 10, height = 5)

  plot_sweep_param(
    all_sweep_results = sweep_results,
    methods = methods,
    sweep_param       = sp,
    base_vals         = base_vals,
    cols = cols
  )
  
  dev.off()
  
  magick::image_write(magick::image_read_pdf(figure_path, density = 300), 
                      path = gsub("\\.pdf$", ".png", figure_path), format = "png")
  
}

all_figures_path <- paste0(figures_dir, "/meta-analysis_coverage-at-alpha-=",
                           base_vals$alpha, ".pdf")
pdftools::pdf_combine(input = figure_paths, output = all_figures_path)


#### plots for presentation ####
source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")
# identity, cs, and ar(1) corrmats
p <- 20
rho_cs <- 0.5
rho_ar1 <- 0.8
R_titles <- list(
  iden = bquote("Identity:" ~ r[ij] == 0 ~~ "(" * i != j * ")"),
  cs = bquote(
    "Compound symmetry:" ~ 
      r[ij] == rho ~~ "(" * i != j * ")" ~~ " =" ~~ .(rho_cs)
  ),
  ar1 = bquote(
    "AR(1):" ~ 
      r[ij] == rho^{ group("|", i - j, "|") } ~~ 
      " =" ~~ .(rho_ar1)^{ group("|", i - j, "|") }
  )
)

Rs <- list(
  iden = diag(1, nrow = p, ncol = p),
  cs = {
    M <- matrix(rho_cs, nrow = p, ncol = p)
    diag(M) <- 1
    M
  },
  ar1 = {
    idx <- seq_len(p)
    D <- abs(outer(idx, idx, "-"))
    rho_ar1^D
  }
)

#add sign-flipped version
flip_vec <- rep(c(1, -1), length.out = p)
S <- diag(flip_vec)
Rs$altflip_ar1 <- S %*% Rs$ar1 %*% S
R_titles$altflip_ar1 <- bquote(
  "" ~
    r[ij] == .(rho_ar1)^{ group("|", i - j, "|") } ~
    "*" ~ s[i] ~ s[j] ~
    "," * "\n" ~
    s[k] %in% "{+1,-1}" ~
    "(" * s[k] == (-1)^{k-1} * ")"
)

for(Rn in names(R_titles)){
  figure_path <- paste0(subfigures_dir, "/presentations_matrix-", Rn, ".pdf")
  cairo_pdf(filename = figure_path, width = 5, height = 5)
  par(mar = c(1,1,6,1))
  my_heatmap(Rs[[Rn]], diag_matters = T, asp = 1, plot_diagonal_labels = F, 
             legend_title = expression(r[ij]), plot_labels = F, 
             plot_numbers = T, plot_guiding_lines = F, mat_size_rule = "abs", 
             number_col = "white", outer_border = T, main = R_titles[[Rn]])
  dev.off()
  magick::image_write(magick::image_read_pdf(figure_path, density = 300), 
                      path = gsub("\\.pdf$", ".png", figure_path), format = "png")
}


#some extra tiny plots
pv <- rbeta(1E4, 1, 1.05)
hist(pv, xlab = "Uniform(0,1) p-values", main = "", freq = F, col = "white", border = "white")
rect(0,0,1,1, col = "grey")


qmax <- 0.999
maxexp <- qexp(qmax, 0.5)
curve(dexp(x, rate = 0.5), xlim = c(0,maxexp), xlab = "x", ylab = "density", lwd = 2, 
      main = "exponential(rate = 0.5)", cex.main = 1)
x <- seq(0, maxexp, length.out = 1E3)
d <- dexp(x, rate = 0.5)
polygon(c(x, rev(x)), c(d, rep(0,length(x))), col = "grey")

for(np in 1:3){
  maxchi2 <- qchisq(qmax, df = 2 * np)
  curve(dchisq(x, df = 2 * np), xlim = c(0,maxchi2), xlab = "x", ylab = "density", lwd = 2, 
        main = paste0("χ²(df = ", 2 * np, ")"), cex.main = 1)
  x <- seq(0, maxchi2, length.out = 1E3)
  d <- dchisq(x, df = 2 * np)
  polygon(c(x, rev(x)), c(d, rep(0,length(x))), col = "grey")
}


curve(dnorm(x, mean = 1), xlim = c(-4,4), xlab = "x", ylab = "density", lwd = 2, 
      main = "normal(1,1)", cex.main = 1)
x <- seq(-4, 4, length.out = 1E3)
d <- dnorm(x, 1, 1)
polygon(c(x, rev(x)), c(d, rep(0,length(x))), col = "grey")


np <- 5000
hist(replicate(1E4, log10(hmean(runif(np)))), freq = F, breaks = 120,
     main = paste0("num. p-vals = ", np))
abline(h = 1, lty = 2, lwd = 2)
text(par("usr")[2], 1, labels = "1", pos = 4, xpd = NA)

#### extract dependence ####

#can we get z matrix from single dataset?
nrep <- 2E2
n <- 100
p <- 5E2
r <- 0.9
R <- outer(1:p, 1:p, function(i, j) r^abs(i - j))

#constant params
b <- rnorm(p) * 0
x <- rnorm(n)

#variable params
sim_dat <- function(){
  e <- simulate_X(n, p, rep(0, p), "ar1", r)
  y <- outer(x, b, "*") + e
  return(y)
}
sim_res <- function(){
  y <- sim_dat()
  out <- as.data.frame(do.call(rbind, lapply(1:p, function(i) summary(lm(y[,i] ~ x))$coefficients[2,c(1,4)])))
  return(out)
}

outs <- parallel::mclapply(1:nrep, function(nri) sim_res(), mc.cores = 8)
ests <- do.call(cbind, lapply(outs, function(res) res$Estimate))
pvals <- do.call(cbind, lapply(outs, function(res) res$`Pr(>|t|)`))
cor_ests <- cor(t(ests))
cor_pvals <- cor(t(pvals))
cor_log10pvals <- cor(t(log10(pvals)))
cor_zs <- cor(t(qnorm(pvals)))
cor_szs <- cor(t(sign(ests) * qnorm(pvals)))

plot(R[upper.tri(R)], cor_ests[upper.tri(cor_ests)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

plot(R[upper.tri(R)], cor_pvals[upper.tri(cor_pvals)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

plot(R[upper.tri(R)], cor_log10pvals[upper.tri(cor_log10pvals)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

plot(R[upper.tri(R)], cor_zs[upper.tri(cor_zs)], pch = 19, col = adjustcolor(1, 0.1))
abline(0,1, lwd = 2, lty = 2, col = 2)

cor(R[upper.tri(R)], cor_ests[upper.tri(cor_ests)])
cor(R[upper.tri(R)], cor_pvals[upper.tri(cor_pvals)])
cor(R[upper.tri(R)], cor_log10pvals[upper.tri(cor_log10pvals)])
cor(R[upper.tri(R)], cor_zs[upper.tri(cor_zs)])

# estimate covariance of null z-statistics via permutation of x
n <- 100
p <- 5E2
r <- 0.9
R <- outer(1:p, 1:p, function(i, j) r^abs(i - j))

#constant params
b <- rnorm(p) * 0
x <- rnorm(n)
y_obs <- sim_dat()

#estimate covariance matrix
null_est <- estimate_null_covariance(
  y         = y_obs,
  x         = x,
  B         = 500,
  transform = "z_from_p"
)

t_obs <- compute_tstats(y_obs, x)
df    <- nrow(y_obs) - 2
p_two <- 2 * pt(-abs(t_obs), df = df)
z_obs <- sign(t_obs) * qnorm(p_two / 2, lower.tail = FALSE)

combine_stouffer_perm(z_obs, null_est$stats_mat)

#quick check of local correlation thing
n <- 1E7
rho <- 0.95
X <- simulate_X(B = n, m = 2, mu_vec = c(0,0), 
                corr_structure = "cs", rho = rho)
cor(X)[1,2]
nbreaks <- 10
eps <- 1E-6
breaks <- seq(-5, 5, length.out = nbreaks + 1)
midp <- breaks[-1] - diff(breaks)/2
cutX1 <- as.numeric(cut(X[,1], breaks))
cutX2 <- as.numeric(cut(X[,2], breaks))
outX1 <- do.call(rbind, lapply(1:nbreaks, function(i){
  pass <- (cutX1==i)
  pass[is.na(pass)] <- F
  if(sum(pass) > 2){
    return(data.frame(cor = cor(X[pass,])[1,2],
                      n = sum(pass)))
  } else {
    return(data.frame(cor = NA, n = NA))
  }
}))
outX2 <- do.call(rbind, lapply(1:nbreaks, function(i){
  pass <- (cutX2==i)
  pass[is.na(pass)] <- F
  if(sum(pass) > 2){
    return(data.frame(cor = cor(X[pass,])[1,2],
                      n = sum(pass)))
  } else {
    return(data.frame(cor = NA, n = NA))
  }
}))
outXb <- do.call(rbind, lapply(1:nbreaks, function(i){
  pass <- (cutX1==i) & (cutX2==i)
  pass[is.na(pass)] <- F
  if(sum(pass) > 2){
    return(data.frame(cor = cor(X[pass,])[1,2],
                      n = sum(pass)))
  } else {
    return(data.frame(cor = NA, n = NA))
  }
}))
plot(midp, outX1$cor, type = "l", ylim = range(c(outX1$cor,
                                                 outX2$cor,
                                                 outXb$cor,
                                                 rho)),
     col = "red")
lines(midp, outXb$cor, col = "purple")
lines(midp, outX2$cor, col = "blue")
abline(h=rho, lty=2, col=adjustcolor(1,0.5), lwd = 2)

#### quick tests ####
# m <- 1E3
# rho <- 0.9
# mu_vec <- rep(1, m)
# sigma2 <- 1
# corr_structure <- c("cs", "ar1")[1]
# x <- simulate_X(1, 
#            m = m, 
#            mu_vec = mu_vec, 
#            rho = rho, 
#            corr_structure = corr_structure) * sqrt(sigma2)
# p <- pnorm(x)
# hist(p, breaks = 0:20/20)
# est <- estimate_gaussian_struct_cs(p)
# est$coef

#### end ####
