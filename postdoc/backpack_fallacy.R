#simulate data
n <- 1E7
true_mean <- rcauchy(n)
n_sample <- 8^(-1+1:5) * 2
sample_stats <- do.call(rbind, lapply(n_sample, function(i) data.frame(n = i,
                                                                       true_mean = true_mean,
                                                                       true_var = 1,
                                                                       sample_mean = rnorm(n = n, mean = true_mean, sd = 1 / sqrt(i)),
                                                                       sample_var = rchisq(n, i-1) / (i-1))))

#initialize 2-tailed comparisons
comparisons <- do.call(rbind, lapply(seq_along(n_sample), function(i) n*(i-1) + data.frame(s1 = sample(1:n), s2 = sample(1:n))))
comparisons$n1 <- sample_stats$n[comparisons$s1]
comparisons$n2 <- sample_stats$n[comparisons$s2]
comparisons$true_diff_means <- sample_stats$true_mean[comparisons$s1] - sample_stats$true_mean[comparisons$s2]

#perform t-tests
pooled_vars <- ((sample_stats$n[comparisons$s1] - 1) * sample_stats$sample_var[comparisons$s1] + 
  (sample_stats$n[comparisons$s2] - 1) * sample_stats$sample_var[comparisons$s2]) / 
  (sample_stats$n[comparisons$s1] + sample_stats$n[comparisons$s2] - 2)
SEs <- sqrt(pooled_vars * (1 / sample_stats$n[comparisons$s1] + 1 / sample_stats$n[comparisons$s2]))
sample_diff_means <- sample_stats$sample_mean[comparisons$s1] - sample_stats$sample_mean[comparisons$s2]
ts <- sample_diff_means / SEs
dfs <- sample_stats$n[comparisons$s1] + sample_stats$n[comparisons$s2] - 2
log_pvals_1tail <- pt(abs(ts), dfs, lower.tail = FALSE, log.p = TRUE)
comparisons$log10_pvals <- (log1p(exp(log_pvals_1tail)) + log_pvals_1tail) * log10(exp(1))
comparisons <- comparisons[order(comparisons$log10_pvals, decreasing = T),]
comparisons <- comparisons[comparisons$true_diff_means != 0,]
comparisons$abs_true_diff_means <- abs(comparisons$true_diff_means)

#compare t-test output
nc <- nrow(comparisons)
comparisons2 <- data.frame(n1 = comparisons$n1[1:(nc-1)], 
                           n2 = comparisons$n1[2:(nc)])
comparisons2$bigger_sample <- as.integer(comparisons2$n2 > comparisons2$n1) + 1
comparisons2$name <- ""
comparisons2$name[comparisons2$bigger_sample == 1] <- paste0(comparisons2$n2[comparisons2$bigger_sample == 1], 
                                                             " vs. ", 
                                                             comparisons2$n1[comparisons2$bigger_sample == 1])
comparisons2$name[comparisons2$bigger_sample == 2] <- paste0(comparisons2$n1[comparisons2$bigger_sample == 2], 
                                                             " vs. ", 
                                                             comparisons2$n2[comparisons2$bigger_sample == 2])

comparisons2$diff_true_effect <- comparisons$abs_true_diff_means[1:(nc-1)] - comparisons$abs_true_diff_means[2:(nc)]
comparisons2$logpval <- -(comparisons$log10_pvals[1:(nc-1)] + comparisons$log10_pvals[2:(nc)]) / 2
comparisons2$t1_bigger_effect <- sign(comparisons2$diff_true_effect) == 1
comparisons2$t1_bigger_sample <- comparisons2$n1 > comparisons2$n2
comparisons2 <- comparisons2[comparisons2$n1 != comparisons2$n2,]
comparisons2$smaller_samp_bigger_effect <- comparisons2$t1_bigger_effect != comparisons2$t1_bigger_sample

comparisons2 <- split(comparisons2, comparisons2$name)
sapply(comparisons2, nrow)

#drag sliding window
ww <- 1
windows <- 0:500/5
logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (1 + exp(x))
line_vals <- parallel::mclapply(comparisons2, function(d){
  prop_smaller_bigger <- sapply(windows, function(wi){
    subd <- d[(d$logpval > (wi - ww)) & (d$logpval < (wi + ww)),]
    logit(mean(subd$smaller_samp_bigger_effect))
  })
}, mc.cores = 8)
inf_hit <- sapply(line_vals, function(x) min(which(x == Inf)))
largest_p <- windows[max(inf_hit)]
  
#### plot results ####
par(mar = c(6,6,2,6), mfrow = c(2,1), xpd = NA)
cols <- c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(3, "Set2"))[1:10]

#first plot
yr <- range(invlogit(unlist(line_vals)[abs(unlist(line_vals)) != Inf]), na.rm = T)
plot(NULL, xlim = c(0,largest_p), 
     ylim = yr,
     axes = F, xlab = "", ylab = "")
axis(side = 1, at = pretty(c(0,largest_p)), line = 2)
axis(side = 4, at = pretty(yr), line = 1)
mtext("Pr(smaller sample has\nlarger true effect)", side = 4, line = 4.5)
mtext(latex2exp::TeX("-log$_{10}$(p-value)"), side = 1, line = 4.5)

#draw lines
ypos_labs <- invlogit(sapply(line_vals, head, 1))
repfac <- 300
ypos_labs <- FField::FFieldPtRep(cbind(0, ypos_labs) * repfac)$y / repfac
for(i in 1:length(line_vals)){
  lines(x = windows, invlogit(line_vals[[i]]), col = cols[i])
  text(y = ypos_labs[i], x = -0.5, labels = names(line_vals)[i], pos = 2, col = cols[i])
  segments(-0.6, ypos_labs[i], 0, invlogit(line_vals[[i]][1]), lty = 2, col = cols[i])
}

#second plot
yr <- range((unlist(line_vals)[abs(unlist(line_vals)) != Inf]), na.rm = T)
plot(NULL, xlim = c(0,largest_p), 
     ylim = yr,
     axes = F, xlab = "", ylab = "")
axis(side = 1, at = pretty(c(0,largest_p)), line = 1.5)
axis(side = 4, at = pretty(yr), line = 1)
mtext("log-odds(smaller sample has\nlarger true effect)", side = 4, line = 4.5)
mtext(latex2exp::TeX("-log$_{10}$(p-value)"), side = 1, line = 4.5)


#draw lines
ypos_labs <- (sapply(line_vals, head, 1))
repfac <- 4
ypos_labs <- FField::FFieldPtRep(cbind(0, ypos_labs) * repfac)$y / repfac
ypos_labs <- ypos_labs + (ypos_labs - mean(ypos_labs)) / 10
for(i in 1:length(line_vals)){
  yvals <- line_vals[[i]]
  monotonic_up_until <- min(which(sign(diff(yvals)) == -1))
  lines(x = windows[1:monotonic_up_until], line_vals[[i]][1:monotonic_up_until], col = cols[i])
  text(y = ypos_labs[i], x = -0.5, labels = names(line_vals)[i], pos = 2, col = cols[i])
  segments(-0.6, ypos_labs[i], 0, (line_vals[[i]][1]), lty = 2, col = cols[i])
  arrows(x0 = windows[monotonic_up_until-1], 
         y0 = yvals[monotonic_up_until-1], 
         x1 = windows[monotonic_up_until], 
         y1 = yvals[monotonic_up_until], length = 0.1, angle = 30, col = cols[i])
}

#label legend
pu <- par("usr")
text(x = pu[1] - diff(pu[1:2])/10, y = pu[4], labels = "sample size\nin each group\nin each t-test", pos = 4)
arrows(x0 = pu[1] - diff(pu[1:2])/30, 
       y0 = pu[4] - diff(pu[3:4])/8, 
       x1 = pu[1] - diff(pu[1:2])/15, 
       y1 = pu[4] - diff(pu[3:4])/2.5, 
       length = 0.1, angle = 30, col = 1, lwd = 2)
arrows(x0 = pu[1] - diff(pu[1:2])/30, 
       y0 = pu[4] + diff(pu[3:4])/6, 
       x1 = pu[1] - diff(pu[1:2])/15, 
       y1 = pu[4] + diff(pu[3:4])/1.75, 
       length = 0.1, angle = 30, col = 1, lwd = 2)
text(x = pu[1] + diff(pu[1:2])/10, y = pu[4] + diff(pu[3:4])/20, 
     labels = paste0("(so eg \"2 vs 16\" means that t-test #1 had n=2 in each group\n",
     "                                      and t-test #2 has n=16 in each group)"), pos = 4)


#horiz axis is -log10(p-value), vertical axis is prop small sample effect > large sample effect
#different lines for each n_sample comparison