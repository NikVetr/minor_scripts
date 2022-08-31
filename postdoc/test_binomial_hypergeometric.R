run_tests <- function(x, y, n, nbeta = 1000, bayes_prior = c(1,1), CI = 0.95, calculate_coverage = F){
  
  #create contingency table
  d <- matrix(c(length(intersect(x,y)), length(setdiff(x,y)), length(setdiff(y,x)), n - length(union(x,y))), 2, 2)
  
  #find posterior dist of difference in p
  rb1 <- rbeta(nbeta, bayes_prior[1] + d[1,1], bayes_prior[2] + d[2,1])
  rb2 <- rbeta(nbeta, bayes_prior[1] + d[1,2], bayes_prior[2] + d[2,2])
  rbd <- sort(rb1 - rb2)
  
  #fisher's exact test
  FET_out <- fisher.test(d, conf.level = CI, alternative = "two.sided")
  
  #chi2 test
  chi2_out <- suppressWarnings(chisq.test(d,correct=F))
  
  #calculate coverage
  if(calculate_coverage){
    CIs <- 1:99/100
    FET_confint <- suppressWarnings(CIs[min(which(sapply(CIs, function(confint) sum(c(fisher.test(d, conf.level = confint)$conf.int) > 1) == 1)))])
    if(is.na(FET_confint)){FET_confint <- 1}
    Bayes_confint <- abs(which.min(abs(rbd)) / nbeta * 2 - 1)
  }
  
  #return results
  c(bayes = mean(rbd > 0), 
    fisher = FET_out$p.value,
    chi2 = chi2_out$p.value,
    in_bayes_CI = ifelse(calculate_coverage, Bayes_confint, sum(quantile(rbd,c((1-CI)/2,(1+CI)/2)) > 0) == 1),
    in_fisher_CI = ifelse(calculate_coverage, FET_confint, sum(c(FET_out$conf.int) > 1) == 1))
  
}


sim_2x2_discr_unif_sample <- function(n = 1000, s = NA, ...){

  #sample how many hits are in each population
  if(all(is.na(s))){s = sample(0:n, 2, T)}
  
  #sample hits for each population
  x <- sample(1:n, size = s[1], replace = F)
  y <- sample(1:n, size = s[2], replace = F)
  
  run_tests(x = x, y = y, n = n, ...)
  
}


sim_2x2_binom_sample <- function(n = 1000, shapes1 = c(1,1), shapes2 = c(1,1), ...){
  
  #sample hits for each population
  x <- which(rbinom(n = n, size = 1, prob = rbeta(1,shapes1[1],shapes1[2])) == 1)
  y <- which(rbinom(n = n, size = 1, prob = rbeta(1,shapes2[1],shapes2[2])) == 1)
  
  run_tests(x = x, y = y, n = n, ...)
  
}

teatasting_out <- data.frame(t(replicate(1E4, sim_2x2_discr_unif_sample(n=1E3, s = c(500, 500)))))
discr_unif_sample_out <- data.frame(t(replicate(1E4, sim_2x2_discr_unif_sample(n=1E3))))
binom_sample_out <- data.frame(t(replicate(1E4, sim_2x2_binom_sample(n=1E3))))

par(mfrow = c(3,3))

mean(teatasting_out$in_bayes_CI)
mean(teatasting_out$in_fisher_CI)
hist(teatasting_out$bayes, breaks = 10, freq = F, xlab = "posterior probability positive difference", main = "Tea Tasting -- Bayes")
hist(teatasting_out$fisher, breaks = 10, freq = F, xlab = "FET p-value", main = "Tea Tasting -- FET")
hist(teatasting_out$chi2, breaks = 10, freq = F, xlab = "Chi2 p-value", main = "Tea Tasting -- Chi2")


mean(discr_unif_sample_out$in_bayes_CI)
mean(discr_unif_sample_out$in_fisher_CI)
hist(discr_unif_sample_out$bayes, breaks = 10, freq = F, xlab = "posterior probability positive difference", main = "Discrete Uniform Tea Tasting -- Bayes")
hist(discr_unif_sample_out$fisher, breaks = 10, freq = F, xlab = "FET p-value", main = "Discrete Uniform Tea Tasting -- FET")
hist(discr_unif_sample_out$chi2, breaks = 10, freq = F, xlab = "Chi2 p-value", main = "Discrete Uniform Tea Tasting -- Chi2")

mean(binom_sample_out$in_bayes_CI)
mean(binom_sample_out$in_fisher_CI)
hist(binom_sample_out$bayes, breaks = 10, freq = F, xlab = "posterior probability positive difference", main = "Binomial Sample -- Bayes")
hist(binom_sample_out$fisher, breaks = 10, freq = F, xlab = "FET p-value", main = "Binomial Sample -- FET")
hist(binom_sample_out$chi2, breaks = 10, freq = F, xlab = "Chi2 p-value", main = "Binomial Sample -- Chi2")
