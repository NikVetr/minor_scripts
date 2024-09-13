# Required libraries
library(exact2x2)
library(Barnard)
library(ggplot2)
library(gridExtra)
library(data.table)


# function to simulate data
simulate_2x2 <- function(margin_fix, total = NA, prob = NA){
  
  if(is.na(total)){
    total <- sample(20:200, 1)
  }
  
  if(is.na(prob)){
    prob <- runif(1, 0.2, 0.8)
  }
  
  
  if (margin_fix == "none") {
    indiv_counts <- data.table(rbinom(total, size = 1, prob = prob),
                               rbinom(total, size = 1, prob = prob))
  } else if (margin_fix == "one") {
    row_totals <- round(prob * total)
    row_totals <- c(row_totals, total - row_totals)
    indiv_counts <- data.table(c(rep(0, row_totals[1]), rep(1, row_totals[2])),
                            rbinom(total, size = 1, prob = prob))
  } else if (margin_fix == "both") {
    row_totals <- round(prob * total)
    row_totals <- c(row_totals, total - row_totals)
    indiv_counts <- data.table(c(rep(0, row_totals[1]), rep(1, row_totals[2])), #only need to sample one row
                               sample(c(rep(0, row_totals[1]), rep(1, row_totals[2]))))
  }
  
  aggregate_counts <- indiv_counts[, .N, by = .(V1, V2)]
  table_2x2 <- matrix(NA, nrow = 2, ncol = 2)
  table_2x2[as.matrix(aggregate_counts[,1:2])+1] <- unlist(aggregate_counts[,3])
  return(table_2x2)
}

# Function to perform simulations
mcprint <- function(...){
  system(sprintf('printf "%s"', paste0(..., collapse="")))
}
simulate_tests <- function(n_sims, alpha_levels, test_fun, margin_fix, total = NA, prob = NA) {
  
  p_values <- unlist(parallel::mclapply(1:n_sims, function(i){
    if(i %% ceiling(n_sims / 10) == 0) mcprint(paste0(i, " "))
    return(test_fun(simulate_2x2(margin_fix = margin_fix, total = total, prob = prob)))
    }, mc.cores = 8))
  
  
  
  false_positive_rates <- sapply(alpha_levels, function(alpha) {
    
    mean(p_values < alpha)
  })
  
  return(false_positive_rates)
}

# Fisher's test function
fisher_test_fun <- function(table_2x2) fisher.test(table_2x2)$p.value

# Barnard's test function
barnard_test_fun <- function(table_2x2){
  junk_output <- capture.output(pval <- Barnard::barnard.test(
    table_2x2[1,1], 
    table_2x2[1,2], 
    table_2x2[2,1], 
    table_2x2[2,2])$p.value[2])
  return(pval)
}

chi2_test_fun <- function(table_2x2){
  chisq.test(table_2x2, correct = F)$p.val
}


# Boschloo's test function
boschloo_test_fun <- function(table_2x2) boschloo(x1 = table_2x2[1,1], 
                                                  n1 = sum(table_2x2[,1]), 
                                                  x2 = table_2x2[1,2], 
                                                  n2 = sum(table_2x2[,2]))$p.value

# specify conditions
alpha_levels <- 1:10/100
n_sims <- 500
totals <- round(4*2^(1:6))
conditions <- c("none", "one", "both")
tests <- list(fisher = fisher_test_fun, 
              barnard = barnard_test_fun, 
              boschloo = boschloo_test_fun,
              chi2 = chi2_test_fun)

results <- list()
for (cond in conditions) {
  cat(paste0("\n", cond, "\n"))
  for (test_name in names(tests)) {
    cat(paste0("\n", test_name, "\n"))
    results[[test_name]][[cond]] <- simulate_tests(n_sims, 
                                                   alpha_levels, 
                                                   tests[[test_name]], 
                                                   cond,
                                                   total = 1000,
                                                   prob = 0.5)
  }
  cat("\n")
}

# Plot results
plots <- list()
for (test_name in names(tests)) {
  for (cond in conditions) {
    data <- data.frame(alpha = alpha_levels, false_positive_rate = results[[test_name]][[cond]])
    p <- ggplot(data, aes(x = alpha, y = false_positive_rate)) +
      geom_line() +
      ggtitle(paste(test_name, "test\n# margins fixed: ", cond)) +
      xlab("Alpha Significance Threshold") +
      ylab("False + Rate") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      theme_minimal()
    plots[[paste(test_name, cond, sep = "_")]] <- p
  }
}

# Arrange plots in a 3x3 grid
grid.arrange(plots$fisher_none, plots$fisher_one, plots$fisher_both,
             plots$barnard_none, plots$barnard_one, plots$barnard_both,
             plots$boschloo_none, plots$boschloo_one, plots$boschloo_both,
             plots$chi2_none, plots$chi2_one, plots$chi2_both, nrow = length(tests))
