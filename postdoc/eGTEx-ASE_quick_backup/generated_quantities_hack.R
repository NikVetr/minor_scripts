
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### generated quantities ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# CONCENTRATION PARAMS //
#////////////////////////

#construct mean gtex effect
conc_gtex_split <- out[["conc_gtex_split"]] <- rep(NA, 2)
conc_gtex_split[1] <- out[["conc_gtex_split"]][1] <- conc_gtex / 2
conc_gtex_split[2] <- out[["conc_gtex_split"]][2] <- -conc_gtex / 2

#construct gene-specific gtex effect
conc_gene_gtex_diff_split <- out[["conc_gene_gtex_diff_split"]] <- matrix(NA, nrow = n_gene, ncol = 2)
conc_gene_gtex_diff_split[,1] <- out[["conc_gene_gtex_diff_split"]][,1] <- conc_gene_gtex_diff / 2
conc_gene_gtex_diff_split[,2] <- out[["conc_gene_gtex_diff_split"]][,2] <- -conc_gene_gtex_diff / 2
conc_gene_gtex_diff_split_vec <- out[["conc_gene_gtex_diff_split_vec"]] <- rep(NA, n)
for (i in 1:n) {
  conc_gene_gtex_diff_split_vec[i] <- out[["conc_gene_gtex_diff_split_vec"]][i] <- conc_gene_gtex_diff_split[gene_i[i], gtex[i]]
}

#construct cumulative grade effect
conc_grade_csum <- out[["conc_grade_csum"]] <- rep(NA, n_grade)
conc_grade_csum[1] <- out[["conc_grade_csum"]][1] <- 0.0
conc_grade_csum[2] <- out[["conc_grade_csum"]][2] <- 0.0
for (i in 3:n_grade) {
  conc_grade_csum[i] <- out[["conc_grade_csum"]][i] <- conc_grade_csum[i-1] + conc_grade[i-2] * sd_conc_grade + mu_conc_grade
}

#construct cumulative grade x gene effect
conc_grade_gene_csum <- out[["conc_grade_gene_csum"]] <- array(NA, dim = c(n_gene, n_grade))
for (i in 1:n_gene) {
  conc_grade_gene_csum[i,1] <- out[["conc_grade_gene_csum"]][i,1] <- 0.0
  conc_grade_gene_csum[i,2] <- out[["conc_grade_gene_csum"]][i,2] <- 0.0
}
for (k in 3:n_grade) {
  for (i in 1:n_gene) {
    conc_grade_gene_csum[i,k] <- out[["conc_grade_gene_csum"]][i,k] <- conc_grade_gene_csum[i,k-1] + conc_grade_gene[i,k-2] * sd_conc_grade_gene[k-2] * 0.5 + mu_conc_grade_gene[i] * sd_mu_conc_grade_gene * 0.5 
  }
}
conc_grade_gene_vec <- out[["conc_grade_gene_vec"]] <- rep(NA, n)
for (i in 1:n) {
  conc_grade_gene_vec[i] <- out[["conc_grade_gene_vec"]][i] <- conc_grade_gene_csum[gene_i[i], grade_k[i]]
}

#construct positive concentration vector
conc <- out[["conc"]] <- rep(NA, n)
conc <- out[["conc"]] <- 100 + exp(log_conc_base + 4 + conc_gene[gene_i] * sd_conc_gene * 0.5 + conc_indiv[indiv_j] * sd_conc_indiv * 0.5 + conc_gtex_split[gtex] + conc_grade_csum[grade_k] + conc_grade_gene_vec + conc_gene_gtex_diff_split_vec * sd_conc_gtex * 0.5)

#///////////////////
# LOCATION PARAMS //
#///////////////////

#construct mean gtex effect
loc_gtex_split <- out[["loc_gtex_split"]] <- rep(NA, 2)
loc_gtex_split[1] <- out[["loc_gtex_split"]][1] <- loc_gtex / 2
loc_gtex_split[2] <- out[["loc_gtex_split"]][2] <- -loc_gtex / 2

#construct gene-specific gtex effect
loc_gene_gtex_diff_split <- out[["loc_gene_gtex_diff_split"]] <- matrix(NA, nrow = n_gene, ncol = 2)
loc_gene_gtex_diff_split[,1] <- out[["loc_gene_gtex_diff_split"]][,1] <- loc_gene_gtex_diff / 2
loc_gene_gtex_diff_split[,2] <- out[["loc_gene_gtex_diff_split"]][,2] <- -loc_gene_gtex_diff / 2
loc_gene_gtex_diff_split_vec <- out[["loc_gene_gtex_diff_split_vec"]] <- rep(NA, n)
for (i in 1:n) {
  loc_gene_gtex_diff_split_vec[i] <- out[["loc_gene_gtex_diff_split_vec"]][i] <- loc_gene_gtex_diff_split[gene_i[i], gtex[i]]
}

#construct cumulative grade effect
loc_grade_csum <- out[["loc_grade_csum"]] <- rep(NA, n_grade)
loc_grade_csum[1] <- out[["loc_grade_csum"]][1] <- 0.0
loc_grade_csum[2] <- out[["loc_grade_csum"]][2] <- 0.0
for (i in 3:n_grade) {
  loc_grade_csum[i] <- out[["loc_grade_csum"]][i] <- loc_grade_csum[i-1] + loc_grade[i-2] * sd_loc_grade + mu_loc_grade
}

#construct cumulative grade x gene effect
loc_grade_gene_csum <- out[["loc_grade_gene_csum"]] <- array(NA, dim = c(n_gene, n_grade))
for (i in 1:n_gene) {
  loc_grade_gene_csum[i,1] <- out[["loc_grade_gene_csum"]][i,1] <- 0.0
  loc_grade_gene_csum[i,2] <- out[["loc_grade_gene_csum"]][i,2] <- 0.0
}
for (k in 3:n_grade) {
  for (i in 1:n_gene) {
    loc_grade_gene_csum[i,k] <- out[["loc_grade_gene_csum"]][i,k] <- loc_grade_gene_csum[i,k-1] + loc_grade_gene[i,k-2] * sd_loc_grade_gene[k-2] * 0.5 + mu_loc_grade_gene[i] * sd_mu_loc_grade_gene * 0.5 
  }
}
loc_grade_gene_vec <- out[["loc_grade_gene_vec"]] <- rep(NA, n)
for (i in 1:n) {
  loc_grade_gene_vec[i] <- out[["loc_grade_gene_vec"]][i] <- loc_grade_gene_csum[gene_i[i], grade_k[i]]
}

#construct positive location vector
loc <- out[["loc"]] <- rep(NA, n)
loc <- out[["loc"]] <- 0.6 + plogis(logit_loc_base / 2 + 0 + loc_gene[gene_i] * sd_loc_gene * 0.5 + loc_indiv[indiv_j] * sd_loc_indiv * 0.5  + loc_gtex_split[gtex] + loc_grade_csum[grade_k] + loc_grade_gene_vec + loc_gene_gtex_diff_split_vec * sd_loc_gtex * 0.5) / 10 * 3

#//////////////////
# MIXTURE PARAMS //
#//////////////////

#// logodds AB ////

#construct mean gtex effect
logodds_AB_gtex_split <- out[["logodds_AB_gtex_split"]] <- rep(NA, 2)
logodds_AB_gtex_split[1] <- out[["logodds_AB_gtex_split"]][1] <- logodds_AB_gtex / 2
logodds_AB_gtex_split[2] <- out[["logodds_AB_gtex_split"]][2] <- -logodds_AB_gtex / 2

#construct gene-specific gtex effect
logodds_AB_gene_gtex_diff_split <- out[["logodds_AB_gene_gtex_diff_split"]] <- matrix(NA, nrow = n_gene, ncol = 2)
logodds_AB_gene_gtex_diff_split[,1] <- out[["logodds_AB_gene_gtex_diff_split"]][,1] <- logodds_AB_gene_gtex_diff / 2
logodds_AB_gene_gtex_diff_split[,2] <- out[["logodds_AB_gene_gtex_diff_split"]][,2] <- -logodds_AB_gene_gtex_diff / 2
logodds_AB_gene_gtex_diff_split_vec <- out[["logodds_AB_gene_gtex_diff_split_vec"]] <- rep(NA, n)
for (i in 1:n) {
  logodds_AB_gene_gtex_diff_split_vec[i] <- out[["logodds_AB_gene_gtex_diff_split_vec"]][i] <- logodds_AB_gene_gtex_diff_split[gene_i[i], gtex[i]]
}

#construct cumulative grade effect
logodds_AB_grade_csum <- out[["logodds_AB_grade_csum"]] <- rep(NA, n_grade)
logodds_AB_grade_csum[1] <- out[["logodds_AB_grade_csum"]][1] <- 0.0
logodds_AB_grade_csum[2] <- out[["logodds_AB_grade_csum"]][2] <- 0.0
for (i in 3:n_grade) {
  logodds_AB_grade_csum[i] <- out[["logodds_AB_grade_csum"]][i] <- logodds_AB_grade_csum[i-1] + logodds_AB_grade[i-2] * sd_logodds_AB_grade + mu_logodds_AB_grade
}

#construct cumulative grade x gene effect
logodds_AB_grade_gene_csum <- out[["logodds_AB_grade_gene_csum"]] <- array(NA, dim = c(n_gene, n_grade))
for (i in 1:n_gene) {
  logodds_AB_grade_gene_csum[i,1] <- out[["logodds_AB_grade_gene_csum"]][i,1] <- 0.0
  logodds_AB_grade_gene_csum[i,2] <- out[["logodds_AB_grade_gene_csum"]][i,2] <- 0.0
}
for (k in 3:n_grade) {
  for (i in 1:n_gene) {
    logodds_AB_grade_gene_csum[i,k] <- out[["logodds_AB_grade_gene_csum"]][i,k] <- logodds_AB_grade_gene_csum[i,k-1] + logodds_AB_grade_gene[i,k-2] * sd_logodds_AB_grade_gene[k-2] * 0.5  + mu_logodds_AB_grade_gene[i] * sd_mu_logodds_AB_grade_gene * 0.5  
  }
}
logodds_AB_grade_gene_vec <- out[["logodds_AB_grade_gene_vec"]] <- rep(NA, n)
for (i in 1:n) {
  logodds_AB_grade_gene_vec[i] <- out[["logodds_AB_grade_gene_vec"]][i] <- logodds_AB_grade_gene_csum[gene_i[i], grade_k[i]]
}

#get mixing probability
prob_AB_vec <- out[["prob_AB_vec"]] <- rep(NA, n)
prob_AB_vec <- out[["prob_AB_vec"]] <- plogis(logodds_AB_base + 0 + logodds_AB_gene[gene_i] * sd_logodds_AB_gene * 0.5 + logodds_AB_indiv[indiv_j] * sd_logodds_AB_indiv * 0.5 + logodds_AB_gtex_split[gtex] + logodds_AB_grade_csum[grade_k] + logodds_AB_grade_gene_vec + logodds_AB_gene_gtex_diff_split_vec * sd_logodds_AB_gtex * 0.5)

#// conditional logodds LoH ////

#construct mean gtex effect
cond_logodds_loh_gtex_split <- out[["cond_logodds_loh_gtex_split"]] <- rep(NA, 2)
cond_logodds_loh_gtex_split[1] <- out[["cond_logodds_loh_gtex_split"]][1] <- cond_logodds_loh_gtex / 2
cond_logodds_loh_gtex_split[2] <- out[["cond_logodds_loh_gtex_split"]][2] <- -cond_logodds_loh_gtex / 2

#construct gene-specific gtex effect
cond_logodds_loh_gene_gtex_diff_split <- out[["cond_logodds_loh_gene_gtex_diff_split"]] <- matrix(NA, nrow = n_gene, ncol = 2)
cond_logodds_loh_gene_gtex_diff_split[,1] <- out[["cond_logodds_loh_gene_gtex_diff_split"]][,1] <- cond_logodds_loh_gene_gtex_diff / 2
cond_logodds_loh_gene_gtex_diff_split[,2] <- out[["cond_logodds_loh_gene_gtex_diff_split"]][,2] <- -cond_logodds_loh_gene_gtex_diff / 2
cond_logodds_loh_gene_gtex_diff_split_vec <- out[["cond_logodds_loh_gene_gtex_diff_split_vec"]] <- rep(NA, n)
for (i in 1:n) {
  cond_logodds_loh_gene_gtex_diff_split_vec[i] <- out[["cond_logodds_loh_gene_gtex_diff_split_vec"]][i] <- cond_logodds_loh_gene_gtex_diff_split[gene_i[i], gtex[i]]
}

#construct cumulative grade effect
cond_logodds_loh_grade_csum <- out[["cond_logodds_loh_grade_csum"]] <- rep(NA, n_grade)
cond_logodds_loh_grade_csum[1] <- out[["cond_logodds_loh_grade_csum"]][1] <- 0.0
cond_logodds_loh_grade_csum[2] <- out[["cond_logodds_loh_grade_csum"]][2] <- 0.0
for (i in 3:n_grade) {
  cond_logodds_loh_grade_csum[i] <- out[["cond_logodds_loh_grade_csum"]][i] <- cond_logodds_loh_grade_csum[i-1] + cond_logodds_loh_grade[i-2] * sd_cond_logodds_loh_grade + mu_cond_logodds_loh_grade
}

#construct cumulative grade x gene effect
cond_logodds_loh_grade_gene_csum <- out[["cond_logodds_loh_grade_gene_csum"]] <- array(NA, dim = c(n_gene, n_grade))
for (i in 1:n_gene) {
  cond_logodds_loh_grade_gene_csum[i,1] <- out[["cond_logodds_loh_grade_gene_csum"]][i,1] <- 0.0
  cond_logodds_loh_grade_gene_csum[i,2] <- out[["cond_logodds_loh_grade_gene_csum"]][i,2] <- 0.0
}
for (k in 3:n_grade) {
  for (i in 1:n_gene) {
    cond_logodds_loh_grade_gene_csum[i,k] <- out[["cond_logodds_loh_grade_gene_csum"]][i,k] <- cond_logodds_loh_grade_gene_csum[i,k-1] + cond_logodds_loh_grade_gene[i,k-2] * sd_cond_logodds_loh_grade_gene[k-2] * 0.5 + mu_cond_logodds_loh_grade_gene[i] * sd_mu_cond_logodds_loh_grade_gene * 0.5 
  }
}
cond_logodds_loh_grade_gene_vec <- out[["cond_logodds_loh_grade_gene_vec"]] <- rep(NA, n)
for (i in 1:n) {
  cond_logodds_loh_grade_gene_vec[i] <- out[["cond_logodds_loh_grade_gene_vec"]][i] <- cond_logodds_loh_grade_gene_csum[gene_i[i], grade_k[i]]
}

#get mixing probability
cond_prob_loh_vec <- out[["cond_prob_loh_vec"]] <- rep(NA, n)
cond_prob_loh_vec <- out[["cond_prob_loh_vec"]] <- plogis(cond_logodds_loh_base + 0 + cond_logodds_loh_gene[gene_i] * sd_cond_logodds_loh_gene * 0.5 + cond_logodds_loh_indiv[indiv_j] * sd_cond_logodds_loh_indiv * 0.5 + cond_logodds_loh_gtex_split[gtex] + cond_logodds_loh_grade_csum[grade_k] + cond_logodds_loh_grade_gene_vec + cond_logodds_loh_gene_gtex_diff_split_vec * sd_cond_logodds_loh_gtex * 0.5)

#///////////////
# COMBINATION //

#///////////////////////
#overall tumor effects//
#///////////////////////

conc_grade_csum_dev <- out[["conc_grade_csum_dev"]] <- rep(NA, n_grade-1) 
conc_grade_csum_dev <- out[["conc_grade_csum_dev"]] <- conc_gtex_split[2] + conc_grade_csum[2:n_grade] 
conc_grade_csum_abs_mean <- out[["conc_grade_csum_abs_mean"]] <- rep(NA, n_grade-1)
conc_grade_csum_abs_mean <- out[["conc_grade_csum_abs_mean"]] <- 100 + exp(log_conc_base + 4 + conc_grade_csum_dev)

loc_grade_csum_dev <- out[["loc_grade_csum_dev"]] <- rep(NA, n_grade-1) 
loc_grade_csum_dev <- out[["loc_grade_csum_dev"]] <- loc_gtex_split[2] + loc_grade_csum[2:n_grade] 
loc_grade_csum_abs_mean <- out[["loc_grade_csum_abs_mean"]] <- rep(NA, n_grade-1)
loc_grade_csum_abs_mean <- out[["loc_grade_csum_abs_mean"]] <- 0.6 + plogis(logit_loc_base / 2 + 0 + loc_grade_csum_dev) / 10 * 3

logodds_AB_grade_csum_dev <- out[["logodds_AB_grade_csum_dev"]] <- rep(NA, n_grade-1) 
logodds_AB_grade_csum_dev <- out[["logodds_AB_grade_csum_dev"]] <- logodds_AB_gtex_split[2] + logodds_AB_grade_csum[2:n_grade] 
logodds_AB_grade_csum_abs_mean <- out[["logodds_AB_grade_csum_abs_mean"]] <- rep(NA, n_grade-1)
logodds_AB_grade_csum_abs_mean <- out[["logodds_AB_grade_csum_abs_mean"]] <- plogis(logodds_AB_base + 0 + logodds_AB_grade_csum_dev)

cond_logodds_loh_grade_csum_dev <- out[["cond_logodds_loh_grade_csum_dev"]] <- rep(NA, n_grade-1) 
cond_logodds_loh_grade_csum_dev <- out[["cond_logodds_loh_grade_csum_dev"]] <- cond_logodds_loh_gtex_split[2] + cond_logodds_loh_grade_csum[2:n_grade] 
cond_logodds_loh_grade_csum_abs_mean <- out[["cond_logodds_loh_grade_csum_abs_mean"]] <- rep(NA, n_grade-1)
cond_logodds_loh_grade_csum_abs_mean <- out[["cond_logodds_loh_grade_csum_abs_mean"]] <- plogis(cond_logodds_loh_base + 0 + cond_logodds_loh_grade_csum_dev)

#///////////////////
#gene-wise effects//
#///////////////////

conc_grade_gene_csum_dev <- out[["conc_grade_gene_csum_dev"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)
conc_grade_gene_csum_dev_abs_mean <- out[["conc_grade_gene_csum_dev_abs_mean"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)

loc_grade_gene_csum_dev <- out[["loc_grade_gene_csum_dev"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)
loc_grade_gene_csum_dev_abs_mean <- out[["loc_grade_gene_csum_dev_abs_mean"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)

logodds_AB_grade_gene_csum_dev <- out[["logodds_AB_grade_gene_csum_dev"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)
logodds_AB_grade_gene_csum_dev_abs_mean <- out[["logodds_AB_grade_gene_csum_dev_abs_mean"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)

cond_logodds_loh_grade_gene_csum_dev <- out[["cond_logodds_loh_grade_gene_csum_dev"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)
cond_logodds_loh_grade_gene_csum_dev_abs_mean <- out[["cond_logodds_loh_grade_gene_csum_dev_abs_mean"]] <- matrix(NA, nrow = n_gene, ncol =  n_grade-1)

for (i in 1:n_gene) {
  conc_grade_gene_csum_dev[i,] <- out[["conc_grade_gene_csum_dev"]][i,] <- ((conc_grade_gene_csum[i,2:n_grade]) + conc_grade_csum_dev)
  conc_grade_gene_csum_dev_abs_mean[i,] <- out[["conc_grade_gene_csum_dev_abs_mean"]][i,] <- (100 + exp(log_conc_base + 4 + conc_grade_gene_csum_dev[i,]))
  
  loc_grade_gene_csum_dev[i,] <- out[["loc_grade_gene_csum_dev"]][i,] <- ((loc_grade_gene_csum[i,2:n_grade]) + loc_grade_csum_dev)
  loc_grade_gene_csum_dev_abs_mean[i,] <- out[["loc_grade_gene_csum_dev_abs_mean"]][i,] <- (0.6 + plogis(logit_loc_base / 2 + 0 + loc_grade_gene_csum_dev[i,]) / 10 * 3)
  
  logodds_AB_grade_gene_csum_dev[i,] <- out[["logodds_AB_grade_gene_csum_dev"]][i,] <- ((logodds_AB_grade_gene_csum[i,2:n_grade]) + logodds_AB_grade_csum_dev)
  logodds_AB_grade_gene_csum_dev_abs_mean[i,] <- out[["logodds_AB_grade_gene_csum_dev_abs_mean"]][i,] <- (plogis(logodds_AB_base + 0 + logodds_AB_grade_gene_csum_dev[i,]))
  
  cond_logodds_loh_grade_gene_csum_dev[i,] <- out[["cond_logodds_loh_grade_gene_csum_dev"]][i,] <- ((cond_logodds_loh_grade_gene_csum[i,2:n_grade]) + cond_logodds_loh_grade_csum_dev)
  cond_logodds_loh_grade_gene_csum_dev_abs_mean[i,] <- out[["cond_logodds_loh_grade_gene_csum_dev_abs_mean"]][i,] <- (plogis(cond_logodds_loh_base + 0 + cond_logodds_loh_grade_gene_csum_dev[i,]))
}