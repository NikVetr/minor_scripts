#specifications
setwd("~/repos/egtex-ase")
sample_prop_tail_to_remove <- 0

#### Load Packages ####
library(abind)
library(cmdstanr)
library(future)
library(posterior)
library(tidyverse)
library(fs)
library(data.table)
library(broom)
library(VGAM)
library(matrixStats)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(S4Vectors)
library(GenomicFeatures)
library(GenomicRanges)

#### Set up future multisession ####
if(!inherits(plan(), "multisession")){
  plan(multisession)  
}

#### specify functions ####
source("scripts/functions.R")
source("~/repos/Stan2R/R/functions.R")
setNm <- function(x) setNames(x,x)
invlogit <- function(x){exp(x) / (1 + exp(x))}

#### create useful variables ####
stan_model_dir <- "Stan/models/"
stan_data_dir <- "Stan/data/"
stan_output_dir <- "Stan/output/"
stan_progress_dir <- "Stan/progress/"
tissue_codes_df <- fread("tissue_codes.txt")
tissue_codes <- setNames(tissue_codes_df$code, tissue_codes_df$name)

#### inspect GTEx-gene model ####

#first, let's extract the params for meta-analysis
stan_output_dir <- "Stan/output/genes/"
stan_progress_dir <- "Stan/progress/genes/"

#specify model name params
use_rnaseq <- F
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}

base_fit <- T
model_index <- 4
model_names <- list(
  "ase_param-expansion_1gene", #1
  "ase_param-expansion_by-gene_eQTL", #2
  "ase_param-expansion_by-gene_raw-params_eQTL", #3
  "ase_param-expansion_by-gene_raw-params_resid-indivs_eQTL" #4
)

model_name <- model_names[[model_index]]
sample_prop_tail_to_remove <- 0
model_path <- paste0(stan_model_dir, model_name, ".stan")
model_path_no_eQTLHets <- gsub("eQTL", "no-eQTL", model_path)
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
cat(model_path)

#set MCMC inputs
tissue_names <- setNames(names(tissue_codes), tissue_codes)
tissue_codes_to_use <- tissue_codes
gene_moments <- rhats <- ess_bulks <- list()
gene_params_base <- c("mu" , "sigma" , "logodds_AB" , "cond_logodds_loh")
gene_params_ext <- c("base", "eQTL_diff")
gene_params <- apply(expand.grid(gene_params_base, gene_params_ext), 1, paste0, collapse = "_")
all_output_files <- list.files(stan_output_dir)

for(tiss_j in tissue_codes_to_use){
  
  print(tiss_j)
  
  base_file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                           ifelse(use_rnaseq, "_rnaseq", ""))
  file_name <- paste0(ifelse(base_fit, "1st-fit_", ""), base_file_name)  
  
  found_files <- all_output_files[startsWith(all_output_files, file_name)]
  escaped_file_name <- gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", file_name))
  gsub_regex <- paste0(escaped_file_name, "_|\\.cmdStanR\\.summ|\\.cmdStanR\\.fit|\\.stanfit|.summ")
  genes_fit <- unique(gsub(pattern = gsub_regex, replacement = "", x = found_files))
  
  if(length(genes_fit) > 0){
    
    for(gi in 1:length(genes_fit)){
      
      gene_symbol <- genes_fit[gi]
      gene_file_name <- paste0(file_name, "_", gene_symbol)
      
      #load in input data
      dat_path <- paste0(stan_data_dir, "genes/", gene_file_name, ".json")
      if(file.exists(dat_path)){
        dat_tiss <- jsonlite::fromJSON(dat_path)  
      }
      
      #load and process mcmc output
      load(file = paste0(stan_output_dir, gene_file_name,".cmdStanR.fit"))
      samps <- data.frame(as_draws_df(fit$draws()))
      gene_samps <- samps[,intersect(colnames(samps), gene_params)]
      gene_moments <- c(gene_moments, list(data.frame(apply(gene_samps, 2, function(x) c(mean = mean(x),
                                                                                         var = var(x), 
                                                                                         skew = e1071::skewness(x),
                                                                                         kurtosis = e1071::kurtosis(x),
                                                                                         gr0 = mean(x>0)), 
                                                            simplify = F), 
                                                      gene = gene_symbol, tissue = tiss_j)))
      
      load(file = paste0(stan_output_dir, gene_file_name, ".cmdStanR.summ"))
      
      ess_bulks <- c(ess_bulks, setNames(min(summ$ess_bulk), 
                                         summ$variable[which.min(summ$ess_bulk)]))
      rhats <- c(rhats, setNames(max(summ$rhat), 
                                 summ$variable[which.max(summ$rhat)]))
    }
  }
  
}

sort(unlist(ess_bulks))
sort(unlist(rhats), T)
split_gene_moments <- split(gene_moments, unlist(lapply(gene_moments, ncol)))
split_gene_moments <- lapply(split_gene_moments, do.call, what = rbind)
ncols_gene_moments <- sapply(split_gene_moments, ncol)
missing_cols <- setdiff(colnames(split_gene_moments[[which.max(ncols_gene_moments)]]),
                        colnames(split_gene_moments[[which.min(ncols_gene_moments)]]))
split_gene_moments[[which.min(ncols_gene_moments)]] <- cbind(
  split_gene_moments[[which.min(ncols_gene_moments)]], lapply(setNm(missing_cols), function(x) 
    rep(NA, nrow(split_gene_moments[[which.min(ncols_gene_moments)]]))
  )
)

#recompile data and split by type
gene_moments <- rbind(split_gene_moments[[1]], split_gene_moments[[2]][,colnames(split_gene_moments[[1]])])
gene_moments <- split(gene_moments, gsub(pattern = "[0-9]", "", rownames(gene_moments)))

#some quick visual checks
hist(gene_moments$mean$mu_base, breaks = 100)

hist(exp(gene_moments$mean$sigma_base), breaks = 20)
hist(gene_moments$mean$sigma_base / sqrt(gene_moments$var$sigma_base), breaks = 20)

hist(invlogit(gene_moments$mean$logodds_AB_base), breaks = 100)
hist(gene_moments$mean$logodds_AB_base / sqrt(gene_moments$var$logodds_AB_base), breaks = 100)

hist(invlogit(gene_moments$mean$cond_logodds_loh_base), breaks = 100)
hist(gene_moments$mean$cond_logodds_loh_base / sqrt(gene_moments$var$cond_logodds_loh_base), breaks = 100)

#look at eQTL effects (state 2 is when they are het for eQTL,
#and the diff gets turned into [-*_eQTL_diff/2, *_eQTL_diff/2])
#so it represents how much more (if +) or less (if -) the eQTL hets are than the eQTL homs
hist(gene_moments$gr$mu_eQTL_diff, breaks = 100)
hist(gene_moments$gr$sigma_eQTL_diff, breaks = 100)
hist(gene_moments$gr$logodds_AB_eQTL_diff, breaks = 100)
hist(gene_moments$gr$cond_logodds_loh_eQTL_diff, breaks = 100)

hist(gene_moments$mean$mu_eQTL_diff, breaks = 100)
hist(gene_moments$mean$sigma_eQTL_diff, breaks = 100)
hist(gene_moments$mean$logodds_AB_eQTL_diff, breaks = 100)
gene_moments$mean[gene_moments$mean$logodds_AB_eQTL_diff > 0.7 & !is.na(gene_moments$mean$logodds_AB_eQTL_diff),]
hist(gene_moments$mean$cond_logodds_loh_eQTL_diff, breaks = 100)


# now let's prep files for meta-analysis #
gene_dists <- list()
gene_params_base <- c("mu" , "sigma" , "logodds_AB" , "cond_logodds_loh")
gene_params_ext <- c("base", "eQTL_diff")
gene_params <- apply(expand.grid(gene_params_base, gene_params_ext), 1, paste0, collapse = "_")
all_output_files <- list.files(stan_output_dir)

for(tiss_j in tissue_codes_to_use){
  
  print(tiss_j)
  
  base_file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                           ifelse(use_rnaseq, "_rnaseq", ""))
  file_name <- paste0(ifelse(base_fit, "1st-fit_", ""), base_file_name)  
  
  found_files <- all_output_files[startsWith(all_output_files, file_name)]
  escaped_file_name <- gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", file_name))
  gsub_regex <- paste0(escaped_file_name, "_|\\.cmdStanR\\.summ|\\.cmdStanR\\.fit|\\.stanfit|.summ")
  genes_fit <- unique(gsub(pattern = gsub_regex, replacement = "", x = found_files))
  
  if(length(genes_fit) > 0){
    
    for(gi in 1:length(genes_fit)){
      
      gene_symbol <- genes_fit[gi]
      gene_file_name <- paste0(file_name, "_", gene_symbol)
      
      #load in input data
      dat_path <- paste0(stan_data_dir, "genes/", gene_file_name, ".json")
      if(file.exists(dat_path)){
        dat_tiss <- jsonlite::fromJSON(dat_path)  
      }
      
      #load and process mcmc output
      load(file = paste0(stan_output_dir, gene_file_name,".cmdStanR.fit"))
      samps <- data.frame(as_draws_df(fit$draws()))
      params_found <- intersect(colnames(samps), gene_params)
      gene_samps <- samps[,params_found]
      
      #if there is no eQTL param listed, sample one
      found_hets <- any(grepl("eQTL", params_found))
      if(!found_hets){
        params_needed <- setdiff(gene_params, params_found)
        normal_RVs <- matrix(rnorm(nrow(gene_samps) * length(params_needed)),
                             nrow(gene_samps), length(params_needed))
        prior_SDs <- c(mu_eQTL_diff = 0.1,
                       sigma_eQTL_diff = 1,
                       logodds_AB_eQTL_diff = 0.5,
                       cond_logodds_loh_eQTL_diff = 0.5)
        prior_gene_samps <- normal_RVs * rep(prior_SDs, each = nrow(gene_samps))
        colnames(prior_gene_samps) <- names(prior_SDs)
        gene_samps <- cbind(gene_samps, prior_gene_samps)
      }
      
      #reconstruct the means in hets vs non-hets
      if(found_hets){
        new_gene_params <- cbind(gene_samps[,paste0(gene_params_base, "_base")] - 
                                   gene_samps[,paste0(gene_params_base, "_eQTL_diff")] / 2,
                                 gene_samps[,paste0(gene_params_base, "_base")] + 
                                   gene_samps[,paste0(gene_params_base, "_eQTL_diff")] / 2
        )
      } else {
        new_gene_params <- cbind(gene_samps[,paste0(gene_params_base, "_base")],
                                 gene_samps[,paste0(gene_params_base, "_base")] + 
                                   gene_samps[,paste0(gene_params_base, "_eQTL_diff")]) 
      }
      colnames(new_gene_params) <- paste0(colnames(new_gene_params), rep(c(".homo", ".het"), each = 4))
      
      gene_means <- apply(new_gene_params, 2, mean)
      gene_cov <- cov(new_gene_params)
      gene_L <- t(chol(gene_cov))
      
      gene_dists[[tiss_j]][[gene_symbol]] <- list(gene = gene_symbol, 
                                                    gene_means = gene_means, 
                                                    gene_L = gene_L)
    }
  }
}

#munge to tissue-specific lists
genewise_tissue_dats <- lapply(gene_dists, function(x){
  list(x_est = lapply(x, function(y) y$gene_means),
       x_err_L = lapply(x, function(y) y$gene_L),
       k = length(x[[1]]$gene_means),
       n = length(x))
})

path_errprop_multivariate <- "Stan/models/errprop_multivariate.stan"
mod_errprop_multivariate <- cmdstan_model(path_errprop_multivariate)

#fit first stage of error propagation model
nchains <- 4
niter <- 5E3
adapt_delta <- 0.95
max_treedepth <- 12
thin <- 1
init_mcmc <- 1
fits_errprop_multivariate <- lapply(1:length(genewise_tissue_dats), function(i){
  cat(paste0("(", i, "/", length(genewise_tissue_dats),") "))
  
  sink(tempfile())
  fit <- mod_errprop_multivariate$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                                         data = genewise_tissue_dats[[i]], 
                                         parallel_chains = nchains, 
                                         adapt_delta = adapt_delta,
                                         refresh = 100, max_treedepth = max_treedepth,
                                         thin = thin, init = init_mcmc)
  sink()
  
  return(fit)
  
})

#examine mcmc diagnostics
summs_errprop_multivariate <- lapply(1:n_groups_1, function(i){
  cat(paste0("(", i, "/", n_groups_1,") "))
  summ_sub <- fits_errprop_multivariate[[i]]$summary()
  c(ess_bulk = min(summ_sub$ess_bulk, na.rm = T), rhat = max(summ_sub$rhat, na.rm = T))
})
summs_errprop_multivariate

#extract samples and process
samps_errprop_multivariate <- lapply(1:n_groups_1, function(i){
  cat(paste0("(", i, "/", n_groups_1,") "))
  samps <- data.frame(as_draws_df(fits_errprop_multivariate[[i]]$draws()))
  samps[,grepl("_mu", colnames(samps))]
})

#fit second stage of error propagation model
multiv_input_recurs <- lapply(1:n_groups_1, function(i){
  
  cat(paste0("(", i, "/", n_groups_1,") "))
  
  list(
    x_est = apply(samps_errprop_multivariate[[i]], 2, mean),
    x_err_L = t(chol(cov(samps_errprop_multivariate[[i]])))
  )
  
})

dat_errprop_multivariate_recurs <- list(
  x_est = lapply(multiv_input_recurs, function(foo) foo$x_est),
  x_err_L = lapply(multiv_input_recurs, function(foo) foo$x_err_L),
  n = length(multiv_input_recurs),
  k = length(multiv_input_recurs[[1]]$x_est)
)

fit_errprop_multivariate_recurs <- mod_errprop_multivariate$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                                                                   data = dat_errprop_multivariate_recurs, parallel_chains = nchains, 
                                                                   adapt_delta = adapt_delta,
                                                                   refresh = 100, max_treedepth = max_treedepth,
                                                                   thin = thin, init = init_mcmc)

summ_errprop_multivariate_recurs <- fit_errprop_multivariate_recurs$summary()
summ_errprop_multivariate_recurs[order(summ_errprop_multivariate_recurs$ess_bulk),]
summ_errprop_multivariate_recurs[order(summ_errprop_multivariate_recurs$rhat, decreasing = T),]
samps_errprop_multivariate_recurs <- data.frame(as_draws_df(fit_errprop_multivariate_recurs$draws()))

#TODO re-fit first stage model using empirical hyperpriors from second stage

#so take posterior mean of the mean of the hyper hyper distribution 
#and posterior mean of the covariance of the hyper hyper distribution 
#and set it as the prior for the first stage 

path_errprop_multivariate_empirical <- "Stan/models/errprop_multivariate_empirical.stan"
mod_errprop_multivariate_empirical <- cmdstan_model(path_errprop_multivariate_empirical)

#then take *that* estimate for the tissue mean
#and refit all the gene-wise estimates again with it
#can compare all the tissue means to each other in plots
#and also compare all the gene-wise estimates in other plots


#### inspect GTEx model ####
init_from_prior <- F
tiss_j <- tissue_codes[4]
use_EB <- F
use_flat_MFVB <- F
model_index <- 17
model_names <- list(
  "ase-gtex_inf-priors_no-indiv", #1
  "ase-gtex_inf-priors", #2
  "ase-gtex", #3
  "ase-gtex_centered", #4
  "ase-gtex_offset-conc", #5
  "ase-gtex_nearly-balanced", #6
  "ase-gtex_nearly-balanced_LoH", #7
  "ase-gtex_nearly-balanced_LoH_nested", #8
  "ase-gtex_nearly-balanced_LoH_noIndiv", #9
  "ase-gtex_nearly-balanced_LoH_nested_noIndiv", #10
  "ase-gtex_nearly-balanced_LoH_indivXgene", #11
  "ase-gtex_nearly-balanced_LoH_nested_indivXgene", #12
  "ase-gtex_nearly-balanced_LoH_nested_centered", #13
  "ase-gtex_nearly-balanced_LoH_nested_conc-scale-bounds", #14
  "ase-gtex_nearly-balanced_LoH_nested_sample-probs", #15
  "ase-gtex_nearly-balanced_LoH_nested_sample-probs-truncnorm", #16
  "ase-gtex_nearly-balanced_LoH_nested", #17
  "ase-gtex_nearly-balanced_LoH_nested_tighter-priors" #18
)

model_name <- model_names[[model_index]]
sample_prop_tail_to_remove <- 0

model_path <- paste0(stan_model_dir, model_name, ".stan")
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")

model_blocks <- parse_stan_blocks(clean_stan(model_lines))
model_blocks <- names(model_blocks)[sapply(model_blocks, length) != 0]
setNm <- function(x) setNames(x, x)
model_objects <- lapply(setNm(model_blocks), retrieve_block_objects, stan_code = model_lines)

#### inspect all GTEx diagnostics ####

#chug through all the model fits and figure out which fit well and which did not
use_rnaseq <- F
summlist <- lapply(setNm(tissue_codes), function(x) NA)
names(summlist) <- tissue_codes
for(tiss_j in tissue_codes){
  
  cat(paste0(tiss_j, " "))
  file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                      ifelse(init_from_prior, "_prior-init", "_reg-init"),
                      ifelse(use_EB,
                             ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                             "_FB"),
                      ifelse(use_rnaseq, "_rnaseq", ""))
  summ_file_path <- paste0(stan_output_dir, file_name, ".cmdStanR.summ")
  if(file.exists(summ_file_path)){
    load(summ_file_path)
    badsumm <- summ[summ$ess_bulk < 100 | summ$rhat > 1.01,]
    summlist[[tiss_j]] <- badsumm[order(badsumm$ess_bulk),]
  }
  
}

summlist <- summlist[sapply(summlist, function(x) "data.frame" %in% class(x))]
badsummlist <- summlist[sapply(summlist, nrow) > 1]
badsummlist 

#### extract GTEx genesets ####
run_extraction <- F

if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
# tissue_codes <- setdiff(tissue_codes, names(badsummlist))


if(run_extraction){
  
  #initialize container variables
  genefits <- lapply(setNm(tissue_codes), function(x) NA)
  intercepts <- lapply(setNm(tissue_codes), function(x) NA)
  foundfit <- sapply(setNm(tissue_codes), function(x) F)
  names(summlist) <- tissue_codes
  
  #params of interest
  params_of_interest <- list(
    genewise = c("conc_gene", "sd_conc_gene", "loc_gene", "sd_loc_gene", "logodds_AB_gene",
                 "sd_logodds_AB_gene", "cond_logodds_loh_gene", "sd_cond_logodds_loh_gene"),
    intercepts = c("log_conc_base", "logit_loc_base", "logodds_AB_base", "cond_logodds_loh_base")
  )
  
  #chugg through mcmc output
  for(tiss_j in tissue_codes){
    
    cat(paste0(tiss_j, " "))
    file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                        ifelse(init_from_prior, "_prior-init", "_reg-init"),
                        ifelse(use_EB,
                               ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                               "_FB"),
                        ifelse(use_rnaseq, "_rnaseq", ""))
    fit_file_path <- paste0(stan_output_dir, file_name, ".cmdStanR.fit")
    
    #process, if file exists
    if(file.exists(fit_file_path)){
      
      #load in data and dictionary
      load(fit_file_path)
      foundfit[tiss_j] <- T
      genes_tiss <- fread(file = paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_gene-map.csv"))
      
      #read in and subset mcmc draws
      samps <- as.data.table(as_draws_df(out$draws(variables = unlist(params_of_interest))))
      samps_by_param <- lapply(unlist(params_of_interest), function(var_name) 
        subset_samps(var_name = var_name, samps = samps) 
      )
      names(samps_by_param) <- unlist(params_of_interest)
      
      #now transform these (manually :/, but much more efficiently)
      transformed_params <- list(genewise = list(), intercepts = list())
      nsamp <- nrow(samps_by_param$conc_gene)
      ngene <- ncol(samps_by_param$conc_gene)
      
      # real shapes_ab = (20 + exp(6 + log_conc_allelic_balance)) * 0.5;
      # real loc_loh = 0.95 + inv_logit(logit_loc_loh) / 20;
      # real shape_loh_alpha = (20 + exp(6 + log_conc_loh)) * loc_loh;
      # real shape_loh_beta = (20 + exp(6 + log_conc_loh)) * (1-loc_loh);
      
      #gene-wise
      transformed_params[["genewise"]][["conc_gene"]] <- samps_by_param$conc_gene * matrix(unlist(samps_by_param$sd_conc_gene), 
                                                                             nrow = nsamp, 
                                                                             ncol = ngene, byrow = F) * 0.5
      transformed_params[["genewise"]][["loc_gene"]] <- samps_by_param$loc_gene * matrix(unlist(samps_by_param$sd_loc_gene), 
                                                                             nrow = nrow(samps_by_param$loc_gene), 
                                                                             ncol = ncol(samps_by_param$loc_gene), byrow = F) * 0.5
      transformed_params[["genewise"]][["logodds_AB_gene"]] <- samps_by_param$logodds_AB_gene * matrix(unlist(samps_by_param$sd_logodds_AB_gene), 
                                                                             nrow = nrow(samps_by_param$logodds_AB_gene), 
                                                                             ncol = ncol(samps_by_param$logodds_AB_gene), byrow = F) * 0.5
      transformed_params[["genewise"]][["cond_logodds_loh_gene"]] <- samps_by_param$cond_logodds_loh_gene * matrix(unlist(samps_by_param$sd_cond_logodds_loh_gene), 
                                                                             nrow = nrow(samps_by_param$cond_logodds_loh_gene), 
                                                                             ncol = ncol(samps_by_param$cond_logodds_loh_gene), byrow = F) * 0.5
      
      #intercepts
      transformed_params[["intercepts"]][["conc_base"]] <- 10 + exp(unlist(samps_by_param$log_conc_base) + 1)
      transformed_params[["intercepts"]][["loc_base"]] <- 0.60 + plogis(unlist(samps_by_param$logit_loc_base) / 2) / 10 * 3
      transformed_params[["intercepts"]][["prob_AB_base"]] <- plogis(unlist(samps_by_param$logodds_AB_base))
      transformed_params[["intercepts"]][["cond_prob_loh_base"]] <- plogis(unlist(samps_by_param$cond_logodds_loh_base))
      transformed_params[["intercepts"]][["prob_loh_base"]] <- ((1-transformed_params[["intercepts"]][["prob_AB_base"]]) * 
                                                  transformed_params[["intercepts"]][["cond_prob_loh_base"]])
      
      #gene-wise absolute
      transformed_params[["genewise"]][["abs_conc_gene"]] <- 10 + exp(transformed_params[["genewise"]][["conc_gene"]] + matrix(unlist(samps_by_param$log_conc_base), 
                                                                                         nrow = nsamp, 
                                                                                         ncol = ngene, byrow = F) + 1)
      transformed_params[["genewise"]][["abs_loc_gene"]] <- data.frame(0.6 + plogis(as.matrix(transformed_params[["genewise"]][["loc_gene"]] + matrix(unlist(samps_by_param$logit_loc_base), 
                                                                                                     nrow = nsamp, 
                                                                                                     ncol = ngene, byrow = F) / 2)) / 10 * 3)
      transformed_params[["genewise"]][["abs_prob_AB_gene"]] <- data.frame(plogis(as.matrix(transformed_params[["genewise"]][["logodds_AB_gene"]] + matrix(unlist(samps_by_param$logodds_AB_base), 
                                                                                                      nrow = nsamp, 
                                                                                                      ncol = ngene, byrow = F))))
      transformed_params[["genewise"]][["abs_cond_prob_loh_gene"]] <- data.frame(plogis(as.matrix(transformed_params[["genewise"]][["cond_logodds_loh_gene"]] + matrix(unlist(samps_by_param$cond_logodds_loh_base), 
                                                                        nrow = nsamp, 
                                                                        ncol = ngene, byrow = F))))
      transformed_params[["genewise"]][["abs_prob_loh_gene"]] <- data.frame(((1-transformed_params[["genewise"]][["abs_prob_AB_gene"]]) * 
                                                      transformed_params[["genewise"]][["abs_cond_prob_loh_gene"]]))
      
      transformed_params$genewise <- lapply(transformed_params$genewise, function(genewise_params){
        colnames(genewise_params) <- genes_tiss$gene_symbol
        return(genewise_params)
      })
      
      genefits[[tiss_j]] <- transformed_params$genewise 
      intercepts[[tiss_j]] <- transformed_params$intercepts 
      
    }
  }
  
  genefits <- genefits[names(foundfit[foundfit])]
  intercepts <- intercepts[names(foundfit[foundfit])]
  focal_samps <- list(genefits = genefits, intercepts = intercepts)
  save(focal_samps, file = paste0(stan_output_dir, "focal-samps_", model_name, "_", sample_prop_tail_to_remove, 
         ifelse(init_from_prior, "_prior-init", "_reg-init"),
         ifelse(use_EB,
                ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                "_FB"),
         ifelse(use_rnaseq, "_rnaseq", ""), ".RData"))

}


#### extract GTEx genesets, gen-quan-block style ####

#append generated quantities block
generated_quantities <- "

generated quantities{
  
  //tissue intercepts
  real conc_base = 10 + exp(log_conc_base + 1);
  real loc_base = 0.60 + inv_logit(logit_loc_base / 2) / 10 * 3;
  real prob_AB_base = inv_logit(logodds_AB_base);
  real cond_prob_loh_base = inv_logit(cond_logodds_loh_base);
  real prob_loh_base = (1 - prob_AB_base) * cond_prob_loh_base;

  //gene-wise values
  vector[n_gene] conc_gene_abs_dev = conc_gene * sd_conc_gene * 0.5;
  vector[n_gene] conc_gene_abs = 10 + exp(log_conc_base + 1 + conc_gene_abs_dev);

  vector[n_gene] loc_gene_abs_dev = loc_gene * sd_loc_gene * 0.5;
  vector[n_gene] loc_gene_abs = 0.60 + inv_logit(logit_loc_base / 2 + 0 + loc_gene_abs_dev) / 10 * 3;

  vector[n_gene] logodds_AB_gene_abs_dev = logodds_AB_gene * sd_logodds_AB_gene * 0.5;
  vector[n_gene] prob_AB_gene_abs = inv_logit(logodds_AB_base + 0 + logodds_AB_gene_abs_dev);

  vector[n_gene] cond_logodds_loh_gene_abs_dev = cond_logodds_loh_gene * sd_cond_logodds_loh_gene * 0.5;
  vector[n_gene] cond_prob_loh_gene_abs = inv_logit(cond_logodds_loh_base + 0 + cond_logodds_loh_gene_abs_dev);
  
  vector[n_gene] prob_loh_gene_abs;
  vector[n_gene] prob_slab_gene_abs;
  for(i in 1:n_gene){
    prob_loh_gene_abs[i] = (1-prob_AB_gene_abs[i]) * cond_prob_loh_gene_abs[i];
    prob_slab_gene_abs[i] = (1-prob_AB_gene_abs[i]) * (1-cond_prob_loh_gene_abs[i]);
  }

}
"

# # Specify the new model with the generated quantities block
model_string_generated_quantities <- paste0(model_string, generated_quantities, collapse = "\n")
gq_model_path <- write_stan_file(model_string_generated_quantities)
# gen_quant_model <- cmdstan_model(gen_quant_model_path)
# 
# # Load the fitted model
# file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
#                     ifelse(init_from_prior, "_prior-init", "_reg-init"),
#                     ifelse(use_EB,
#                            ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
#                            "_FB"),
#                     ifelse(use_rnaseq, "_rnaseq", ""))
# fit_file_path <- paste0(stan_output_dir, file_name, ".cmdStanR.fit")
# load(fit_file_path)
# fit <- out
# rm(out)
# 
# # load the data
# dat <- jsonlite::fromJSON(paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, ".json"))
# 
# # Generate the quantities
# gen_quant_fit <- gen_quant_model$generate_quantities(
#   fitted_params = fit$output_files(),
#   data = dat, parallel_chains = fit$num_chains()
# )
# 
# # Access the generated quantities
# generated_quantities <- as.data.table(as_draws_df(gen_quant_fit$draws()))
# 
# #check these generate the same transformed params (they should!)
# # plot(unlist(genefits[[tiss_j]]$conc_gene[,1]), c(generated_quantities$`conc_gene_abs_dev[1]`))
# # plot(unlist(genefits[[tiss_j]]$loc_gene[,1]), c(generated_quantities$`loc_gene_abs_dev[1]`))
# # plot(unlist(genefits[[tiss_j]]$logodds_AB_gene[,1]), c(generated_quantities$`logodds_AB_gene_abs_dev[1]`))
# # plot(unlist(genefits[[tiss_j]]$cond_logodds_loh_gene[,1]), c(generated_quantities$`cond_logodds_loh_gene_abs_dev[1]`))
# # plot(unlist(genefits[[tiss_j]]$abs_conc_gene[,1]), c(generated_quantities$`conc_gene_abs[1]`))
# # plot(unlist(genefits[[tiss_j]]$abs_loc_gene[,1]), c(generated_quantities$`loc_gene_abs[1]`))
# # plot(unlist(genefits[[tiss_j]]$abs_prob_AB_gene[,1]), c(generated_quantities$`prob_AB_gene_abs[1]`))
# # plot(unlist(genefits[[tiss_j]]$abs_cond_prob_loh_gene[,1]), c(generated_quantities$`cond_prob_loh_gene_abs[1]`))
# # plot(unlist(genefits[[tiss_j]]$abs_prob_loh_gene[,1]), c(generated_quantities$`prob_loh_gene_abs[1]`))
# #they all do!

#chug through these like we did before, but using the genquan code

if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
# tissue_codes <- setdiff(tissue_codes, names(badsummlist))

#initialize container variables
genefits <- lapply(setNm(tissue_codes), function(x) NA)
intercepts <- lapply(setNm(tissue_codes), function(x) NA)
foundfit <- sapply(setNm(tissue_codes), function(x) F)
names(summlist) <- tissue_codes

# Specify the new model with the generated quantities block
gq_model_path <- write_stan_file(model_string_generated_quantities)
gq_model <- cmdstan_model(gq_model_path)
gq_var_names <- retrieve_block_objects(stan_code = readLines(gq_model_path), block = "generated quantities")$var_name
gq_gene_params <- gq_var_names[grepl("gene", gq_var_names)]
gq_tiss_params <- gq_var_names[!grepl("gene", gq_var_names)]
focal_samps_path <- paste0(stan_output_dir, "focal-samps_", model_name, "_", sample_prop_tail_to_remove, 
       ifelse(init_from_prior, "_prior-init", "_reg-init"),
       ifelse(use_EB,
              ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
              "_FB"),
       ifelse(use_rnaseq, "_rnaseq", ""), ".RData")
#chugg through mcmc output
if(file.exists(focal_samps_path)){
  
  load(focal_samps_path)
  
} else {
  
  for(tiss_j in tissue_codes){
    
    cat(paste0(tiss_j, " "))
    file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                        ifelse(init_from_prior, "_prior-init", "_reg-init"),
                        ifelse(use_EB,
                               ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                               "_FB"),
                        ifelse(use_rnaseq, "_rnaseq", ""))
    fit_file_path <- paste0(stan_output_dir, file_name, ".cmdStanR.fit")
    
    #process, if file exists
    if(file.exists(fit_file_path)){
      
      #load in data and dictionary
      load(fit_file_path)
      fit <- out
      rm(out)
      foundfit[tiss_j] <- T
      genes_tiss <- fread(file = paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_gene-map.csv"))
      
      # load the data
      dat <- jsonlite::fromJSON(paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_all.json"))
      
      # Generate the quantities
      sink(tempfile()) 
      gq_fit <- gq_model$generate_quantities(
        fitted_params = fit$output_files(),
        data = dat, parallel_chains = fit$num_chains()
      )
      sink()
      
      # Access the generated quantities
      gq <- as.data.table(as_draws_df(gq_fit$draws()))
      
      #read in and subset mcmc draws
      gq_genes_by_param <- lapply(setNm(gq_gene_params), function(var_name){ 
        gene_vals <- subset_samps(var_name = var_name, samps = gq)
        colnames(gene_vals) <- setNames(genes_tiss$gene_symbol, genes_tiss$index)[sub(".*\\.", "", colnames(gene_vals))]
        gene_vals
      })
      
      gq_tiss_by_param <- lapply(setNm(gq_tiss_params), function(var_name){ 
        subset_samps(var_name = var_name, samps = gq)
      })
      
      #now record these
      transformed_params <- list(genewise = list(), intercepts = list())
      transformed_params[["genewise"]] <- gq_genes_by_param
      transformed_params[["intercepts"]] <- gq_tiss_by_param
      
      #and put them in the larger output variable
      genefits[[tiss_j]] <- transformed_params$genewise 
      intercepts[[tiss_j]] <- transformed_params$intercepts 
      
    }
  }
  
  genefits <- genefits[names(foundfit[foundfit])]
  intercepts <- intercepts[names(foundfit[foundfit])]
  focal_samps <- list(genefits = genefits, intercepts = intercepts)
  save(focal_samps, file = focal_samps_path)
  
}

#### plotting GTEx comparisons ####

#get more convenient variables
intercepts <- focal_samps$intercepts
genefits <- focal_samps$genefits
rm(focal_samps); gc()

param_names <- setNm(names(intercepts$`adipose-visceral-(omentum)`))
ntiss <- length(tissue_codes)

hmp_intercepts <- lapply(param_names, function(pn){
  prob_sign_diff <- pmean_diff <- diag(ntiss) * 0
  for(i in 1:(length(tissue_codes)-1)){
    cat(paste0(i, " "))
    tiss_i <- tissue_codes[i]
    vals_i <- intercepts[[tiss_i]][[pn]]
    for(j in (i+1):length(tissue_codes)){
      tiss_j <- tissue_codes[j]
      vals_j <- intercepts[[tiss_j]][[pn]]
      diff_vals <- unlist(vals_i - vals_j)
      pmean_diff[i,j] <- mean(diff_vals)
      prob_sign_diff[i,j] <- mean(diff_vals > 0)
      if(sign(pmean_diff[i,j]) == -1){
        prob_sign_diff[i,j] <- 1 - prob_sign_diff[i,j]   
      }
    }
  }
  prob_sign_diff <- prob_sign_diff + t(prob_sign_diff)
  pmean_diff <- pmean_diff + t(pmean_diff)
  
  return(list(pmean_diff = pmean_diff,
              prob_sign_diff = prob_sign_diff))
})

# hist(hmp_intercepts$conc_base$pmean_diff[upper.tri(diag(ntiss))], breaks = 100)
# hist(hmp_intercepts$conc_base$prob_sign_diff[upper.tri(diag(ntiss))], breaks = 100)
# plot(hmp_intercepts$conc_base$pmean_diff[upper.tri(diag(ntiss))], 
#      hmp_intercepts$conc_base$prob_sign_diff[upper.tri(diag(ntiss))], 
#      pch = 19)

# for(pn in param_names){
#   hist(hmp_intercepts[[pn]]$prob_sign_diff[upper.tri(diag(ntiss))], breaks = 100, main = pn)
# }
# 
# for(pn in param_names){
#   hist(hmp_intercepts[[pn]]$pmean_diff[upper.tri(diag(ntiss))], breaks = 100, main = pn)
# }

# hmp_intercepts[[pn]]$prob_sign_diff

#### heatmap of intercept diffs ####
source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")

pn <- param_names[5]
tiss_order <- order(sapply(intercepts, function(x) mean(unlist(x[[pn]]))))
pmean_diff <- hmp_intercepts[[pn]]$pmean_diff[tiss_order, tiss_order] #matrix for cell colors
prob_sign_diff <- hmp_intercepts[[pn]]$prob_sign_diff[tiss_order, tiss_order] #matrix for cell sizes
tiss_names <- tissue_codes[tiss_order]
diag(prob_sign_diff) <- 0.5
# my_heatmap(mat_cols = pmean_diff, mat_sizes = prob_sign_diff, dim_names = tiss_names)



#### GTEx intercept ridge plots ####
source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")
param_names <- setNm(names(intercepts$`adipose-visceral-(omentum)`))

#specify plot metadata
fancy_names <- list(
  loc_base = "location of ASE",
  conc_base = "log(concentration) residual ASE",
  prob_AB_base = "Pr(AB)",
  cond_prob_loh_base = "Pr(LoH | ¬AB)",
  prob_loh_base = "Pr(LoH)"
)
params_to_plot <- param_names[-c(1,5)][c(1,3,2)]
space_betw_plots <- 1.5
decr_vec <- c(F, T, F)
reorder_each_panel <- T

#initialize graphical device
grDevices::cairo_pdf(filename = paste0("figures/multitiss-ridgelines", ifelse(reorder_each_panel, "_reordered", "_same-order"), ".pdf"), width = 1000 / 72, height = 600 / 72)
par(xpd = NA, mar = c(3, 14, 3, 0), fig = c(0, 1, 0, 1))
plot.new()
plot.window(xlim = c(0,(length(params_to_plot)-1) * space_betw_plots + 1), 
            ylim = c(1, length(intercepts)))

#plot ridges and retain info on bounds
bounds <- list()
for(i in seq_along(params_to_plot)){
  pn <- params_to_plot[i]
  d <- lapply(intercepts, function(x) setNames(unlist(x[[pn]]), rep("", nrow(x[[pn]]))))
  hard_bounds_range <- quantile(sapply(d, quantile, c(0.01, 0.99)), c(0.01, 0.99))
  bounds[[i]] <- my_ridgeline(d = d, scale_vertical = 0.4, kern_adj = 2, 
                              decreasing = decr_vec[i], xlab = fancy_names[pn], 
                              space_scale_ylabs = 1.2, hard_bounds_range = hard_bounds_range,
                              ylabs = ifelse(i == 1, T, F),
                              panel_coords = c(c(0,1) + i * space_betw_plots - space_betw_plots, 
                                               c(0,length(tiss_order))),
                              order_inputs = ifelse2(i==1 || reorder_each_panel, 
                                                     NULL, 
                                                     match(bounds[[1]]$name, names(d))))
}

#connect ridges with bounds
connect_ridges(bounds)

dev.off()

#### heatmap of tissue corrs x genes #### 
tiss_names <- setNm(names(genefits))
param_names <- setNm(names(genefits[[1]]))
nsamp <- nrow(genefits[[1]][[1]])
tiss_genes <- lapply(genefits, function(x) colnames(x[[1]]))
all_genes <- setNm(unique(unlist(tiss_genes)))
shared_genes <- setNm(Reduce(intersect, tiss_genes))
n_samps_to_use <- 10
samp_inds <- sort(sample(nsamp, n_samps_to_use))
mean_of_corrmat <- T
genewise_corrmats <- lapply(param_names, function(pn){
  print(pn)
  
  if(mean_of_corrmat){
    #calculate mean of corrmats
    sample_corrmats <- parallel::mclapply(samp_inds, function(i){
      gene_x_tiss_iter <- t(as.matrix(do.call(rbind, lapply(tiss_names, function(tn){
        genefits_sub <- genefits[[tn]][[pn]][1,]
        missing_genes <- all_genes[!(all_genes %in% colnames(genefits_sub))]
        genefits_sub <- cbind(genefits_sub,
                              t(as.data.frame(setNames(rep(NA, length(missing_genes)), missing_genes))))
        genefits_sub[,..all_genes]
      }))))
      colnames(gene_x_tiss_iter) <- tiss_names
      sample_cormat <- cor(gene_x_tiss_iter, use = "pairwi")
      return(sample_cormat)
    }, mc.cores = 4)
    out_corrmat <- Reduce("+", sample_corrmats) / n_samps_to_use
  } else {
    #calculate corrmat of means
    gene_x_tiss_iter <- t(as.matrix(do.call(rbind, parallel::mclapply(tiss_names, function(tn){
      genefits_sub <- apply(genefits[[tn]][[pn]], 2, mean)
      missing_genes <- all_genes[!(all_genes %in% colnames(genefits_sub))]
      genefits_sub <- c(genefits_sub, setNames(rep(NA, length(missing_genes)), missing_genes))
      genefits_sub[all_genes]
    }, mc.cores = 4))))
    colnames(gene_x_tiss_iter) <- tiss_names
    out_cormat <- cor(gene_x_tiss_iter, use = "pairwi")
    return(out_corrmat)
  }
})
sapply(genewise_corrmats, function(x) mean(abs(x[upper.tri(x)])))

#### heatmaps plots ####

source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")
for(i in seq_along(genewise_corrmats)){
  pn <- names(genewise_corrmats)[i]; print(pn)
  grDevices::cairo_pdf(filename = paste0("figures/gene_corrmats/", pn,".pdf"), 
                       width = 675 / 72, height = 500 / 72)
  my_heatmap(genewise_corrmats[[i]], reorder_mat = T, plot_diagonal_labels = F, symmetric_cols = F)
  dev.off()
}
pdftools::pdf_combine(input = paste0("figures/gene_corrmats/", names(genewise_corrmats),".pdf"),
                        output = "figures/all_gene_corrmats_combined.pdf")

#### MDS / PCOA plot? #####
tiss_dists <- sqrt(Reduce("+", lapply(param_names, function(pn) hmp_intercepts[[pn]]$pmean_diff^2 / 
         sd((hmp_intercepts[[pn]]$pmean_diff^2)[upper.tri(diag(ntiss))]))) / length(param_names))

MDSpts <- PCoApts <- cmdscale(tiss_dists)
stdMDSpts <- smacof::mds(tiss_dists, init = MDSpts, itmax = 1E6, eps = 1E-10)$conf

tissue_categories <- c(
  "adipose-subcutaneous" = "Adipose",
  "adipose-visceral-(omentum)" = "Adipose",
  "adrenal-gland" = "Endocrine",
  "artery-aorta" = "Cardiovascular",
  "artery-coronary" = "Cardiovascular",
  "artery-tibial" = "Cardiovascular",
  "bladder" = "Urinary",
  "brain-amygdala" = "Brain",
  "brain-anterior-cingulate-cortex-(ba24)" = "Brain",
  "brain-caudate-(basal-ganglia)" = "Brain",
  "brain-cerebellar-hemisphere" = "Brain",
  "brain-cerebellum" = "Brain",
  "brain-cortex" = "Brain",
  "brain-frontal-cortex-(ba9)" = "Brain",
  "brain-hippocampus" = "Brain",
  "brain-hypothalamus" = "Brain",
  "brain-nucleus-accumbens-(basal-ganglia)" = "Brain",
  "brain-putamen-(basal-ganglia)" = "Brain",
  "brain-spinal-cord-(cervical-c-1)" = "Brain",
  "brain-substantia-nigra" = "Brain",
  "breast-mammary-tissue" = "Breast",
  "colon-sigmoid" = "Digestive",
  "colon-transverse" = "Digestive",
  "esophagus-gastroesophageal-junction" = "Digestive",
  "esophagus-mucosa" = "Digestive",
  "esophagus-muscularis" = "Digestive",
  "fallopian-tube" = "Reproductive",
  "heart-atrial-appendage" = "Cardiovascular",
  "heart-left-ventricle" = "Cardiovascular",
  "kidney-cortex" = "Urinary",
  "kidney-medulla" = "Urinary",
  "liver" = "Digestive",
  "lung" = "Respiratory",
  "minor-salivary-gland" = "Digestive",
  "muscle-skeletal" = "Muscular",
  "nerve-tibial" = "Brain",
  "ovary" = "Reproductive",
  "pancreas" = "Endocrine",
  "pituitary" = "Endocrine",
  "prostate" = "Reproductive",
  "skin-not-sun-exposed-(suprapubic)" = "Skin",
  "skin-sun-exposed-(lower-leg)" = "Skin",
  "small-intestine-terminal-ileum" = "Digestive",
  "spleen" = "Immune",
  "stomach" = "Digestive",
  "testis" = "Reproductive",
  "thyroid" = "Endocrine",
  "uterus" = "Reproductive",
  "vagina" = "Reproductive",
  "whole-blood" = "Blood"
)


cols <- c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(12, "Set3"))[1:length(unique(tissue_categories))]
cols <- setNames(cols, unique(tissue_categories))
all(names(tissue_categories) %in% tissue_codes)
all(tissue_codes %in% names(tissue_categories))

par(mar = c(0,0,2,5), xpd = T)
dist_thresh <- quantile(tiss_dists[upper.tri(tiss_dists)], 0.2)
edges_to_draw <- which((tiss_dists) < dist_thresh, arr.ind = T)
edges_to_draw <- edges_to_draw[apply(edges_to_draw, 1, diff) != 0,]

corr_edge_col <- 1
corr_edge_weight <- 1
plot(1,1,xlim = range(PCoApts[,1]), ylim = range(PCoApts[,2]), 
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
segments(x0 = PCoApts[edges_to_draw[,1],1], y0 = PCoApts[edges_to_draw[,1],2], 
         x1 = PCoApts[edges_to_draw[,2],1], y1 = PCoApts[edges_to_draw[,2],2], 
         col = corr_edge_col, lwd = corr_edge_weight, xpd = NA)
points(x = PCoApts[,1], y = PCoApts[,2], pch = 19, xpd = NA, col = cols[tissue_categories[tissue_codes]], cex = 2)
text(x = PCoApts[,1], y = PCoApts[,2], labels = tissue_codes, xpd = NA, col = cols[tissue_categories[tissue_codes]])

plot(NULL,xlim = range(stdMDSpts[,1]), ylim = range(stdMDSpts[,2]), 
     col = "white", xaxt = "n", yaxt = "n", frame.plot = FALSE, xlab = "", ylab = "")
segments(x0 = stdMDSpts[edges_to_draw[,1],1], y0 = stdMDSpts[edges_to_draw[,1],2], 
         x1 = stdMDSpts[edges_to_draw[,2],1], y1 = stdMDSpts[edges_to_draw[,2],2], 
         col = corr_edge_col, lwd = corr_edge_weight, xpd = NA)
points(x = stdMDSpts[,1], y = stdMDSpts[,2], pch = 19, xpd = NA, 
       col = cols[tissue_categories[tissue_codes]], cex = 2)
text(x = stdMDSpts[,1], y = stdMDSpts[,2], labels = tissue_codes, xpd = NA, col = cols[tissue_categories[tissue_codes]])


#### compare rnaseq to mmpcrseq ####

#load in the processed model fits
if(!exists("focal_samps_rnaseq")){
  load(file = paste0(stan_output_dir, "focal-samps_", model_name, "_", sample_prop_tail_to_remove, 
                     ifelse(init_from_prior, "_prior-init", "_reg-init"),
                     ifelse(use_EB,
                            ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                            "_FB"), 
                     "_rnaseq.RData"))
  focal_samps_rnaseq <- focal_samps
  rm(focal_samps)
}

if(!exists("focal_samps_mmpcrseq")){
  load(file = paste0(stan_output_dir, "focal-samps_", model_name, "_", sample_prop_tail_to_remove, 
                     ifelse(init_from_prior, "_prior-init", "_reg-init"),
                     ifelse(use_EB,
                            ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                            "_FB"), 
                     ".RData"))
  focal_samps_mmpcrseq <- focal_samps
  rm(focal_samps)
}

param_of_interest <- "logodds_AB_gene_abs_dev"
tissue_of_interest <- names(focal_samps_rnaseq$genefits)[1]
tiss_names <- names(focal_samps_rnaseq$genefits)

posterior_means <- list(rnaseq = apply(focal_samps_rnaseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, mean), 
                        mmpcrseq = apply(focal_samps_mmpcrseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, mean))
shared_genes <- intersect(names(posterior_means$rnaseq), names(posterior_means$mmpcrseq))
plot(posterior_means$rnaseq[shared_genes], posterior_means$mmpcrseq[shared_genes])
cor(posterior_means$rnaseq[shared_genes], posterior_means$mmpcrseq[shared_genes])

posterior_prec_diff <- parallel::mclapply(setNm(tiss_names), function(tiss_j){
  mmpcrsamps <- focal_samps_mmpcrseq$genefits[[tiss_j]][[param_of_interest]]
  rnasamps <-  focal_samps_rnaseq$genefits[[tiss_j]][[param_of_interest]]
  pprec_mmpcrsamps <- 1 / apply(mmpcrsamps, 2, var)
  pprec_rnasamps <- 1 / apply(rnasamps, 2, var)
  shared_genes <- intersect(names(pprec_mmpcrsamps), names(pprec_rnasamps))
  prec_diff <- pprec_mmpcrsamps[shared_genes] - pprec_rnasamps[shared_genes]
  prec_diff
}, mc.cores = 4)

#####

grDevices::cairo_pdf(filename = paste0("figures/posterior_precision_difference_in_", param_of_interest, ".pdf"), 
                     width = 500 / 72, height = 500 / 72)
par(xpd = NA, mar = c(3, 14, 3, 0), fig = c(0, 1, 0, 1))
plot.new()
plot.window(xlim = c(0,1), 
            ylim = c(1, length(posterior_prec_diff)))
my_ridgeline(posterior_prec_diff, scale_vertical = 0.25, space_scale_ylabs = 1.2,
             panel_coords = c(c(0.2,1), c(0,length(posterior_prec_diff))), resc_to_maxh = T)

dev.off()
  
#####

list(rnaseq = apply(focal_samps_rnaseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, var), 
                        mmpcrseq = apply(focal_samps_mmpcrseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, var))


tissue_pm_corrs <- parallel::mclapply(setNm(names(focal_samps_rnaseq$genefits[[1]])), function(param_of_interest){
  mcprint(paste0("(", param_of_interest, ")"))
  sapply(setNm(names(focal_samps_rnaseq$genefits)), function(tissue_of_interest){
    posterior_means <- list(rnaseq = apply(focal_samps_rnaseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, mean), 
                            mmpcrseq = apply(focal_samps_mmpcrseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, mean))
    shared_genes <- intersect(names(posterior_means$rnaseq), names(posterior_means$mmpcrseq))
    cor(posterior_means$rnaseq[shared_genes], posterior_means$mmpcrseq[shared_genes])
  })
}, mc.cores = 4)

par(mfrow = c(3,3))
for(i in names(tissue_pm_corrs)){
  hist(tissue_pm_corrs[[i]], main = i)
  abline(v = 0, col = 2, lwd = 2)
}

hist(tissue_pm_corrs[["logodds_AB_gene"]], 
     main = latex2exp::TeX("Pearson's \\textit{r} in Gene-wise Posterior Means for \\textbf{log-odds(Allelic Balance)} across Tissues between RNASeq vs mmPCRSeq"), 
     freq = F, xlab = "Pearson's Correlation Coefficient", cex.main = 0.85)

tissue_posterior_vars <- parallel::mclapply(setNm(names(focal_samps_rnaseq$genefits[[1]])), function(param_of_interest){
  lapply(setNm(names(focal_samps_rnaseq$genefits)), function(tissue_of_interest){
    posterior_vars <- list(rnaseq = apply(focal_samps_rnaseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, var), 
                            mmpcrseq = apply(focal_samps_mmpcrseq$genefits[[tissue_of_interest]][[param_of_interest]], 2, var))
    shared_genes <- intersect(names(posterior_vars$rnaseq), names(posterior_vars$mmpcrseq))
    return(data.frame(rnaseq = unlist(posterior_vars$rnaseq[shared_genes]), 
                      mmpcrseq = unlist(posterior_vars$mmpcrseq[shared_genes])))
  })
}, mc.cores = 8)

par(mfrow = c(3,3))
for(i in names(tissue_pm_corrs)){
  plot(do.call(rbind, tissue_posterior_vars[[i]]), main = i, col = adjustcolor(1, 0.1), pch = 19)
  abline(a = 0, b = 1, col = 2, lwd = 2)
}

par(mfrow = c(1,1), mar = c(5,5,6,8))
var1_gr_var2_prop <- mean(apply(do.call(rbind, tissue_posterior_vars$logodds_AB_gene), 1, diff) < 0)
plot(do.call(rbind, tissue_posterior_vars$logodds_AB_gene), 
     main = latex2exp::TeX("Variance of Marginal Posterior Distributions of Gene-wise \\textbf{log-odds(Allelic Balance)}"), 
     cex.main = 1.25,
     col = adjustcolor(1, 0.1), pch = 19, xlab = "RNASeq", ylab = "mmPCRSeq")
text(x = par("usr")[2], y = par("usr")[4] - diff(par("usr")[3:4])/50, 
     labels = paste0("Pr(mmPCRSeq < RNASeq) ≈ ", round(var1_gr_var2_prop, 3)), pos = 2)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
legend(x = par("usr")[2] + diff(par("usr")[1:2])/100, y = par("usr")[4], col = c(adjustcolor(1, 0.1), 2), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
       legend = c("same gene", "1-to-1 line"), xpd = NA,bty = "n")


plot(do.call(rbind, lapply(tissue_posterior_vars$logodds_AB_gene, function(x) apply(x, 2, function(y) exp(mean(log(y)))))), 
     main = latex2exp::TeX("Tissue-Specific Geometric Mean Variance\nof Marginal Posterior Distributions of Gene-wise \\textbf{log-odds(Allelic Balance)}"), 
     cex.main = 1.25,
     col = adjustcolor(1, 0.5), pch = 19, xlab = "RNASeq", ylab = "mmPCRSeq")
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
legend(x = par("usr")[2] + diff(par("usr")[1:2])/100, y = par("usr")[4], col = c(adjustcolor(1, 0.1), 2), lwd = c(NA, 2), lty = c(NA, 2), pch = c(19, NA), 
       legend = c("same gene", "1-to-1 line"), xpd = NA,bty = "n")



#some examples for a random tissue (1st alphabetically), adipose

hist(genefits$`adipose-subcutaneous`$conc_gene$JAK1)
hist(apply(genefits$`adipose-subcutaneous`$conc_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$conc_gene, 1, sd))
hist(apply(genefits$`adipose-subcutaneous`$abs_conc_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$abs_conc_gene, 2, sd))

hist(genefits$`adipose-subcutaneous`$loc_gene$JAK1)
hist(apply(genefits$`adipose-subcutaneous`$loc_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$loc_gene, 2, sd))
hist(apply(genefits$`adipose-subcutaneous`$abs_loc_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$abs_loc_gene, 2, sd))

hist(genefits$`adipose-subcutaneous`$logodds_AB_gene$JAK1)
hist(apply(genefits$`adipose-subcutaneous`$logodds_AB_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$logodds_AB_gene, 2, sd))
hist(apply(genefits$`adipose-subcutaneous`$abs_prob_AB_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$abs_prob_AB_gene, 2, sd))

hist(genefits$`adipose-subcutaneous`$cond_logodds_loh_gene$JAK1)
hist(apply(genefits$`adipose-subcutaneous`$cond_logodds_loh_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$cond_logodds_loh_gene, 2, sd))
hist(apply(genefits$`adipose-subcutaneous`$abs_prob_loh_gene, 2, mean))
hist(apply(genefits$`adipose-subcutaneous`$abs_prob_loh_gene, 2, sd))


hist(samps$sd_conc_gene)
hist(apply(samps_by_param$conc_gene, 2, sd))


#### generate intercepts comparison figures ####
# dev.off()

#high level params
target_param <- names(intercepts$`adipose-subcutaneous`)[4]

#density plot
all_cols <- c('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499', '#DDDDDD')
param_ranges <- do.call(rbind, lapply(intercepts, function(xint) c(min(xint[[target_param]]), quantile(xint[[target_param]], 0.999))))
params_range <- c(min(param_ranges[,1]), max(param_ranges[,2]))
densities <- lapply(intercepts, function(xint) density(xint[[target_param]], from = params_range[1], to = params_range[2]))
density_maxima <- sapply(densities, function(dens) max(dens$y))
density_maxima_locs <- sapply(densities, function(dens) dens$x[which.max(dens$y)])
par(mar = c(0,6,18,1))
layout(mat = matrix(1:2, nrow=2), height = c(4,2))
label_locs_x <- seq(from = params_range[1], to = params_range[2], length.out = length(intercepts))
label_loc_y <- max(density_maxima) * 1.05
label_order <- rank(density_maxima_locs, ties.method = "first")
tiss_cols <- all_cols[c(1:5, rep(10, length(label_order)-9), 6:9)]
plot(x = NULL, y = NULL, xlim = params_range, ylim = c(0, max(density_maxima)), ylab = "density", xlab = target_param, frame = F)
for(i in seq_along(densities)){
  polygon(x = c(densities[[i]]$x, rev(densities[[i]]$x)), 
          y = c(densities[[i]]$y, rep(0, length(densities[[i]]$y))), 
          col = adjustcolor(tiss_cols[label_order[i]], 0.05), xpd = NA, border = tiss_cols[label_order[i]])
  
  text(labels = names(densities)[i], x = label_locs_x[label_order[i]] - strheight(names(densities)[i]) / 1.5 * xyrat(), 
       y = label_loc_y, srt = 90, density_maxima_locs[i], pos = 4, xpd = NA, col = tiss_cols[label_order[i]])
  segments(x0 = density_maxima_locs[i], x1 = label_locs_x[label_order[i]], 
           y0 = density_maxima[i], y1 = label_loc_y,lty = 2, 
           xpd = NA, col = tiss_cols[label_order[i]])
}

#quantiles plot
par(mar = c(5,6,4,1))
plot(x = NULL, y = NULL, xlim = params_range, ylim = c(0, 1), ylab = "quantile", xlab = target_param, frame = F)
for(i in seq_along(densities)){
  lines(x = densities[[i]]$x[-1] - diff(densities[[i]]$x)/2, 
        y = cumsum(diff(densities[[i]]$x) * (densities[[i]]$y[-1] - diff(densities[[i]]$y)/2)), col = tiss_cols[label_order[i]])
  # text(labels = names(densities)[i], x = label_locs_x[label_order[i]] - strheight(names(densities)[i]) / 1.5 * xyrat(), 
  #      y = label_loc_y, srt = 90, density_maxima_locs[i], pos = 4, xpd = NA)
  # segments(x0 = density_maxima_locs[i], x1 = label_locs_x[label_order[i]], y0 = density_maxima[i], y1 = label_loc_y,lty = 2, xpd = NA)
}

#### generate genesets of interest ####
dev.off()
gene_params_thresh <- c(conc_gene = 0, loc_gene = 0, logodds_AB_gene = 0, cond_logodds_loh_gene = 0)
gene_params_thresh <- gene_params_thresh[-1]


prop_priority_genes <- 0.1
n_priority_genes <- 20
prob_greater_than <- function(x, thresh = 0) mean(x > thresh)
top_gene_effects <- lapply(setNm(names(gene_params_thresh)), function(gpi){
  gene_summaries <- lapply(genefits, function(xgene){
    gene_means <- apply(xgene[[gpi]], 2, mean)
    gene_prob_pos <- apply(xgene[[gpi]], 2, prob_greater_than, thresh = gene_params_thresh[gpi])
    df_genes <- data.frame(gene_means = gene_means, gene_prob_pos = gene_prob_pos, gene = names(gene_means))
    df_genes <- df_genes[order(df_genes$gene_means),]
    return(df_genes)
  })
  n_genes_per_tiss <- sapply(gene_summaries, nrow)
  
  top_genes <- lapply(setNm(names(gene_summaries)), function(gsi){
    n_priority_genes <- floor(n_genes_per_tiss[[gsi]] * prop_priority_genes)
    df_genes <- data.frame(rbind(cbind(gene = head(gene_summaries[[gsi]]$gene, n_priority_genes), dir = "neg"), 
                                 cbind(gene = tail(gene_summaries[[gsi]]$gene, n_priority_genes), dir = "pos")))
    df_genes$tissue <- gsi
    return(df_genes)
  })
  top_genes_df <- do.call(rbind, top_genes)
  top_genes_df$param <- gpi
  return(top_genes_df)
})

multitiss_genes <- lapply(top_gene_effects, function(param_genes){
  return(list(ntiss_pos_top_20 = head(sort(table(param_genes$gene[param_genes$dir == "pos"]), decreasing = T), 10),
              ntiss_neg_top_20 = head(sort(table(param_genes$gene[param_genes$dir == "neg"]), decreasing = T), 10)))
})

topgene_freqs_x_tissue <- list(
  # loc = sort(table(top_gene_effects$loc_gene$gene[top_gene_effects$loc_gene$dir == "pos"])),
  logodds_AB = sort(table(top_gene_effects$logodds_AB_gene$gene[top_gene_effects$logodds_AB_gene$dir == "pos"])),
  cond_logodds_loh = sort(table(top_gene_effects$cond_logodds_loh_gene$gene[top_gene_effects$cond_logodds_loh_gene$dir == "pos"]))
)

gene_tiss_freqs_df <- data.frame(gene = unlist(sapply(topgene_freqs_x_tissue, names)), freq = unlist(topgene_freqs_x_tissue))
top_n_hits <- sapply(split(gene_tiss_freqs_df$freq, gene_tiss_freqs_df$gene), max)
sum(top_n_hits >= 20)
sum(top_n_hits <= 5)
table(top_n_hits)

all_top_genes <- unique(unlist(sapply(topgene_freqs_x_tissue, names)))

#ERRRROOOOR
#gotta still apply these transforms:

# //likewise coerce the shapes of the allelic balance and loh beta
# real shapes_ab = (20 + exp(6 + log_conc_allelic_balance)) * 0.5;
# real loc_loh = 0.95 + inv_logit(logit_loc_loh) / 20;
# real shape_loh_alpha = (20 + exp(6 + log_conc_loh)) * loc_loh;
# real shape_loh_beta = (20 + exp(6 + log_conc_loh)) * (1-loc_loh);

#### GTEx intercepts heatmaps ####


#### inspect individual GTEx output ####
tiss_j <- "brain-caudate-(basal-ganglia)"
tiss_j <- "lung"

#useful variables
file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                    ifelse(init_from_prior, "_prior-init", "_reg-init"),
                    ifelse(use_EB,
                           ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                           "_FB"))

#load in MCMC output
load(paste0(stan_output_dir, file_name, ".cmdStanR.fit"))
nchains <- out$num_chains()
print(nchains)
load(paste0(stan_output_dir, file_name, ".cmdStanR.summ"))

#### inspect EB from GTEx ####

if(use_EB){
  
  #inspect output briefly to eg assess sensitivity to meanfield init
  use_flat_MFVB <- F
  file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove)
  load(paste0(stan_output_dir, file_name, "_MF-VarInf-Runs.RData"))
  
  if(!use_flat_MFVB){
    lps <- unlist(lapply(seq_along(var_outs), function(i) var_outs[[i]]$summary(variables = "lp_approx__")$mean))
    model_params <- retrieve_block_objects(model_lines, block = "parameters")$var_name
    model_scale_params <- model_params[grepl("sd_", model_params)]
    
    sd_means <- data.frame(do.call(rbind, lapply(seq_along(var_outs), 
                                                 function(i) var_outs[[i]]$summary(model_scale_params)$mean)))
    colnames(sd_means) <- model_scale_params
    cbind(sd_means, lps)
    apply(sd_means, 2, mean)
  } else {
    #hrmph, this also is getting stuck at the optima!
    #can try to fit the flattened model and then summarize?
    
    #extract target param samples
    
    sd_means <- data.frame(do.call(rbind, lapply(seq_along(var_outs), function(i) {
      samps_mfvb <- as.data.table(data.frame(as_draws_df(var_outs[[i]]$draws(model_scale_param_targets))))[1:1000,]
      param_mean_var_mfvb <- sapply(model_scale_param_targets, function(scale_param){
        param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_mfvb)))
        param_vars <- apply(param_samps, 1, var)
        mean(param_vars)
      })
      param_sds <- sqrt(param_mean_var_mfvb)
      param_sds
    }))
    )
    colnames(sd_means) <- model_scale_params
    cbind(sd_means, lps)
  }
  
  #PROBLEM here is that the start values are not valid for the truncated normal parameters?
  
}

#### inspect MCMC from GTEx ####
tiss_j <- "lung"

#useful variables
file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                    ifelse(init_from_prior, "_prior-init", "_reg-init"),
                    ifelse(use_EB,
                           ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                           "_FB"))

#load in MCMC output
load(paste0(stan_output_dir, file_name, ".cmdStanR.fit"))
nchains <- out$num_chains()
print(nchains)
load(paste0(stan_output_dir, file_name, ".cmdStanR.summ"))

#check output
summ[order(summ$ess_bulk), c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]
summ[order(summ$rhat, decreasing = T),c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]


#when <out> does not work as anticipated... except of course the chain ids are not
#in the outputted object and I have no way to reconstruct them :/
# outmtdt <- out$metadata()
# csv_files <- list.files(stan_output_dir, pattern = "\\.csv$", full.names = TRUE)
# csv_files <- csv_files[order(file.info(csv_files)$mtime, decreasing = T)][1:12]
# # csv_files <- csv_files[grepl(file_name, csv_files)]
# csv_file_info <- do.call(rbind, strsplit(csv_files, "-"))
# csv_file_info <- split(as.data.frame(csv_file_info), csv_file_info[,6])
# csv_files <- csv_files[grepl(names(csv_file_info)[4], csv_files)]
# out <- as_cmdstan_fit(csv_files)
# summ <- out$summary()

#load in data
dt <- fread(file = paste0(stan_data_dir, tiss_j, "_",  sample_prop_tail_to_remove, ".csv"))
dat <- jsonlite::fromJSON(paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, ".json"))
# dt_removed <- fread(file = paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_removed.csv"))
# dt_dupes <- rbind(dt[dt$gene_indiv %in% dt_removed$gene_indiv,], dt_removed)
genes_tiss <- fread(file = paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_gene-map.csv"))
indivs_tiss <- fread(file = paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_indiv-map.csv"))

#check diagnostics and extract samples
summ[order(summ$ess_bulk), c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]
summ[order(summ$rhat, decreasing = T),c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]
samps <- data.frame(as_draws_df(out$draws()))
plot(samps$lp__, type = "l")
plot(samps$logit_loc_base, type = "l")
plot(samps$log_conc_base, type = "l")
plot(samps$sd_logodds_AB_gene, type = "l")
plot(samps$sd_conc_gene, type = "l")

# plot(samps$sd, type = "l")

#compare eg gene #128 or 75 vs gene #160

gene_index <- 185
obs_inds <- which(dt$gene == genes_tiss$gene_symbol[genes_tiss$index == gene_index])
subdt <- as.data.frame(dt[obs_inds,])
hist(subdt$refCount / subdt$totalCount, breaks = 0:100/100)

plot(samps[,paste0("conc_gene.", gene_index, ".")], type = "l")
pairs(list(gene_conc = samps[,paste0("conc_gene.", gene_index, ".")], 
           log_conc_base = samps$log_conc_base, sd_conc_gene = samps$sd_conc_gene), 
      col = adjustcolor(1, 0.02), pch = 19)

plot(samps[,paste0("loc_gene.", gene_index, ".")], type = "l")
pairs(list(gene_loc = samps[,paste0("loc_gene.", gene_index, ".")], 
           logit_loc_base = samps$logit_loc_base, sd_loc_gene = samps$sd_loc_gene), 
      col = adjustcolor(1, 0.02), pch = 19)

plot(samps[,paste0("logodds_AB_gene.", gene_index, ".")], type = "l")
pairs(list(samps[,paste0("logodds_AB_gene.", gene_index, ".")], 
           samps$logodds_AB_base, samps$sd_logodds_AB_gene), 
      col = adjustcolor(1, 0.02), pch = 19)

pairs(list(conc = samps[,paste0("conc_gene.", gene_index, ".")], 
           loc = samps[,paste0("loc_gene.", gene_index, ".")], 
           logodds_AB = samps[,paste0("logodds_AB_gene.", gene_index, ".")], 
           cond_logodds_LoH = samps[,paste0("cond_logodds_loh_gene.", gene_index, ".")]), 
      col = adjustcolor(1, 0.02), pch = 19)


#process Stan model to evaluate posterior predictive densities / massess
#with a focus on these genes
source("~/repos/Stan2R/R/functions.R")
stan_code <- model_lines <- readLines(model_path)
dat <- dat
chains_to_compare <- c(1,3)
nsamps <- nrow(samps)
n_per_chain <- 100
sample_indices_list <- lapply(chains_to_compare, function(ci){
  sort(sample((nsamps / nchains * (ci-1) + 1):(nsamps / nchains * ci), size = n_per_chain, replace = F))
})
sample_indices <- unlist(sample_indices_list)
ppfits <- parallel::mclapply(seq_along(sample_indices), function(i){
  
  mcprint(paste0(i, " "))
  sample_index <- sample_indices[i]
  r_code <- parse_Stan(stan_code = stan_code, 
                       dat = dat, 
                       samps = as.data.table(samps), 
                       output_file = NA, 
                       sample_index = sample_index, 
                       post_pred_sim = T, 
                       sim = F)
  # open_script(r_code)
  
  stan_env <- new.env()
  eval(parse(text = r_code), envir = stan_env)
  posterior_predictive_mass <- stan_env$out
  
  return(posterior_predictive_mass)
}, mc.cores = 8)

n_iter_per_chain <- sapply(sample_indices_list, length)
iter_per_chain <- lapply(seq_along(n_iter_per_chain), function(i) 1:n_iter_per_chain[i] + ifelse(i != 1, sum(n_iter_per_chain[1:(i-1)]), 0))
sim_varnames <- setNames(names(ppfits[[1]]), names(ppfits[[1]]))
ppfits_per_chain_means <- lapply(iter_per_chain, function(eyes){
  sum_of_params <- lapply(sim_varnames, function(varname){
    Reduce("+", lapply(eyes, function(i) {
      ppfits[[i]][[varname]]
    })) / length(eyes)
  })
})

pairs(rbind(ppfits_per_chain_means[[1]]$mix_prob[obs_inds,1:3*2-1], 
            ppfits_per_chain_means[[2]]$mix_prob[obs_inds,1:3*2-1]),
      col = rep(c(1,2), each = length(obs_inds)), 
      labels = c("Pr(AB)", "Pr(LoH)", "Pr(Slab)"))

pairs(rbind(cbind(ppfits_per_chain_means[[1]]$mix_prob[obs_inds,1:3*2-1], 
                  ppfits_per_chain_means[[1]]$loc[obs_inds], 
                  ppfits_per_chain_means[[1]]$conc[obs_inds]), 
            cbind(ppfits_per_chain_means[[2]]$mix_prob[obs_inds,1:3*2-1], 
                  ppfits_per_chain_means[[2]]$loc[obs_inds], 
                  ppfits_per_chain_means[[2]]$conc[obs_inds])),
      col = rep(c(1,2), each = length(obs_inds)), 
      labels = c("Pr(AB)", "Pr(LoH)", "Pr(Slab)", "Slab Loc", "Slab Conc"))

#probability of being in AB is high when slab loc is 0.9, but low when slab loc is 0.55

#OK, new problem, slab is piling mass near 0 and 1 with absurd concentrations eg 50000, which steals from Pr(LoH)
#will try harder bounds on location first and leave concentration untouched -- otherwise may have to swap an invlogit in for the exp

hist(subdt$refCount / subdt$totalCount, breaks = 0:100/100)
(20 + exp(6 + ppfits[[200]]$log_conc_allelic_balance))
(20 + exp(6 + ppfits[[200]]$log_conc_loh))

ppfits[[1]]$shapes_ab
ppfits[[200]]$shapes_ab
ppfits[[1]]$shape_loh_alpha
ppfits[[200]]$shape_loh_alpha
ppfits[[1]]$shape_loh_beta
ppfits[[200]]$shape_loh_beta



plot(ppfits[[1]]$count[obs_inds],
     ppfits[[200]]$count[obs_inds]); abline(0,1)

ppfits[[1]]$mix_prob[obs_inds,1]

iter_ind <- 200
cols2use <- colgrad(ppfits[[iter_ind]]$count[obs_inds])
orderobs <- order(ppfits[[iter_ind]]$count[obs_inds])
par(mar = c(5,5,3,6))
cexes <- sqrt(subdt$totalCount / max(subdt$totalCount) * 8)
plot(subdt$refCount / subdt$totalCount, 
     ppfits[[iter_ind]]$mix_prob[obs_inds,1],
     col = adjustcolor(cols2use, 0.75), pch = 19, cex = cexes)
# text(subdt$refCount / subdt$totalCount, 
#      ppfits[[iter_ind]]$mix_prob[obs_inds,1],
#      col = adjustcolor(cols2use, 0.75), pch = 19, cex = cexes)
add_continuous_legend(colors = cols2use[orderobs],  scientific = F,
                      labels = ppfits[[iter_ind]]$count[obs_inds][orderobs], 
                      positions = resc01(ppfits[[iter_ind]]$count[obs_inds])[orderobs], 
                      x = par("usr")[2] + diff(par("usr")[1:2])/50, 
                      y = par("usr")[4], xpd = NA, left_below = F)





#PROBLEM!

#my allelic balance coefficient concentration is 1 (exp(6+ppfits[[200]]$log_conc_allelic_balance)), 
#which puts mass near 0 and 1 -- the LoH condition! with low probability

#my slab location is very tightly clustered around 0.5 for all most genes: hist(ppfits[[200]]$loc, breaks = 100)

#LoH has moderate probability, but not high

#these are the opposite of the dynamics I want -- need better enforcement, maybe via
# 1) priors (check prior predictive on component probabilities and locations? maybe it is doing poorly)
# 2) explicit liabilities with truncated constraints, maybe with adaptive sampling of observation probs? only adds dat$n = 7152 more parameters
# 3) modeling concentration as increasing near boundaries?

#can do 1 but then check 2 for misbehavior


pairs(rbind(cbind(loc = ppfits[[1]]$loc[obs_inds], 
                  conc = ppfits[[1]]$conc[obs_inds], 
                  slab_prob = ppfits[[1]]$mix_prob[obs_inds,5]), 
            cbind(loc = ppfits[[200]]$loc[obs_inds], 
                  conc = ppfits[[200]]$conc[obs_inds], 
                  slab_prob = ppfits[[200]]$mix_prob[obs_inds,5])),
      col = rep(c(1,2), each = length(obs_inds)))


#model-based outlier calling?

#summarize flat fit scale empirically

#### inspect MCMC from OVC-vs-GTEx ####
stan_model_dir <- "Stan/models/"
stan_data_dir <- "Stan/data/"
stan_output_dir <- "Stan/output/"
stan_progress_dir <- "Stan/progress/"
tiss_j <- tissue_code <- "ovary"

#load everything isample_prop_tail_to_remove <- 0

#first find model name
group_index <- 2
model_index <- 9
model_names <- list(
  list(
    "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene", #1
    "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid", #2
    "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid-identifiability", #3
    "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid_nearly-balanced", #4
    "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_cumulative_nearly-balanced" #5
  ),
  list(
    "ase-gtex_inf-priors_no-indiv", #1
    "ase-gtex_inf-priors", #2
    "ase-gtex", #3
    "ase-gtex_centered", #4
    "ase-gtex_offset-conc", #5
    "ase-gtex_nearly-balanced", #6
    "ase-gtex_nearly-balanced_LoH", #7
    "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade", #8
    "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc", #9
    "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc-correct_grade_REfs" #10
  )
)
model_name <- model_names[[group_index]][[model_index]]

#then data name
use_rnaseq <- F
data_index <- 4
data_names <- list(
  "OVC-vs-GTEx_0", #1
  "OVC-vs-GTEx_perm_0", #2
  "OVC-vs-GTEx_sim_0", #3
  "OVC-vs-GTEx_0_all", #4, 
  "", #5
  "" #6
)
data_name <- paste0(data_names[[data_index]], ifelse(use_rnaseq, "_rnaseq", ""))


#load in MCMC output
load(paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.fit"))
load(paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.summ"))

#load in data
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
dat <- jsonlite::fromJSON(paste0(stan_data_dir, data_name, ".json"))
d <- fread(paste0(stan_data_dir, "OVC-vs-GTEx", "_", sample_prop_tail_to_remove, ".csv"))
gene_map <- fread(file = paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove, "_gene-map.csv"))
indiv_map <- fread(file = paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove, "_indiv-map.csv"))

#check diagnostics and extract samples
summ[order(summ$ess_bulk), c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]
summ[order(summ$rhat, decreasing = T),c(colnames(summ)[1:2], "rhat", "ess_bulk")][1:20,]
samps <- data.frame(as_draws_df(out$draws()))

#inspect genes w/ failed diagsnostics?
gene_index <- 13
subd <- d[d$gene_i == gene_index,]
hist(subd$refCount / subd$totalCount, breaks = 0:100/100)
par(mfrow = c(length(unique(subd$grade_k)), 1), mar = c(2,4,2,4))
lapply(split(subd$refCount / subd$totalCount, subd$grade_k), hist, breaks = 0:100/100)
table(subd$grade_k)
subd$totalCount
mean(nrow(subd) > sapply(split(d$gene, d$gene_i), length)) #n obs relative to other genes
mean(subd$refCount[subd$grade_k == 1] / subd$totalCount[subd$grade_k == 1] > 0.99)

#some quick plots of interest
hist(samps$logit_loc_base)
hist(0.6 + plogis(samps$logit_loc_base / 2) / 10 * 3)
hist(plogis(samps$logodds_AB_base)) #prob AB
hist((1-plogis(samps$logodds_AB_base)) * plogis(samps$cond_logodds_loh_base)) #prob loh
hist((1-plogis(samps$logodds_AB_base)) * (1-plogis(samps$cond_logodds_loh_base))) #prob ase


hist(samps$mu_logodds_AB_grade)
hist(samps$mu_conc_grade)
hist(samps$mu_cond_logodds_loh_grade)
hist(samps$mu_loc_grade)


hist(samps$logodds_AB_gtex)
hist(samps$conc_gtex)
hist(samps$cond_logodds_loh_gtex)
hist(samps$loc_gtex)
hist(samps$mu_loc_grade)


mean(samps$mu_logodds_AB_grade < 0)
mean(samps$mu_conc_grade < 0)
mean(samps$mu_cond_logodds_loh_grade < 0)
mean(samps$mu_loc_grade < 0)

hist(samps$log_conc_base)
hist(samps$sd_conc_gene)
hist(samps$sd_conc_indiv)
hist(samps$conc_gtex)


plot(samps$lp__, type = "l")
plot(samps$sd_loc_gene, type = "l")
plot(samps$sd_conc_indiv, type = "l")
plot(samps$sd_loc_gene, samps$sd_conc_gene)
plot(samps$sd_loc_indiv, samps$sd_conc_indiv)
plot(samps$sd_loc_grade, samps$sd_conc_grade)
plot(samps$conc_gene.205., samps$loc_gene.205.)
hist(samps$log_conc_allelic_balance)
hist((20 + exp(6 + samps$log_conc_loh)))
hist(plogis(samps$logit_loc_loh)/20+0.95)

# dgtex <- d[d$grade == "GTEx",]
# gene_prop_extreme <- sapply(split(dgtex$refCount / dgtex$totalCount, dgtex$gene), function(x) mean(x > 0.99 | x < 0.01))
# hist(gene_prop_extreme, breaks = 0:100/100)
# quantile(gene_prop_extreme, 0:10/10)
# mean(gene_prop_extreme > 0.5)
# 
# het_prop_extreme <- sapply(split(dgtex$refCount / dgtex$totalCount, dgtex$variantID), function(x) mean(x > 0.99 | x < 0.01))
# hist(het_prop_extreme, breaks = 0:100/100)
# quantile(het_prop_extreme)
# mean(het_prop_extreme > 0.75)
# 
# #hmm, just remove genes, not het sites?
# 
# #or filter out obs with eg > 0.99 het?
# hist(d$refCount / d$totalCount, breaks = 0:100/100)
# hist(0.5 - abs(d$refCount / d$totalCount - 0.5), breaks = 0:500/1000)
# min((0.5 - abs(d$refCount / d$totalCount - 0.5)))
# min(d$refCount)
# min(d$altCount)

#### Generate Quantities, OVC-vs-GTEx ####

model_path <- paste0(stan_model_dir, model_name, ".stan")
stan_code <- model_lines <- paste0(readLines(model_path), collapse = "\n")
model_block_lines <- parse_stan_blocks(stan_code = readLines(model_path))$model
gq_background <- model_block_lines[which(grepl("CONCENTRATION PARAMS", model_block_lines)) : which(grepl("COMBINATION", model_block_lines))]

#append generated quantities block
#want genewise effects of tumor grade vs GTEx
#also overall grade effects
#on conc, loc, logodds_AB, and cond_logodds_loh
#gtex=2 corresponds to the effect in OVC
gq_components <- c("

generated quantities{

", paste0(gq_background, collapse = "\n"),
"
  
  /////////////////////////
  //overall tumor effects//
  /////////////////////////
  
  vector[n_grade-1] conc_grade_csum_dev = conc_gtex_split[2] + conc_grade_csum[2:n_grade]; 
  vector[n_grade-1] conc_grade_csum_abs_mean = 100 + exp(log_conc_base + 4 + conc_grade_csum_dev);

  vector[n_grade-1] loc_grade_csum_dev = loc_gtex_split[2] + loc_grade_csum[2:n_grade]; 
  vector[n_grade-1] loc_grade_csum_abs_mean = 0.6 + inv_logit(logit_loc_base / 2 + 0 + loc_grade_csum_dev) / 10 * 3;

  vector[n_grade-1] logodds_AB_grade_csum_dev = logodds_AB_gtex_split[2] + logodds_AB_grade_csum[2:n_grade]; 
  vector[n_grade-1] logodds_AB_grade_csum_abs_mean = inv_logit(logodds_AB_base + 0 + logodds_AB_grade_csum_dev);

  vector[n_grade-1] cond_logodds_loh_grade_csum_dev = cond_logodds_loh_gtex_split[2] + cond_logodds_loh_grade_csum[2:n_grade]; 
  vector[n_grade-1] cond_logodds_loh_grade_csum_abs_mean = inv_logit(cond_logodds_loh_base + 0 + cond_logodds_loh_grade_csum_dev);
  
  /////////////////////
  //gene-wise effects//
  /////////////////////

  matrix[n_gene, n_grade-1] conc_grade_gene_csum_dev;
  matrix[n_gene, n_grade-1] conc_grade_gene_csum_dev_abs_mean;

  matrix[n_gene, n_grade-1] loc_grade_gene_csum_dev;
  matrix[n_gene, n_grade-1] loc_grade_gene_csum_dev_abs_mean;

  matrix[n_gene, n_grade-1] logodds_AB_grade_gene_csum_dev;
  matrix[n_gene, n_grade-1] logodds_AB_grade_gene_csum_dev_abs_mean;

  matrix[n_gene, n_grade-1] cond_logodds_loh_grade_gene_csum_dev;
  matrix[n_gene, n_grade-1] cond_logodds_loh_grade_gene_csum_dev_abs_mean;

  for(i in 1:n_gene){
    conc_grade_gene_csum_dev[i,] = to_row_vector(to_vector(conc_grade_gene_csum[i,2:n_grade]) + conc_grade_csum_dev);
    conc_grade_gene_csum_dev_abs_mean[i,] = to_row_vector(100 + exp(log_conc_base + 4 + conc_grade_gene_csum_dev[i,]));
  
    loc_grade_gene_csum_dev[i,] = to_row_vector(to_vector(loc_grade_gene_csum[i,2:n_grade]) + loc_grade_csum_dev);
    loc_grade_gene_csum_dev_abs_mean[i,] = to_row_vector(0.6 + inv_logit(logit_loc_base / 2 + 0 + loc_grade_gene_csum_dev[i,]) / 10 * 3);
  
    logodds_AB_grade_gene_csum_dev[i,] = to_row_vector(to_vector(logodds_AB_grade_gene_csum[i,2:n_grade]) + logodds_AB_grade_csum_dev);
    logodds_AB_grade_gene_csum_dev_abs_mean[i,] = to_row_vector(inv_logit(logodds_AB_base + 0 + logodds_AB_grade_gene_csum_dev[i,]));
  
    cond_logodds_loh_grade_gene_csum_dev[i,] = to_row_vector(to_vector(cond_logodds_loh_grade_gene_csum[i,2:n_grade]) + cond_logodds_loh_grade_csum_dev);
    cond_logodds_loh_grade_gene_csum_dev_abs_mean[i,] = to_row_vector(inv_logit(cond_logodds_loh_base + 0 + cond_logodds_loh_grade_gene_csum_dev[i,]));
  }
  
}
")
generated_quantities <- paste0(gq_components, collapse = "\n")

# Specify the new model with the generated quantities block
load(paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.fit"))
fit <- out

model_string_generated_quantities <- paste0(stan_code, generated_quantities, collapse = "\n")
# model_string_generated_quantities <- paste0(stan_code)
gq_model_path <- write_stan_file(model_string_generated_quantities)
gq_model <- cmdstan_model(gq_model_path)
# open_script(model_string_generated_quantities, "stan")
gq_var_names <- setdiff(retrieve_block_objects(stan_code = readLines(gq_model_path), block = "generated quantities")$var_name,
                        retrieve_block_objects(stan_code = readLines(gq_model_path), block = "model")$var_name)
gq_gene_params <- gq_var_names[grepl("gene", gq_var_names)]
gq_grade_params <- gq_var_names[!grepl("gene", gq_var_names)]

gene_map <- fread(file = paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_gene-map.csv"))

# Generate the quantities

#but first correct variable names?
# param_names <- fit$metadata()$model_params
# altformat_param_names <- gsub("\\[", ".", param_names)
# altformat_param_names <- gsub("]", "", altformat_param_names)
# altformat_param_names <- gsub(",", ".", altformat_param_names)
# csv_paths <- fit$output_files()
# for(i in seq_along(csv_paths)){
#   print(csv_paths[i])
#   csv_out <- readLines(csv_paths[i])
#   param_name_line <- which(grepl(altformat_param_names[12], csv_out))
#   csv_param_names <- csv_out[param_name_line]
#   csv_param_names_split <- str_split(csv_param_names, "(?<!\\d),(?!\\d)")[[1]]
#   
#   #forward direction
#   new_csv_param_names_split <- param_names[match(csv_param_names_split, 
#                                                  altformat_param_names)]
#   
#   #backward direction
#   new_csv_param_names_split <- altformat_param_names[match(csv_param_names_split, 
#                                                  param_names)]
#   
#   new_csv_param_names_split[is.na(new_csv_param_names_split)] <- csv_param_names_split[is.na(new_csv_param_names_split)]
#   new_csv_param_names <- paste(new_csv_param_names_split, collapse = ",")
#   csv_out[param_name_line] <- new_csv_param_names
#   writeLines(csv_out, con = csv_paths[i])
# }
#nope, did not work :/

gq_fit <- gq_model$generate_quantities(
  fitted_params = fit$output_files(), data = dat, parallel_chains = out$num_chains()
)

#that did not work :/
#also the Stan2R version is giving errors, probably because of the genquan block or bc of the to_row_vector functions or something

#### manually generate quantities ####
source("~/repos/Stan2R/R/functions.R")
model_path <- paste0(stan_model_dir, model_name, ".stan")
stan_code <- model_lines <- paste0(readLines(model_path), collapse = "\n")
samps <- data.table(as_draws_df(out$draws()))
param_names <- retrieve_block_objects(readLines(model_path), "parameters")$var_name
param_names <- setNm(param_names)
munged_samples_path <- paste0(stan_output_dir, model_name, "_", data_name, "_munged_samps.RData")
if(!file.exists(munged_samples_path)){
  munged_samps <- lapply(param_names, function(pname){
    print(pname)
    
    psamps <- subset_samps(pname, samps)
    msamps <- munge_samps(pname, psamps)
    
    if("array" %in% class(msamps[[1]])){
      outsamps <- abind(msamps, along = 3)
    } else {
      outsamps <- do.call(rbind, msamps)
    }
    
    return(outsamps)
  })
  save(munged_samps, file = munged_samples_path)
} else {
  load(munged_samples_path)
}

stan_code_gq <- readLines("Stan/models/ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc-gq.stan")
stan_code <- readLines("Stan/models/ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc.stan")
r_code <- parse_Stan(stan_code = stan_code_gq, 
           dat = dat, 
           samps = samps, 
           sample_index = 1, 
           post_pred_sim = T, 
           sim = F)
# open_script(r_code)

#### Stan2R version ####

#compute posterior predictive distributions using 
source("~/repos/Stan2R/R/functions.R")

#read from disk the original model?
model_path <- paste0(stan_model_dir, model_name, ".stan")
model_gq_path <- paste0(stan_model_dir, model_name, "-gq.stan")
stan_code <- model_lines <- readLines(model_path)
stan_code_gq <- readLines(model_gq_path)


#or use the generated quantities block?
# model_path <- gq_model_path
# stan_code <- model_lines <- model_string_generated_quantities

#get other useful variables
gq_variables <- setdiff(retrieve_block_objects(stan_code_gq, "model")$var_name,
                        retrieve_block_objects(stan_code, "model")$var_name)
gq_variables <- setNm(gq_variables)
nchains <- out$num_chains()
chains_to_compare <- 1:nchains
nsamps <- nrow(samps)
n_per_chain <- nsamps / nchains
sample_indices_list <- lapply(chains_to_compare, function(ci){
  sort(sample((nsamps / nchains * (ci-1) + 1):(nsamps / nchains * ci), size = n_per_chain, replace = F))
})
sample_indices <- unlist(sample_indices_list)


#process posterior predictive files
ppfits_path <- paste0(stan_output_dir, "posterior_predictive/", model_name, "_", 1, ".RData")
fsamps_path <- paste0(stan_output_dir, "posterior_predictive/", model_name, "_fsamps.RData")
if(file.exists(ppfits_path)){
  if(!file.exists(fsamps_path)){
    ppfits <- lapply(seq_along(sample_indices), function(i){
      if(i %% 500 == 0){mcprint(paste0(i, " "))}
      sample_index <- sample_indices[i]
      load(file = paste0(stan_output_dir, "posterior_predictive/", model_name, "_", sample_index, ".RData"))
      return(posterior_predictive)
    })  
  }
} else {
  ppfits <- parallel::mclapply(seq_along(sample_indices), function(i){
    
    if(i %% 50 == 0){mcprint(paste0(i, " "))}
    sample_index <- sample_indices[i]
    r_code <- parse_Stan(stan_code = stan_code_gq, 
                         dat = dat, 
                         samps = samps, 
                         output_file = NA, 
                         sample_index = sample_index, 
                         post_pred_sim = T, 
                         sim = F)
    # open_script(r_code)
    
    stan_env <- new.env()
    eval(parse(text = r_code), envir = stan_env)
    posterior_predictive <- stan_env$out
    posterior_predictive <- posterior_predictive[gq_variables]
    save(posterior_predictive, file = paste0(stan_output_dir, "posterior_predictive/", model_name, "_", sample_index, ".RData"))
    
    return(posterior_predictive)
    
  }, mc.cores = 12)
}

#process ppfits into appropriate data object
if(file.exists(fsamps_path)){
  load(fsamps_path)
} else {
  fsamps <- lapply(gq_variables, function(vname){
    print(vname)
    subsamps <- lapply(ppfits, function(ppfit){
      ppfit[[vname]]
    })
    var_dim <- max(1, length(dim(subsamps[[1]])))
    subsamps <- abind::abind(subsamps, along = var_dim + 1)
    if(var_dim == 1){
      subsamps <- t(subsamps)
    }
    return(subsamps)
  })
  save(fsamps, file = fsamps_path)
}

nv_string <- (paste0("nullvals <- c(\n", paste0("  ", gq_variables, " = ", 
                          c("", "0")[endsWith(suffix = "_dev", x = gq_variables) + 1], 
                          c("", "NA")[endsWith(suffix = "abs_mean", x = gq_variables) + 1], 
                          c("", "NA")[endsWith(suffix = "abs_gtex", x = gq_variables) + 1], 
                          c("", "0")[endsWith(suffix = "_abs_dev_gtex", x = gq_variables) + 1],
                          collapse = ",\n"), "\n)"))
cat(nv_string)
eval(parse(text = nv_string))

functions <- list(
  "mean" = function(x, ...) mean(x),
  "sd" = function(x, ...) sd(x),
  "CI95" = function(x, ...) quantile(x, probs = c(0.025, 0.975)),
  "q_null" = function(x, nullval = NA, ...) ifelse(is.na(nullval), NA, mean(nullval > x))
)

#metadata for gene IDs and grade IDs
grades <- setNames(object = c(1, 2, 2, 2, 3, 4, 5), 
         nm = c("GTEx", "Benign,\nUnstaged", "B", "X", "I", "II", "III"))
split(names(grades), grades)
genes <- fread(paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove, "_gene-map.csv"))

fsamps_summary <- lapply(gq_variables, function(vname){
  print(vname)
  subsamps <- fsamps[[vname]]
  
  var_dim <- length(dim(subsamps))
  if(var_dim == 2){
    posterior_summary <- lapply(functions, function(f) apply(subsamps, 2, function(x){
      f(x, nullval = nullvals[vname])
    }))
  } else if(var_dim == 3){
    if(grepl(pattern = "_gene", vname) && nrow(subsamps) == length(genes)){
      rownames(subsamps) <- genes$gene_symbol
    }
    posterior_summary <- lapply(functions, function(f) {
      ps <- apply(subsamps, c(1,2), function(x){
        f(x, nullval = nullvals[vname])
      })
      if(length(dim(ps)) == 3){
        return(aperm(ps, c(2, 3, 1)))
      } else {
        return(ps)
      }
    })
  }
})

# hist(fsamps_summary$cond_logodds_loh_grade_gene_csum_abs_dev_base$prob_null[,4], breaks = 100, xlim = c(0,1))
fsamps_summary$loc_grade_csum_dev
fsamps_summary$conc_grade_csum_dev
fsamps_summary$logodds_AB_grade_csum_dev
fsamps_summary$cond_logodds_loh_grade_csum_dev


fss <- fsamps_summary[gq_variables[endsWith(gq_variables, "abs_dev_gtex")]]
fss <- split(fss, c("grade", "gene")[grepl(pattern = "gene", x = names(fss)) + 1])

fss$grade$prob_slab_grade_csum_abs_dev_gtex$q_null
fss$grade$logodds_AB_grade_csum_abs_dev_gtex
sapply(names(fss$grade), function(vname){
  return(list(vname = fss$grade[[vname]]$q_null))
})

sapply(names(fss$grade), function(vname){
  return(list(vname = fss$grade[[vname]]$mean))
})

sapply(names(fss$gene), function(vname){
  return(list(vname = fss$gene[[vname]]$q_null[1,]))
})

hist(fss$gene$prob_slab_grade_gene_csum_abs_dev_gtex$q_null[,4])

#### plotting results ####
grDevices::cairo_pdf(filename = "figures/OVC_vs_GTEx.pdf", width = 400 / 72, height = 600 / 72)
source("~/scripts/minor_scripts/postdoc/geometric_word_labels_function.R") #source function for label plotting

mars <- list(hist = c(6.5,4,3,1), volc = c(4,3,1,1))

#specify ad hoc offsets for gene labels
offsets_x <- array(0, dim = c(4,2,length(vname_roots)), dimnames = list(NULL, c("pos", "neg"), vname_roots))
offsets_x[4, "neg", "logodds_AB"] <- 0.1
offsets_x[2, "pos", "logodds_AB"] <- -0.5
offsets_x[4, "pos", "loc"] <- -0.1

layout_mat_seed <- rbind(c(1,1,1,1), c(2,3,4,5))
layout_mat <- do.call(rbind, lapply(0:2, function(i) layout_mat_seed + i * 5))
layout(mat = layout_mat, heights = rep(c(1.5,1), nrow(layout_mat)/2))

fancy_names <- list(
  logodds_AB = "log-odds(allelic balance)",
  loc = "location of ASE",
  conc = "log(concentration) residual ASE",
  prob_slab = "Pr(ASE)",
  prob_loh = "Pr(LoH)",
  cond_logodds = "conditional log-odds(LoH)"
)

cols_grade <- viridis::rocket(4, begin = 0.3, end = 0.7)
rgb2 <- function(x) rgb(x[1], x[2], x[3])
dcols_grade <- sapply(seq_along(cols_grade), function(i) rgb2((col2rgb(cols_grade[i]) + 0) / 2 / 255))
col_gtex <- "darkslategray3"

vname_roots <- c("logodds_AB", "loc", "conc")
#start iterating through parameters
for(vname_root in vname_roots){

vname_grade <- names(fss$grade)[grepl(vname_root, names(fss$grade))]
vname_gene <- names(fss$gene)[grepl(vname_root, names(fss$gene))]

#retrieve inputs
vsamps_grade <- fsamps[[vname_grade]]
vsumm_grade <- fss$grade[[vname_grade]]
vsamps_gene <- fsamps[[vname_gene]]
vsumm_gene <- fss$gene[[vname_gene]]

#process grade histograms input
nbins <- 50
nbreaks <- nbins + 1
breaks <- seq(min(vsamps_grade), max(vsamps_grade), length.out = nbreaks)
hists_raw <- do.call(cbind, apply(vsamps_grade, 2, function(x) 
  hist(x, breaks = breaks, plot = F)$density, 
  simplify = F))
hists <- data.frame(lb = breaks[-(nbreaks)], 
                    ub = breaks[-1],
                    dens = hists_raw)
hist_polys <- lapply(1:4, function(i){
  data.frame(x = c(breaks[1:nbins], 
                   rep(rev(breaks), each = 2)
                   ),
             y = c(rep(0, nbreaks), 
                   rep(rev(hists_raw[1:nbins,i]), each = 2),
                   0)
             )
})

#process genic volcano plot input
alpha_sig <- 0.95
volc <- lapply(1:4, function(i){
  pmeans <- fss$gene[[vname_gene]]$mean[,i]
  psd <- fss$gene[[vname_gene]]$sd[,i]
  qs <- fss$gene[[vname_gene]]$q_null[,i]
  ppsign <- rep(0, length(pmeans))
  ppsign[sign(pmeans) == -1] <- qs[sign(pmeans) == -1]
  ppsign[sign(pmeans) == 1] <- (1-qs)[sign(pmeans) == 1]
  ppsign[ppsign == 1] <- 1 - 1 / (tail(dim(fsamps[[vname]]), 1))
  psig <- ppsign > alpha_sig
  return(data.frame(pmeans = pmeans, ppsign = ppsign, psd = psd, psig = psig, gene = genes$gene_symbol))
})
nlabs <- ifelse(vname_root == "logodds_AB",5, 3)
volc_genelabs <- lapply(1:4, function(i){
  volci <- volc[[i]]
  if(!any(volci$psig)){
    return(list(pos = NULL, neg = NULL))
  }
  volci <- volci[volci$psig,]
  v <- split(volci, c("neg", "pos")[sign(volci$pmeans) / 2 + 1.5])
  npos <- min(nlabs, ifelse(is.null(nrow(v$pos)), 0, nrow(v$pos)))
  nneg <- min(nlabs, ifelse(is.null(nrow(v$neg)), 0, nrow(v$neg)))
  pos <- NULL
  neg <- NULL
  if(npos > 0){
    pos <- v$pos[order(v$pos$pmeans, decreasing = T)[1:npos],]  
  }
  if(nneg > 0){
    neg <- v$neg[order(v$neg$pmeans, decreasing = F)[1:nneg],]
  }
  return(list(pos = pos, neg = neg))
})

#actual plotting

#plot histogram
par(mar = mars[["hist"]])

plot.new()
plot.window(xlim = range(c(0, breaks)), ylim = range(c(0, hists_raw)))
for(i in 1:4){
  polygon(hist_polys[[i]]$x, 
          hist_polys[[i]]$y, 
          col = adjustcolor(cols_grade[i], 0.5),
          border = dcols_grade[i])
  
  #label blobs
  if(vname_root != "conc"){
    labx <- hist_polys[[i]]$x[which.max(hist_polys[[i]]$y)] + ifelse(i == 3, 0.02, 0)
    laby <- max(hist_polys[[i]]$y) * 1.15
    text(x = labx, y = laby,
         labels = ifelse(i == 1, "Benign / Unstaged", 
                         paste0("Grade ", paste0(rep("I", i-1), collapse = ""))),
         col = dcols_grade[i], xpd = NA)  
  }
  
}

#histogram axes
axis(1, at = pretty(range(c(0, breaks))), labels = pretty(range(c(0, breaks))), line = 0.1)
axis(2, line = 0.1)

mtext(side = 1, 
      text = latex2exp::TeX(paste0(
        "Expected Deviation in \\textbf{", 
        fancy_names[vname_root], 
        "} from GTEx by Grade")
      ), 
      line = 2.5, cex = 0.7)
mtext(side = 2, 
      text = "Density", 
      line = 2.5, cex = 0.7)
abline(v = 0,lty = 2, col = adjustcolor(1,0.5))


#plot gene-specific estimates
par(mar = mars[["volc"]], mgp = c(3, 0.5, 0))

#can maybe also look at pooled mu for the average linear gene-wise effect across grades?
cexpts <- 1.5
incr_order <- sign(mean(sign(diff(vsumm_grade$mean))) + rnorm(1)/1E3) == 1
ifelse2 <- function(test, yes, no) if(test){return(yes)}else{return(no)}
order_genes <- ifelse2(incr_order, 1:4, 4:1)
pmeans_range <- range(unlist(lapply(order_genes, function(i) volc[[i]]$pmeans)))
for(i in order_genes){
  
  plot.new()
  plot.window(xlim = pmeans_range, ylim = range(c(0.5, 1)))
  
  #plot points
  points(x = volc[[i]]$pmeans[volc[[i]]$psig], 
         y = volc[[i]]$ppsign[volc[[i]]$psig], 
         col = adjustcolor(cols_grade[i], 0.5),
         cex = cexpts,
         pch = 19, xpd = NA)
  points(x = volc[[i]]$pmeans[!volc[[i]]$psig], 
         y = volc[[i]]$ppsign[!volc[[i]]$psig], 
         col = adjustcolor(cols_grade[i], 0.05),
         cex = cexpts,
         pch = 19, xpd = NA)
  points(x = volc[[i]]$pmeans[volc[[i]]$psig], 
         y = volc[[i]]$ppsign[volc[[i]]$psig], 
         col = cols_grade[i],
         cex = cexpts, xpd = NA)
  abline(h = alpha_sig ,lty = 2, col = adjustcolor(1,0.5))
  segments(x0 = 0, x1 = 0, y0 = -10, y1 = alpha_sig, lty = 2, col = adjustcolor(1,0.5))
  axis(1, line = 0.2, cex.axis = 0.75)

  if(i==4){
    axis(2, line = 0.2, labels = 5:10/10, at = 5:10/10, cex.axis = 0.75)
    mtext(side = 2, text = "Pr(Sign of Effect)", 
          line = 2, cex = 0.5)
  }
  
  mtext(side = 1, text = "Expected Genic Effect", 
        line = 2.25, cex = 0.5)
  mtext(side = 1, 
        text = latex2exp::TeX(paste0("for \\textbf{", 
          ifelse(i == 1, 
                 "Benign / Unstaged", 
                 paste0("Grade ", paste0(rep("I", i-1), collapse = ""))), 
          "}")
        ), 
        line = ifelse(i == 1, 3.25, 3), cex = 0.5)
  
  
  #label top points
  for(j in c("pos", "neg")){
    vgl <- volc_genelabs[[i]][[j]]
    if(is.null(vgl)){
      next()
    } else {
      par(xpd = NA)
      offset_x <- offsets_x[i, j, vname_root] 
      stacked_word_labels(words = vgl$gene, wcols = rep(dcols_grade[i], nrow(vgl)), 
                          x0 = vgl$pmeans, y0 = vgl$ppsign, 
                          nro = 3, offx = 0.01, variable_priority = T, offset_y = -0.05, 
                          wsh_scale = 0.5, offset_x = offset_x)
      par(xpd = F)
    }
  }
}

}

dev.off()

#### look at gross grade effects ####

#need location effect, concentration effect, logodds_AB, and cond_logodds_loh
#can use *_grade_csum + -*_gtex / 2 (if GTEx = 2, you sample is cancer)
#or *_grade_csum + *_gtex_split[2]
model_vars <- rbind(retrieve_block_objects(stan_code, "model"),
                        retrieve_block_objects(stan_code, "parameters"))
focal_prefixes <- c("conc", "loc", "logodds_AB", "cond_logodds_loh")
focal_suffixes <- c("grade_csum", "gtex_split")
focal_param_subset <- apply(expand.grid(focal_prefixes, focal_suffixes), 1, paste0, collapse = "_")
focal_param_subset <- c(focal_param_subset, model_vars$var_name[grepl("_base", model_vars$var_name)])

ppfit_sub <- lapply(seq_along(ppfits), function(i){
  ppfits[[i]][focal_param_subset]
})

focal_samps <- lapply(setNm(focal_param_subset), function(pname){
  do.call(rbind, lapply(ppfit_sub, function(x) x[[pname]]))
})

grade_effects <- lapply(setNm(focal_prefixes), function(pi){
  t(t(focal_samps[[paste0(pi, "_grade_csum")]][,-1]) + focal_samps[[paste0(pi, "_gtex_split")]][,2])
})

xbreaks_list <- lapply(setNm(names(grade_effects)), function(pi){
  lims <- quantile(grade_effects[[pi]], c(0.001, 0.999))
  breaks <- seq(lims[1], lims[2], length.out = 100)
  breaks
})

interpretations <- setNames(c(
  "\\textit{(negative values mean there's more residual variation / overdispersion in the binomial count)}",
  "\\textit{(positive values mean you're further away from allelic balance)}",
  "\\textit{(negative values mean a lower probability (log-odds) of allelic balance)}",
  "\\textit{(negative values mean a lower probability (log-odds) of LoH and a higher probability of ASE)}"
), focal_prefixes)
fullnames <- c(conc = "ASE Concentration", loc = "ASE Location", 
               logodds_AB = "log-odds(Allelic Balance)",
               cond_logodds_loh = "Conditional log-odds(LoH)")
grades <- setNames(object = c(1, 2, 2, 2, 3, 4, 5), 
                   nm = c("GTEx", "Benign,\nUnstaged", "B", "X", "I", "II", "III"))

var_name <- c("conc", "loc", "logodds_AB", "cond_logodds_loh")[2]
par(mfrow = c(5,1), mar = c(3,4,3,2))
for(i in 1:4){
  
  x <- grade_effects[[var_name]][,i]
  main <- paste0("Mean ", fullnames[var_name]," Effect for \\textbf{Grade ", 
                 paste0(gsub("\\n", " ", names(grades)[grades == i + 1]), collapse = " + "), "} vs. GTEx")
  breaks <- xbreaks_list[[var_name]]
  prxg0 <- mean(x > 0)
  prxg0 <- ifelse(prxg0 > 0.5, paste0("\\textit{Pr(X>0) = ", prxg0, "}"), paste0("\\textit{Pr(X<0) = ", 1 - prxg0, "}"))
  x_sub <- x[x > head(breaks, 1) & x < tail(breaks, 1)]
  hist(x_sub, xlim = range(breaks), breaks = breaks, xlab = paste0(var_name, "[", i, "]"), main = latex2exp::TeX(main))
  legend("topleft", legend = latex2exp::TeX(prxg0), bty = "n")
  abline(v = 0, lty = 2, col = adjustcolor(2, 1), lwd = 2)
  
  if(i == 4){
    text(latex2exp::TeX(interpretations[[var_name]]),
         x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4]), xpd = NA, pos = 1)
    plot.new()
  }
}

#### look at gene-specific grade effects ####
model_vars <- rbind(retrieve_block_objects(stan_code, "model"),
                    retrieve_block_objects(stan_code, "parameters"))
focal_prefixes <- c("conc", "loc", "logodds_AB", "cond_logodds_loh")
focal_suffixes <- c("grade_csum", "gtex_split", "grade_gene_csum", "gene_gtex_diff_split")
focal_param_subset <- apply(expand.grid(focal_prefixes, focal_suffixes), 1,
                            paste0, collapse = "_")

focal_param_subset <- c(focal_param_subset, 
                        paste0("sd_", focal_prefixes, "_gtex"),
                        model_vars$var_name[grepl("_base", model_vars$var_name)])

ppfit_sub <- lapply(seq_along(ppfits), function(i){
  ppfits[[i]][focal_param_subset]
})

focal_samps <- lapply(setNm(focal_param_subset), function(pname){
  
  param_list <- lapply(ppfit_sub, function(x) x[[pname]])
  if (is.matrix(param_list[[1]])) {
    do.call(abind, c(param_list, along = 3))
  } else {
    do.call(rbind, param_list)
  }

})

grade_x_gene_effects <- lapply(setNm(focal_prefixes), function(pi){
  grade_means <- t(
    t(focal_samps[[paste0(pi, "_grade_csum")]][,-1]) + 
      focal_samps[[paste0(pi, "_gtex_split")]][,2]
  )
  
  grade_x_gene <- focal_samps[[paste0(pi, "_grade_gene_csum")]][,-1,]
  gtex_x_gene <- sweep(t(focal_samps[[paste0(pi, "_gene_gtex_diff_split")]][,2,]), 
        1, 
        focal_samps[[paste0("sd_", pi, "_gtex")]] * 0.5, 
        FUN = "*")
  
  gene_devs <- grade_x_gene + aperm(abind(lapply(1:4, function(x) gtex_x_gene), along = 3), c(2,3,1))
  
  gene_effects <- aperm(abind(lapply(1:253, function(x) grade_means), along = 3), c(3,2,1)) + gene_devs
  dim(gene_effects)
  return(apply(gene_effects, c(1, 2), function(x) mean(x < 0)))
})

var_name <- c("conc", "loc", "logodds_AB", "cond_logodds_loh")[3]
par(mfrow = c(4,1), mar = c(3,4,3,2))
for(i in 1:4){
  
  
  x <- grade_x_gene_effects[[var_name]][,i]
  main <- paste0("Prob ", fullnames[var_name]," Gene Effect for \\textbf{Grade ", 
                 paste0(gsub("\\n", " ", names(grades)[grades == i + 1]), collapse = " + "), "} is negative vs. GTEx")
  breaks <- 0:20/20
  x_sub <- x[x > head(breaks, 1) & x < tail(breaks, 1)]
  hist(x_sub, xlim = range(breaks), breaks = breaks, xlab = paste0(var_name, "[", i, "]"), main = latex2exp::TeX(main))
  
}
  
sum(grade_x_gene_effects$logodds_AB[,2] > 0.95)






setNm <- function(x) setNames(x,x)
focal_samps <- lapply(setNm(focal_param_subset), function(pname){
  do.call(rbind, lapply(ppfit_sub, function(x) x[[pname]]))
})



n_iter_per_chain <- sapply(sample_indices_list, length)
iter_per_chain <- lapply(seq_along(n_iter_per_chain), function(i) 1:n_iter_per_chain[i] + ifelse(i != 1, sum(n_iter_per_chain[1:(i-1)]), 0))
sim_varnames <- setNames(names(ppfits[[1]]), names(ppfits[[1]]))
ppfits_per_chain_means <- lapply(iter_per_chain, function(eyes){
  sum_of_params <- lapply(sim_varnames, function(varname){
    Reduce("+", lapply(eyes, function(i) {
      ppfits[[i]][[varname]]
    })) / length(eyes)
  })
})


#simulation code for posterior predictive simulation
#(and also for generating new data to fit)
# subset_samps <- function(include = "", exclude = "", samps){
#   incl_inds <- unique(unlist(lapply(include, function(i) grep(i, colnames(samps)))))
#   if(exclude != ""){
#     excl_inds <- unique(unlist(lapply(exclude, function(i) grep(i, colnames(samps)))))
#     return_inds <- setdiff(incl_inds, excl_inds)  
#   } else {
#     return_inds <- incl_inds
#   }
#   return(samps[,return_inds])
# }


#fake data
# d_sim <- d
# subsamp <- samps[sample(nrow(samps), 1),]
# conc_gtex_split <- c(subsamp$conc_gtex / 2, -subsamp$conc_gtex/2)
# conc_gene_gtex_diff <- subset_samps(var_name = "conc_gene_gtex_diff", samps = subsamp)


#basic plotting of marginal histograms

# overall params #

#location sds
par(mfrow = c(3,1), mar = c(5,6,4,2))
range_x <- range(c(samps$sd_loc_gene, samps$sd_loc_gtex, samps$sd_loc_indiv))
xbreaks <- seq(range_x[1], range_x[2], length.out = 100)
hist(samps$sd_loc_gene, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\mu_i$"),
     main = latex2exp::TeX("Posterior Distribution of Gene SD on Location ($\\sigma^\\mu_i$)"),
     cex.lab = 1.1, cex.main = 1.4)
hist(samps$sd_loc_gtex, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\mu_{i,k}$"),
     main = latex2exp::TeX("Posterior Distribution of Gene x GTEx SD on Location ($\\sigma^\\mu_{i,k}$)"), 
     cex.lab = 1.1, cex.main = 1.4)
hist(samps$sd_loc_indiv, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\mu_i$"),
     main = latex2exp::TeX("Posterior Distribution of Individual SD on Location ($\\sigma^\\mu_j$)"),
     cex.lab = 1.1, cex.main = 1.4)

#concentration sds
par(mfrow = c(3,1), mar = c(5,6,4,2))
range_x <- range(c(samps$sd_conc_gene, samps$sd_conc_gtex, samps$sd_conc_indiv))
xbreaks <- seq(range_x[1], range_x[2], length.out = 100)
hist(samps$sd_conc_gene, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\nu_i$"),
     main = latex2exp::TeX("Posterior Distribution of Gene SD on Concentration ($\\sigma^\\nu_i$)"),
     cex.lab = 1.1, cex.main = 1.4)
hist(samps$sd_conc_gtex, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\nu_{i,k}$"),
     main = latex2exp::TeX("Posterior Distribution of Gene x GTEx SD on Concentration ($\\sigma^\\nu_{i,k}$)"), 
     cex.lab = 1.1, cex.main = 1.4)
hist(samps$sd_conc_indiv, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\nu_i$"),
     main = latex2exp::TeX("Posterior Distribution of Individual SD on Concentration ($\\sigma^\\nu_j$)"),
     cex.lab = 1.1, cex.main = 1.4)

#pointmass mixture probabilities sds
par(mfrow = c(3,1), mar = c(5,6,4,2))
range_x <- range(c(samps$sd_ppm_gene, samps$sd_ppm_gtex, samps$sd_ppm_indiv))
xbreaks <- seq(range_x[1], range_x[2], length.out = 100)
hist(samps$sd_ppm_gene, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\pi_i$"),
     main = latex2exp::TeX("Posterior Distribution of Gene SD on Mixture Probability ($\\sigma^\\pi_i$)"),
     cex.lab = 1.1, cex.main = 1.4)
hist(samps$sd_ppm_gtex, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\pi_{i,k}$"),
     main = latex2exp::TeX("Posterior Distribution of Gene x GTEx SD on Mixture Probability ($\\sigma^\\pi_{i,k}$)"), 
     cex.lab = 1.1, cex.main = 1.4)
hist(samps$sd_ppm_indiv, breaks = xbreaks, xlab = latex2exp::TeX("$\\sigma^\\pi_i$"),
     main = latex2exp::TeX("Posterior Distribution of Individual SD on Mixture Probability ($\\sigma^\\pi_j$)"),
     cex.lab = 1.1, cex.main = 1.4)

hist(samps$sd_ppm_gene, breaks = 100)
hist(samps$sd_ppm_gtex, breaks = 100)
hist(samps$sd_ppm_indiv, breaks = 100)

#average values
hist(invlogit(samps$logit_loc_base) / 2, breaks = 100)
hist(exp(samps$log_conc_base), breaks = 100)
hist(invlogit(samps$logit_ppm_base), breaks = 100)
#hmm, picking up on tiny deviations away from exactly 0.05 eh?
#may want to make the sparsity pointmass into a small beta

#average gtex effects
hist(invlogit(samps$loc_gtex) / 2 - 0.25, breaks = 100, 
     main = paste0("Average GTEx Effect on Location (Pr(X>0) ≈ ", 
                   mean(invlogit(samps$loc_gtex) / 2 - 0.25 > 0), ")"), 
     xlab = "Displacement Away from 0.5 (@ ±0.25 equiv)"); abline(v = 0, col = 2, lwd = 2)

hist(exp(samps$conc_gtex), breaks = 100, 
     main = paste0("Average GTEx Effect on Concentration (Pr(X>1) ≈ ", 
                   mean(samps$conc_gtex > 0), ")"), xlab = "Factor")
abline(v = 1, col = 2, lwd = 2)

# hist(samps$ppm_gtex, breaks = 100,  
#      main = paste0("Average GTEx Effect on log-odds of\ncoming from spike (Pr(X>0) ≈ ", 
#                    mean(samps$ppm_gtex > 0), ")"), xlab = "log-odds"); 
# abline(v = 0, col = 2, lwd = 2)

hist(invlogit(samps$ppm_gtex) - 0.5, breaks = 100,  
     main = paste0("Average GTEx Effect on log-odds of\ncoming from spike (Pr(X>0) ≈ ", 
                   mean(samps$ppm_gtex > 0), ")"), xlab = "log-odds (@ deviation from 0.5 equiv)"); 
abline(v = 0, col = 2, lwd = 2)
#this gives difference in probability between gtex coming from binomial(p = 0.5) and OVC
#from a baseline of 0.5
#ppm_gtex is the effect of ppm(gtex = 1) - ppm(gtex = 2)
#gtex = 2 is ovarian cancer, gtex = 1 is gtex, so when ppm_gtex is positive
#that means that ppm(gtex = 1) is larger than ppm(gtex = 3), or equiv that
#ppm for gtex smamples is larger (more positive) than for OVC samples

#gtex-gene interaction term variation
hist(samps$sd_conc_gtex, breaks = 100)
hist(samps$sd_loc_gtex, breaks = 100)
hist(samps$sd_ppm_gtex, breaks = 100)

#gtex-x-gene differences in ASE
ppi_ppm_samps <- subset_samps(var_name = "ppm_gene_gtex_diff", samps = samps)
ppi_loc_samps <- subset_samps("loc_gene_gtex_diff", samps = samps)
ppi_conc_samps <- subset_samps("conc_gene_gtex_diff", samps = samps)

prob_greater_0 <- function(x) mean(x>0)

par(mfrow = c(3,1), mar = c(5,6,4,2))
ppi_ppm_logodds_samps <- ppi_ppm_samps * samps$sd_ppm_gtex + samps$ppm_gtex
ppm_gtex_probs <- apply(ppi_ppm_logodds_samps, 2, prob_greater_0)
hist(ppm_gtex_probs, breaks = 0:100/100, xlab = "Posterior Probability",
     main = latex2exp::TeX("Probability of Positive Gene-wise effect of GTEx on $\\pi$"))
mean(ppm_gtex_probs > 0.95)
genes[ppm_gtex_probs > 0.95]
paste0(genes[ppm_gtex_probs > 0.95], collapse = ", ")
genes[which.max(ppm_gtex_probs)]
paste0(sample(genes, 10), collapse = ", ")
max(ppm_gtex_probs)

ppi_loc_logodds_samps <- ppi_loc_samps * samps$sd_loc_gtex + samps$loc_gtex
loc_gtex_probs <- apply(ppi_loc_logodds_samps, 2, prob_greater_0)
hist(loc_gtex_probs, breaks = 0:100/100, xlab = "Posterior Probability",
     main = latex2exp::TeX("Probability of Positive Gene-wise effect of GTEx on $\\mu$"))
mean(loc_gtex_probs > 0.95)
max(loc_gtex_probs)

ppi_conc_logodds_samps <- ppi_conc_samps * samps$sd_conc_gtex + samps$conc_gtex
conc_gtex_probs <- apply(ppi_conc_logodds_samps, 2, prob_greater_0)
hist(conc_gtex_probs, breaks = 0:100/100, xlab = "Posterior Probability",
     main = latex2exp::TeX("Probability of Positive Gene-wise effect of GTEx on $\\nu$"))
mean(conc_gtex_probs > 0.95)
max(conc_gtex_probs)

#if we have grade information, let's also look at grade

#gene-wise sd in grade
par(mfrow = c(4,1), mar = c(3,4,3,2))
xlims <- quantile(unlist(subset_samps(var_name = 
                                        c("sd_conc_grade_gene", "sd_loc_grade_gene", "sd_ppm_grade_gene"), samps = samps)),
                  prob = c(1E-3, 1-1E-3))
for(i in 1:4){hist(subset_samps(var_name = "sd_conc_grade_gene", samps = samps)[,i], xlim = xlims)}
for(i in 1:4){hist(subset_samps(var_name = "sd_loc_grade_gene", samps = samps)[,i], xlim = xlims)}
for(i in 1:4){hist(subset_samps(var_name = "sd_ppm_grade_gene", samps = samps)[,i], xlim = xlims)}

#bulk grade effect
all_vals <- unlist(subset_samps(var_name = c("loc_grade", "conc_grade", "ppm_grade"), samps = samps)) * 
  rep(samps$sd_ppm_grade, 3) - 
  rep(samps$ppm_gtex, 3) / 2
raw_xlims <- quantile(all_vals, prob = c(1E-3, 1-1E-3))
nbreaks <- 20
hist_interval <- 10^floor(log10(diff(raw_xlims) / nbreaks))
hist_interval <- hist_interval * round(diff(raw_xlims) / hist_interval / nbreaks)
xlims <- c(floor(raw_xlims[1] / hist_interval), 
           ceiling(raw_xlims[2] / hist_interval)) * hist_interval
breaks <- seq(xlims[1], xlims[2], by = hist_interval)
breaks <- c(rev(seq(xlims[1] - hist_interval, min(all_vals) - hist_interval, by = -hist_interval)),
            breaks,
            seq(xlims[2] + hist_interval, max(all_vals) + hist_interval, by = hist_interval))

par(mfrow = c(5,1), mar = c(3,4,3,2))
for(i in 1:4){
  var_name <- "loc_grade"
  x <- subset_samps(var_name = var_name, samps = samps)[,i] * samps$sd_loc_grade - samps$loc_gtex / 2
  main <- paste0("Mean Location Effect for \\textbf{Grade ", 
                 paste0(gsub("\\n", " ", names(grades)[grades == i + 1]), collapse = " + "), "} vs. GTEx")
  prxg0 <- mean(x > 0)
  prxg0 <- ifelse(prxg0 > 0.5, paste0("\\textit{Pr(X>0) = ", prxg0, "}"), paste0("\\textit{Pr(X<0) = ", 1 - prxg0, "}"))
  hist(x, xlim = xlims, breaks = breaks, xlab = paste0(var_name, "[", i, "]"), main = latex2exp::TeX(main))
  legend("topleft", legend = latex2exp::TeX(prxg0), bty = "n")
  abline(v = 0, lty = 2, col = adjustcolor(2, 1), lwd = 2)
  
  if(i == 4){
    text(latex2exp::TeX("\\textit{(positive values mean you're further away from allelic balance)}"),
         x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4]), xpd = NA, pos = 1)
    plot.new()
  }
}

for(i in 1:4){
  var_name <- "conc_grade"
  x <- subset_samps(var_name = var_name, samps = samps)[,i] * samps$sd_conc_grade - samps$conc_gtex / 2
  main <- paste0("Mean Concentration Effect for \\textbf{Grade ", 
                 paste0(gsub("\\n", " ", names(grades)[grades == i + 1]), collapse = " + "), "} vs. GTEx")
  prxg0 <- mean(x > 0)
  prxg0 <- ifelse(prxg0 > 0.5, paste0("\\textit{Pr(X>0) = ", prxg0, "}"), paste0("\\textit{Pr(X<0) = ", 1 - prxg0, "}"))
  hist(x, xlim = xlims, breaks = breaks, xlab = paste0(var_name, "[", i, "]"), main = latex2exp::TeX(main))
  legend("topleft", legend = latex2exp::TeX(prxg0), bty = "n")
  abline(v = 0, lty = 2, col = adjustcolor(2, 1), lwd = 2)
  
  if(i == 4){
    text(latex2exp::TeX("\\textit{(negative values mean there's more residual variation / overdispersion in the binomial count)}"),
         x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4]), xpd = NA, pos = 1)
    plot.new()
  }
}

for(i in 1:4){
  var_name <- "ppm_grade"
  x <- subset_samps(var_name = var_name, samps = samps)[,i] * samps$sd_ppm_grade - samps$ppm_gtex / 2
  main <- paste0("Mean Mixture Probability Effect for \\textbf{Grade ", 
                 paste0(gsub("\\n", " ", names(grades)[grades == i + 1]), collapse = " + "), "} vs. GTEx")
  prxg0 <- mean(x > 0)
  prxg0 <- ifelse(prxg0 > 0.5, paste0("\\textit{Pr(X>0) = ", prxg0, "}"), paste0("\\textit{Pr(X<0) = ", 1 - prxg0, "}"))
  hist(x, xlim = xlims, breaks = breaks, xlab = paste0(var_name, "[", i, "]"), main = latex2exp::TeX(main))
  legend("topleft", legend = latex2exp::TeX(prxg0), bty = "n")
  abline(v = 0, lty = 2, col = adjustcolor(2, 1), lwd = 2)
  
  if(i == 4){
    text(latex2exp::TeX("\\textit{(negative values mean a lower probability (log-odds) of allelic balance)}"),
         x = mean(par("usr")[1:2]), y = par("usr")[3] - diff(par("usr")[3:4]), xpd = NA, pos = 1)
    plot.new()
  }
}


ppi_ppm_grade_samps <- subset_samps(var_name = "ppm_grade_gene", samps = samps)
ppi_loc_grade_samps <- subset_samps(var_name = "loc_grade_gene", samps = samps)
ppi_conc_grade_samps <- subset_samps(var_name = "conc_grade_gene", samps = samps)

ppi_ppm_grade_samps <- subset_samps(var_name = "sd_ppm_grade_gene", samps = samps)
ppi_loc_grade_samps <- subset_samps(var_name = "sd_loc_grade_gene", samps = samps)
ppi_conc_grade_samps <- subset_samps(var_name = "sd_conc_grade_gene", samps = samps)

ppi_ppm_grade_samps <- subset_samps(var_name = "ppm_grade_gene", samps = samps)
ppi_loc_grade_samps <- subset_samps(var_name = "loc_grade_gene", samps = samps)
ppi_conc_grade_samps <- subset_samps(var_name = "conc_grade_gene", samps = samps)

#gtex-x-gene differences in ASE
ppi_ppm_samps <- subset_samps(var_name = "ppm_gene_gtex_diff", samps = samps)
ppi_loc_samps <- subset_samps("loc_gene_gtex_diff", samps = samps)
ppi_conc_samps <- subset_samps("conc_gene_gtex_diff", samps = samps)

par(mfrow = c(3,1), mar = c(5,6,4,2))
ppi_ppm_logodds_samps <- ppi_ppm_samps * samps$sd_ppm_gtex + samps$ppm_gtex
ppm_gtex_probs <- apply(ppi_ppm_logodds_samps, 2, prob_greater_0)
hist(ppm_gtex_probs, breaks = 0:100/100, xlab = "Posterior Probability",
     main = latex2exp::TeX("Probability of Positive Gene-wise effect of GTEx on $\\pi$"))
mean(ppm_gtex_probs > 0.95)
genes[ppm_gtex_probs > 0.95]
paste0(genes[ppm_gtex_probs > 0.95], collapse = ", ")
genes[which.max(ppm_gtex_probs)]
paste0(sample(genes, 10), collapse = ", ")
max(ppm_gtex_probs)

ppi_loc_logodds_samps <- ppi_loc_samps * samps$sd_loc_gtex + samps$loc_gtex
loc_gtex_probs <- apply(ppi_loc_logodds_samps, 2, prob_greater_0)
hist(loc_gtex_probs, breaks = 0:100/100, xlab = "Posterior Probability",
     main = latex2exp::TeX("Probability of Positive Gene-wise effect of GTEx on $\\mu$"))
mean(loc_gtex_probs > 0.95)
max(loc_gtex_probs)

ppi_conc_logodds_samps <- ppi_conc_samps * samps$sd_conc_gtex + samps$conc_gtex
conc_gtex_probs <- apply(ppi_conc_logodds_samps, 2, prob_greater_0)
hist(conc_gtex_probs, breaks = 0:100/100, xlab = "Posterior Probability",
     main = latex2exp::TeX("Probability of Positive Gene-wise effect of GTEx on $\\nu$"))
mean(conc_gtex_probs > 0.95)
max(conc_gtex_probs)

#to calculate gene-wise differences between gtex and OVC samps, need to
#recreate the parameters, ie add together:
invlogit <- function(x){exp(x) / (1 + exp(x))}
invlogit(logit_ppm_base + 
           ppm_gene[gene_i] * sd_ppm_gene + 
           ppm_gtex_split[gtex] + 
           ppm_gene_gtex_diff_split_vec * sd_ppm_gtex)

#and multiply it by 0,

#and then take that from 1 and multiply it by draws from the beta,
#doing the same with the mean and conc parameters




# hist(samps$loc  )

# 
# #eda
# hist(table(ocvTestCount$variantID))
# hist(table(ocvTestCount$gene), breaks = 100)
# 
# fit = vglm(cbind(altCount, totalCount - altCount) ~ log(totalCount),
#            family = betabinomial(lmu = "identitylink", nsimEIM=100),
#            data = ocvTestCount)
# head(ocvTestCount$totalCount)
# head(ocvTestCount$altCount)
# nrow(ocvTestCount)
# 
# mod <- cmdstan_model("Stan/beta-binomial.stan")
# cat(paste0(readLines("Stan/beta-binomial.stan"), collapse = "\n"))
# d <- list(n = nrow(ocvTestCount),
#           x = ocvTestCount$altCount,
#           k = ocvTestCount$totalCount)
# plot(d$x / d$k, d$k, pch = 19, col = adjustcolor(1, 0.05))
# out <- mod$sample(chains = 4, iter_sampling = 5E2, iter_warmup = 5E2, 
#                   data = d, parallel_chains = 4, adapt_delta = 0.95)
# # summ <- out$summary()
# # summ[order(summ$ess_bulk),]
# # summ[order(summ$rhat, decreasing = T),]
# 
# samps <- data.frame(as_draws_df(out$draws()))
# # mu_samps <- subset_samps(samps = samps, include = "mu")
# # hist(apply(mu_samps, 2, mean))
# samps$shape1.1.


# now let's examine the effects of grade
hist(samps$mu_ppm_grade)
mean(samps$mu_ppm_grade > 0)
hist(samps$mu_loc_grade)
mean(samps$mu_loc_grade > 0)
hist(samps$mu_conc_grade)
mean(samps$mu_conc_grade > 0)

hist(samps$ppm_grade.1. )
hist(samps$ppm_grade.2.)
hist(samps$ppm_grade.3.)  