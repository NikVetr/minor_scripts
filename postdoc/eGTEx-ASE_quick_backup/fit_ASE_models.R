#specifications
setwd("~/repos/egtex-ase")
sample_prop_tail_to_remove <- 0

#### Load Packages ####
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

#### create useful variables ####
stan_model_dir <- "Stan/models/"
stan_data_dir <- "Stan/data/"
stan_output_dir <- "Stan/output/genes/"
stan_progress_dir <- "Stan/progress/genes/"
tissue_codes_df <- fread("tissue_codes.txt")
tissue_codes <- setNames(tissue_codes_df$code, tissue_codes_df$name)

#### fit for GTEx ####
use_rnaseq <- F
base_fit <- F
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
stan_data_genes_dir <- paste0(stan_data_dir, "genes/")


model_index <- 6
model_names <- list(
  "ase_param-expansion_1gene", #1
  "ase_param-expansion_by-gene_eQTL", #2
  "ase_param-expansion_by-gene_raw-params_eQTL", #3
  "ase_param-expansion_by-gene_raw-params_resid-indivs_eQTL", #4
  "ase_param-expansion_by-gene_raw-params_resid-indivs_infer-bounds_eQTL", #5
  "ase_param-expansion_by-gene_raw-params_resid-indivs_eQTL_empirical" #6
)

model_name <- model_names[[model_index]]
sample_prop_tail_to_remove <- 0
model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- cmdstan_model(model_path)
mod_no_eQTLHets <- cmdstan_model(gsub("eQTL", "no-eQTL", model_path))
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
cat(model_path)

#set MCMC inputs
sampling_param_set <- 2
nchains <- 4
niter <- c(5E3, 1E4, 2E4, 2E3)[sampling_param_set]
adapt_delta <- c(0.9, 0.95, 0.95, 0.8)[sampling_param_set]
max_treedepth <- c(10, 12, 15, 10)[sampling_param_set]
thin <- c(2, 4, 8, 1)[sampling_param_set]
init_mcmc <- 1

#get some useful variables
model_blocks <- parse_stan_blocks(clean_stan(model_lines))
model_blocks <- names(model_blocks)[sapply(model_blocks, length) != 0]
setNm <- function(x) setNames(x, x)
model_objects <- lapply(setNm(model_blocks), retrieve_block_objects, stan_code = model_lines)

#start analysis
tissue_names <- setNames(names(tissue_codes), tissue_codes)
simultaneous_sets <- 2
tissue_codes_to_use <- tissue_codes
tissue_sets <- suppressWarnings(split(tissue_codes_to_use, f = 1:simultaneous_sets))

for(tissue_set in tissue_sets){
  
  future({
    
    for(tiss_j in tissue_set){
      
      #specify names for files
      base_file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                          ifelse(use_rnaseq, "_rnaseq", ""))
      base_data_name <- paste0(tiss_j, "_", sample_prop_tail_to_remove)
      
      #read in genes
      gene_map <- fread(file = paste0(stan_data_dir, tiss_j, "_", 
                                      sample_prop_tail_to_remove, "_gene-map.csv"))
      genes <- gene_map$gene_symbol
      
      for(gene_i in genes[1:10]){
        
        #read in the data
        data_name <- paste0(stan_data_genes_dir, base_data_name, "_", gene_i, ".json")
        dat <- jsonlite::fromJSON(data_name)
        
        #check if hets are present and select model accordingly
        if(all(dat$eQTL_het == 1) || all(dat$eQTL_het == 2)){
          model <- mod
        } else {
          model <- mod_no_eQTLHets
        }
        
        #specify output name structure
        file_name <- paste0(ifelse(base_fit, "1st-fit_", ""), base_file_name,  "_", gene_i)
        
        #start analysis
        sink(paste0(stan_progress_dir, file_name, ".txt"))
        
        cat("Beginning Analysis.\n\n")
        
        #actually fit the model
        fit <- model$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                          data = dat_tiss, parallel_chains = nchains, adapt_delta = adapt_delta,
                          refresh = 100, max_treedepth = max_treedepth, output_dir = stan_output_dir,
                          thin = thin, init = init_mcmc)
        
        #compute mcmc diagnostics and report
        summ <- fit$summary()
        
        cat("\n\n")
        print(summ[order(summ$ess_bulk), c("variable", 
                                           "rhat", "ess_bulk", "ess_tail", 
                                           "mean", "sd")])
        
        cat("\n\n")
        print(summ[order(summ$rhat, decreasing = T), c("variable", 
                                                       "rhat", "ess_bulk", "ess_tail", 
                                                       "mean", "sd")])
        
        sink()
        
        #save output files
        save(fit, file = paste0(stan_output_dir, file_name,".cmdStanR.fit"))
        save(summ, file = paste0(stan_output_dir, file_name, ".cmdStanR.summ"))
        
        #write data to disk too for later reference
        fwrite(dt_sub, file = paste0(stan_data_dir, "genes/", file_name, ".csv"))
        write_stan_json(dat_tiss, file = paste0(stan_data_dir, "genes/", file_name, ".json"), 
                        always_decimal = FALSE)
        
      }
      
    }
    
  })
}

#### fit for GTEx OLD ####
use_rnaseq <- F
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
use_EB <- F
use_flat_MFVB <- F
init_from_prior <- F
gene_datsub <- T
model_index <- 21
lead_eQTL <- T
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
  "ase-gtex_nearly-balanced_LoH_nested_tighter-priors", #18
  "ase_param-expansion", #19
  "ase_param-expansion_1gene", #20
  "ase_param-expansion_by-gene_eQTL" #21
)
model_name <- model_names[[model_index]]
sample_prop_tail_to_remove <- 0
if(!exists("deQTL")){
  deQTL <- fread("countsGtex_lead-eqtl-genotype.csv")
  deQTL <- deQTL[deQTL$genotype == "0/1", c("tissue_site", "individual", "gene")]
  deQTL <- split(deQTL, deQTL$tissue_site)
}
model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- cmdstan_model(model_path)
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
cat(model_path)

#set MCMC inputs
sampling_param_set <- 1
nchains <- 4
niter <- c(5E3, 1E4, 2E4, 2E3)[sampling_param_set]
adapt_delta <- c(0.9, 0.8, 0.8, 0.8)[sampling_param_set]
max_treedepth <- c(10, 15, 15, 10)[sampling_param_set]
thin <- c(2, 4, 8, 1)[sampling_param_set]

init_vals_all <- list(
  log_conc_base = 1, logit_loc_base = 0, logit_ppm_base = 0,
  sd_loc_gene = 0.5,  sd_loc_indiv = 0.5,
  sd_conc_gene = 0.5, sd_conc_indiv = 0.5,
  sd_ppm_gene = 0.5, sd_ppm_indiv = 0.5,
  weight_slab_base = 0, sd_weight_slab_gene = 0.5, sd_weight_slab_indiv = 0.5,
  weight_loh_base = 0, sd_weight_loh_gene = 0.5, sd_weight_loh_indiv = 0.5
)

model_blocks <- parse_stan_blocks(clean_stan(model_lines))
model_blocks <- names(model_blocks)[sapply(model_blocks, length) != 0]
setNm <- function(x) setNames(x, x)
model_objects <- lapply(setNm(model_blocks), retrieve_block_objects, stan_code = model_lines)
init_vals <- init_vals_all[intersect(names(init_vals_all), model_objects$parameters$var_name)]

#convert model to one with fixed scale parameters
model_params <- retrieve_block_objects(model_lines, block = "parameters")$var_name
model_scale_params <- model_params[grepl("sd_", model_params)]
model_scale_param_targets <- setNames(sub("sd_", "", model_scale_params), model_scale_params)
if(use_EB){
  
  if(use_flat_MFVB){
    flattened_model <- flatten_model(model_lines, set_as = 10)
    flat_mod <- cmdstan_model(write_stan_file(flattened_model))
  }
  
  empirical_model <- flatten_model(model_lines, put_params_in_data = T)
  empirical_mod <- cmdstan_model(write_stan_file(empirical_model))
  empirical_init_vals <- init_vals[!names(init_vals) %in% model_scale_params]
  n_ADVI <- 32
  
}

tissue_names <- setNames(names(tissue_codes), tissue_codes)
simultaneous_sets <- 2
# tissue_codes_to_use <- tissue_codes[c(2,4)]
tissue_codes_to_use <- tissue_codes
# tissue_codes_to_use <- c("brain-caudate-(basal-ganglia)", "brain-cortex", "lung", "testis")
# tissue_codes_to_use <- c("lung")
tissue_sets <- suppressWarnings(split(tissue_codes_to_use, f = 1:simultaneous_sets))

for(tissue_set in tissue_sets){
  
  future({
    
    for(tiss_j in tissue_set){
      
      file_name <- paste0(model_name, "_", tiss_j, "_", sample_prop_tail_to_remove, 
                          ifelse(init_from_prior, "_prior-init", "_reg-init"),
                          ifelse(use_EB,
                                 ifelse(use_flat_MFVB, "_EB-flat", "_EB-multilevel"), 
                                 "_FB"),
                          ifelse(use_rnaseq, "_rnaseq", ""))
      files_found <- grepl(paste0(file_name, ".cmdStanR.fit"), list.files(stan_output_dir))
      dat_tiss <- jsonlite::fromJSON(paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_all.json"))
      tiss_deQTL <- deQTL[[tissue_names[tiss_j]]]
      
      # if("binom_probs" %in% model_objects$parameters$var_name){
      #   init_vals$binom_probs <- cbind(rep(0.5, dat_tiss$n), rep(0.9999, dat_tiss$n), rep(0.75, dat_tiss$n))  
      # }
      
      gene_name <- ""
      if(gene_datsub){
        
        gene_i <- 2
        dt <- as.data.frame(fread(file = paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_all.csv")))
        dt <- dt[dt$gene_i == gene_i,]
        gene_name <- paste0("_", dt$gene[1])
        
        #get het indiv
        het_indiv <- unique(tiss_deQTL$individual[tiss_deQTL$gene ==  dt$gene[1]])
        dt$eQTL_het <- dt$individual %in% het_indiv
        
        
        indivs <- unique(dt$individual)
        
        #recompile data
        dat_tiss <- list(
          n = nrow(dt),
          n_indiv = length(unique(indivs)),
          indiv_j = match(dt$individual, indivs),
          count = dt$altCount,
          total = dt$totalCount,
          eQTL_het = dt$eQTL_het + 1
        )
      }
      file_name <- paste0(file_name, gene_name)
      
      sink(paste0(stan_progress_dir, file_name, ".txt"))
      
      cat("Beginning Analysis.\n\n")
      # process Stan file for prior predictive simulation
      if(init_from_prior){
        cat("sampling from prior for chains: ")
        
        init_mcmc <- lapply(1:nchains, function(i) {
          cat(paste0(" ", i))
          r_code <- parse_Stan(model_lines, 
                               dat = dat_tiss, 
                               samps = NA, 
                               output_file = NA, 
                               sample_index = NA, 
                               post_pred_sim = F, 
                               sim = TRUE)
          
          cat(paste0("... got code..."))
          
          #evaluate code
          stan_env <- new.env()
          eval(parse(text = r_code), envir = stan_env)
          prior_predictive_sim <- stan_env$out
          
          cat(paste0("... got prior sample..."))
          
          
          return(prior_predictive_sim[model_objects$parameters$var_name])
        })
        
        cat("\nprior samples obtained!\n\n")
        
      } else {
        
        init_mcmc <- lapply(1:nchains, function(i) init_vals)
        init_mcmc <- 1
        
      }
      
      sink()
      
      
      sink(paste0(stan_progress_dir, file_name, ".txt"), append = T)
      
      if(use_EB){
        
        #going to just use standard meanfield ADVI for now?
        cat(paste0("Now performing ADVI.", "\n\n"))
        
        sink()
        
        var_outs <- parallel::mclapply(1:n_ADVI, function(iter_i) {
          
          sink(paste0(stan_progress_dir, file_name, ".txt"), append = T)
          cat(paste0(iter_i, " "))
          sink()
          
          var_out <- tryCatch({
            if(use_flat_MFVB){
              flat_mod$variational(data = dat_tiss, init = 0.1,
                                   algorithm = "meanfield",
                                   iter = 1E5, output_samples = 5E3)
            } else {
              mod$variational(data = dat_tiss,
                              algorithm = "meanfield",
                              iter = 1E5, output_samples = 5E3)
            }
            
          },error = function(cond) {NA},
          warning = function(cond) {NA})
          
          return(var_out)
          
        }, mc.cores = nchains)
        
        
        #resume progress file
        sink(paste0(stan_progress_dir, file_name, ".txt"), append = T)
        cat(paste0("Finished ADVI.", "\n"))
        
        var_outs <- var_outs[unlist(lapply(seq_along(var_outs), function(i){
          tryCatch("draws_matrix" %in% class(var_outs[[i]]$draws()), 
                   error = function(cond) F, warning = function(cond) F)
        }))]
        save(var_outs, file = paste0(stan_output_dir, file_name, "_MF-VarInf-Runs.RData"))
        
        #report success rate
        cat(paste0(length(var_outs), " / ", n_ADVI, "runs were successful.\n"))
        
        #select best fit and extract scale parameters
        best_fit_index <- which.max(lapply(seq_along(var_outs), function(i) var_outs[[i]]$summary(variables = "lp_approx__")$mean))
        var_out <- var_outs[[best_fit_index]]
        
        if(use_flat_MFVB){
          
          #extract target param samples
          samps_mfvb <- as.data.table(data.frame(as_draws_df(var_out$draws(model_scale_param_targets))))
          
          #summarize flat fit scale empirically
          param_mean_var_mfvb <- sapply(model_scale_param_targets, function(scale_param){
            param_samps <- do.call(rbind, munge_samps(scale_param, subset_samps(scale_param, samps_mfvb)))
            param_vars <- apply(param_samps, 1, var)
            mean(param_vars)
          })
          param_sds <- sqrt(param_mean_var_mfvb)
          
        } else {
          
          param_sds <- var_out$summary(variables = model_scale_params)
          param_sds <- setNames(param_sds$mean, model_scale_params)
          
        }
        
        #now fit the final model via MCMC
        empirical_dat_tiss <- c(dat_tiss, param_sds)
        out <- empirical_mod$sample(chains = nchains, iter_sampling = 3E3, iter_warmup = 3E3, 
                                    data = empirical_dat_tiss, parallel_chains = nchains, adapt_delta = 0.9, 
                                    refresh = 100, max_treedepth = 10, output_dir = stan_output_dir,
                                    thin = 2, threads_per_chain = 1, 
                                    init = lapply(1:nchains, function(i) empirical_init_vals))
        
      } else {
        
        out <- mod$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                          data = dat_tiss, parallel_chains = nchains, adapt_delta = adapt_delta,
                          refresh = 100, max_treedepth = max_treedepth, output_dir = stan_output_dir,
                          thin = thin, init = init_mcmc)
      }
      
      summ <- out$summary()
      cat("\n\n")
      print(summ[order(summ$ess_bulk), c("variable", 
                                         "rhat", "ess_bulk", "ess_tail", 
                                         "mean", "sd")])
      cat("\n\n")
      print(summ[order(summ$rhat, decreasing = T), c("variable", 
                                                     "rhat", "ess_bulk", "ess_tail", 
                                                     "mean", "sd")])
      sink()
      
      save(out, file = paste0(stan_output_dir, file_name,".cmdStanR.fit"))
      save(summ, file = paste0(stan_output_dir, file_name, ".cmdStanR.summ"))
      
    }
    
  })
}


#### fit for OVC-vs-GTEx ####

#load in model
use_rnaseq <- F
base_fit <- T
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
stan_data_gene_dir <- paste0(stan_data_dir, "genes/")

model_index <- 2
model_names <- list(
  "ovc-ase_param-expansion_by-gene_raw-params_resid-indivs", #1
  "ovc-ase_param-expansion_by-gene_raw-params_resid-indivs_empirical" #2
)

model_name <- model_names[[model_index]]
sample_prop_tail_to_remove <- 0
model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- cmdstan_model(model_path)
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
cat(model_path)

#load in data
data_index <- 1
well_sampled_genes <- F
remove_dupe_hetsites <- F
min_n_per_grade <- 3
data_names <- list(
  "OVC-vs-GTEx_0", #1
  "OVC-vs-GTEx_perm_0", #2
  "OVC-vs-GTEx_sim_0", #3
  "OVC-vs-GTEx_0_all" #4
)
data_name <- paste0(data_names[[data_index]],
                    ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                    ifelse(use_rnaseq, "_rnaseq", ""),
                    ifelse(remove_dupe_hetsites, "", "_all"))
dat <- jsonlite::fromJSON(paste0(stan_data_dir, data_name, ".json"))
dt <- fread(paste0(stan_data_dir, data_name, ".csv"))
genes <- fread(paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove,
       ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
       "_gene-map.csv"))$gene_symbol

#set MCMC inputs
sampling_param_set <- 2
nchains <- 4
niter <- c(5E3, 1E4, 2E4, 2E3)[sampling_param_set]
adapt_delta <- c(0.9, 0.95, 0.98, 0.8)[sampling_param_set]
max_treedepth <- c(10, 12, 15, 10)[sampling_param_set]
thin <- c(2, 4, 8, 1)[sampling_param_set]

model_blocks <- parse_stan_blocks(clean_stan(model_lines))
model_blocks <- names(model_blocks)[sapply(model_blocks, length) != 0]
setNm <- function(x) setNames(x, x)
model_objects <- lapply(setNm(model_blocks), retrieve_block_objects, stan_code = model_lines)

future({

  base_file_name <- paste0(model_name, "_", sample_prop_tail_to_remove, 
                           ifelse(use_rnaseq, "_rnaseq", ""))
  
  for(gene_i in genes[1:20]){
    
    file_name <- paste0(ifelse(base_fit, "1st-fit_", ""), base_file_name,  "_", gene_i)  
    data_path <- paste0(stan_data_gene_dir, "OVC-vs-GTEx", "_", 
                        sample_prop_tail_to_remove,
                        ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""),
                        ifelse(remove_dupe_hetsites, "", "_all"), 
                        "_", gene_i, ".json")
    gdat <- jsonlite::fromJSON(data_path)

    sink(paste0(stan_progress_dir, file_name, ".txt"))
    
    cat("Beginning Analysis.\n\n")
    init_mcmc <- 1
    
    #actually fit the model
    fit <- mod$sample(chains = nchains, iter_sampling = niter/2, iter_warmup = niter/2,
                        data = gdat, parallel_chains = nchains, adapt_delta = adapt_delta,
                        refresh = 100, max_treedepth = max_treedepth, output_dir = stan_output_dir,
                        thin = thin, init = init_mcmc)
    
    #compute mcmc diagnostics and report
    summ <- fit$summary()
    
    cat("\n\n")
    print(summ[order(summ$ess_bulk), c("variable", 
                                       "rhat", "ess_bulk", "ess_tail", 
                                       "mean", "sd")])
    
    cat("\n\n")
    print(summ[order(summ$rhat, decreasing = T), c("variable", 
                                                   "rhat", "ess_bulk", "ess_tail", 
                                                   "mean", "sd")])
    
    sink()
    
    #save output files
    save(fit, file = paste0(stan_output_dir, file_name,".cmdStanR.fit"))
    save(summ, file = paste0(stan_output_dir, file_name, ".cmdStanR.summ"))
    
    #write data to disk too for later reference
    fwrite(dt_sub, file = paste0(stan_data_dir, "genes/", file_name, ".csv"))
    write_stan_json(dat_tiss, file = paste0(stan_data_dir, "genes/", file_name, ".json"), 
                    always_decimal = FALSE)
    
  }
  
})


#### fit for OVC-vs-GTEx OLD ####

#just to the gtex <-> ovc diff params for now

#thought: maybe I can accommodate pseudoreplication by iterating through all
#the different het sites for each individual x gene, and then average over them 
#when computing the likelihood in proportion to their relative #  of observations?
#for that individual x gene

#set mcmc params
sampling_param_set <- 4
nchains <- 4
niter <- c(5E3, 1E4, 2E4, 2E3)[sampling_param_set]
adapt_delta <- c(0.8, 0.9, 0.95, 0.8)[sampling_param_set]
max_treedepth <- c(10, 12, 15, 10)[sampling_param_set]
thin <- c(2, 4, 8, 1)[sampling_param_set]

#load in the model
use_rnaseq <- F
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
model_index <- 27
model_names <- list(
  "beta-binomial-joint-conc", #1
  "beta-binomial-joint-conc-loc", #2
  "beta-binomial-joint-conc-loc-noncentered", #3
  "beta-binomial-joint-conc-loc-noncentered-simple", #4, does not fit
  "beta-binomial-joint-loc", #5
  "beta-binomial-joint-conc-loc-indiv", #6
  "beta-binomial-simple-mixture", #7
  "beta-binomial-simplest-mixture", #8
  "beta-binomial-simpler-mixture", #9
  "beta-binomial-joint-conc-loc-indiv-pointmix", #10
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered", #11
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_grade", #12
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene", #13
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid", #14
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid-identifiability", #15
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid_nearly-balanced", #16
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_cumulative_nearly-balanced", #17
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade", #18
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc", #19
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc-correct_grade_REfs", #20
  "ovc-ase-gtex_IID_nearly-balanced_LoH_nested-tumor_grade", #21
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc_wider-init", #22
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc_IID-Diff-Split", #23
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc_wider-init_cumul-add", #24
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-less_inform_conc", #25
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-less_inform_conc-no_indiv", #26
  "ovc-ase_param-expansion" #27
)

model_name <- model_names[[model_index]]
model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- cmdstan_model(model_path, cpp_options = list(STAN_THREADS = T, stan_threads = TRUE))
model_string <- paste0(readLines(model_path), collapse = "\n")
cat(model_path)
# cat(model_string)

#load in the data
data_index <- 1
well_sampled_genes <- T
remove_dupe_hetsites <- F
min_n_per_grade <- 3
data_names <- list(
  "OVC-vs-GTEx_0", #1
  "OVC-vs-GTEx_perm_0", #2
  "OVC-vs-GTEx_sim_0", #3
  "OVC-vs-GTEx_0_all" #4
)
data_name <- paste0(data_names[[data_index]],
                    ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                    ifelse(use_rnaseq, "_rnaseq", ""),
                    ifelse(remove_dupe_hetsites, "", "_all"))
dat <- jsonlite::fromJSON(paste0(stan_data_dir, data_name, ".json"))

#specify chain parameters
init_vals_all <- list(log_conc_base = 2, logit_loc_base = 0,
                      logit_ppm_base = 0,
                      
                      sd_loc_gene = 0.5, sd_loc_grade = 0.5, sd_loc_gtex = 0.5,
                      sd_loc_indiv = 0.5, sd_loc_gene_x_indiv = 0.5,
                      
                      sd_conc_gene = 0.5, sd_conc_grade = 0.5, sd_conc_gtex = 0.5,
                      sd_conc_indiv = 0.5, sd_conc_gene_x_indiv = 0.5,
                      
                      sd_ppm_gene = 0.5, sd_ppm_grade = 0.5, 
                      sd_ppm_gtex = 0.5, sd_ppm_indiv = 0.5,
                      
                      mu_loc_grade = 0, mu_conc_grade = 0, mu_ppm_grade = 0,
                      loc_gtex = 0, conc_gtex = 0, ppm_gtex = 0,
                      
                      sd_mu_conc_grade_gene = 0.5, sd_mu_loc_grade_gene = 0.5, sd_mu_ppm_grade_gene = 0.5
                      
)
params_present <- sapply(paste0(names(init_vals_all), ";"), function(param_name)
  grepl(pattern = param_name, x = model_string))
init_vals <- init_vals_all[params_present]

# files_found <- grepl(paste0(model_name, "-[0-9]+"), list.files(stan_output_dir))
files_found <- grepl(paste0(model_name, "_", data_name, ".cmdStanR.fit"), list.files(stan_output_dir))

# if(!any(files_found)){


future({
  
  sink(paste0(stan_progress_dir, model_name, "_", data_name, ".txt"))
  out <- mod$sample(chains = nchains, iter_sampling = niter, iter_warmup = niter, 
                    data = dat, parallel_chains = nchains, adapt_delta = adapt_delta, 
                    refresh = 20, max_treedepth = max_treedepth, output_dir = stan_output_dir,
                    thin = thin, threads_per_chain = 1, 
                    init = 1)
  summ <- out$summary()
  cat("\n\n")
  print(summ[order(summ$ess_bulk), c("variable", 
                                     "rhat", "ess_bulk", "ess_tail", 
                                     "mean", "sd")])
  cat("\n\n")
  print(summ[order(summ$rhat, decreasing = T), c("variable", 
                                                 "rhat", "ess_bulk", "ess_tail", 
                                                 "mean", "sd")])
  sink()
  
  save(out, file = paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.fit"))
  save(summ, file = paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.summ"))
  
})


