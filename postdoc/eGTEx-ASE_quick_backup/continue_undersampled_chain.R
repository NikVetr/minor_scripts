setwd("~/repos/egtex-ase")
sample_prop_tail_to_remove <- 0

#### Load Packages ####
library(abind)
library(cmdstanr)
library(future)
library(posterior)
library(data.table)

#### Set up future multisession ####
if(!inherits(plan(), "multisession")){
  plan(multisession)  
}

#### specify functions ####
source("scripts/functions.R")
source("~/repos/Stan2R/R/functions.R")
setNm <- function(x) setNames(x,x)

#### inspect MCMC from OVC-vs-GTEx ####
stan_model_dir <- "Stan/models/"
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


#load in MCMC i/o
#load in data
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
dat <- jsonlite::fromJSON(paste0(stan_data_dir, data_name, ".json"))
model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- cmdstan_model(model_path, cpp_options = list(STAN_THREADS = T, stan_threads = TRUE))
load(paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.fit"))
fit <- out

# load(paste0(stan_output_dir, model_name, "_", data_name, ".cmdStanR.summ"))
# print(summ[order(summ$ess_bulk), c("variable", 
#                                    "rhat", "ess_bulk", "ess_tail", 
#                                    "mean", "sd")])
# print(summ[order(summ$rhat, decreasing = T), c("variable", 
#                                                "rhat", "ess_bulk", "ess_tail", 
#                                                "mean", "sd")])


#run continuation
sampling_param_set <- 1
nchains <- 4
niter <- c(5E3, 1E4, 2E4)[sampling_param_set]
adapt_delta <- c(0.9, 0.8, 0.8)[sampling_param_set]
max_treedepth <- c(10, 15, 15)[sampling_param_set]
thin <- c(2, 4, 8)[sampling_param_set]

fit_metadata <- fit$metadata()
step_size <- fit_metadata$step_size_adaptation
last_draws <- fit$draws(variables = NULL) # Extract all draws
init_values <- lapply(1:dim(last_draws)[2], function(chain) {
  as.list(last_draws[1, chain, , drop = TRUE])
})

efit <- mod$sample(chains = nchains, iter_sampling = niter, iter_warmup = 0, 
                      data = dat, parallel_chains = nchains, adapt_delta = adapt_delta, 
                      refresh = 100, max_treedepth = max_treedepth, output_dir = stan_output_dir,
                      thin = thin, threads_per_chain = 1, step_size = step_size, 
                      init = init_values)

save(efit, file = paste0(stan_output_dir, model_name, "_", data_name, "_extra.cmdStanR.fit"))


