library(rstan)
library(data.table)

# specify number of cores
num_cores <- 128
options(mc.cores = num_cores)


#### functions ####
setNm <- function(x) setNames(x, x)

#### create useful variables ####
stan_model_dir <- "Stan/models/"
stan_data_dir <- "Stan/data/"
stan_output_dir <- "Stan/output/genes/"
stan_progress_dir <- "Stan/progress/genes/"
tissue_codes_df <- fread("tissue_codes.txt")
tissue_codes <- setNames(tissue_codes_df$code, tissue_codes_df$name)
tissue_names <- setNames(names(tissue_codes), tissue_codes)

#### fit for GTEx ####

#specify fitting settings
base_fit <- T
redo_failed_chains <- F
use_rnaseq <- F
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/mmpcrseq/"
}
stan_data_genes_dir <- paste0(stan_data_dir, "genes/")

#load in the model
model_index <- 1
model_names <- list(
  "ovc-ase_param-expansion_by-gene_raw-params_resid-indivs" #1
)

model_name <- model_names[[model_index]]
sample_prop_tail_to_remove <- 0
model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- stan_model(file = model_path)

#set MCMC inputs
sampling_param_set <- 2
nchains <- 4
niter <- c(5E3, 1E4, 2E4, 2E3)[sampling_param_set]
adapt_delta <- c(0.9, 0.95, 0.95, 0.8)[sampling_param_set]
max_treedepth <- c(10, 12, 15, 10)[sampling_param_set]
thin <- c(2, 4, 8, 1)[sampling_param_set]
init_mcmc <- 1

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
base_data_name <- paste0(stan_data_genes_dir, "OVC-vs-GTEx", "_", 
                         sample_prop_tail_to_remove,
                         ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""),
                         ifelse(remove_dupe_hetsites, "", "_all"))
base_file_name <- paste0(model_name, "_", sample_prop_tail_to_remove, 
                         ifelse(use_rnaseq, "_rnaseq", ""))
genes <- fread(paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove,
                      ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                      "_gene-map.csv"))$gene_symbol

#start analysis
parallel::mclapply(genes, function(gene_i){

  #read in the data
  data_name <- paste0(base_data_name, "_", gene_i, ".json")
  dat <- jsonlite::fromJSON(data_name)
  
  #skip genes with only one individual 
  if(dat$n_indiv == 1){
    next()
  }
  
  #specify output name structure
  file_name <- paste0(ifelse(base_fit, "1st-fit_", ""), base_file_name,  "_", gene_i)
  
  #start analysis
  sink(paste0(stan_progress_dir, file_name, ".txt"))
  
  cat("Beginning Analysis.\n\n")
  
  #actually fit the model
  fit <- sampling(
    object = model,
    data = dat,
    chains = nchains,
    iter = niter,
    warmup = niter / 2,
    thin = thin,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    refresh = 100,
    init = init_mcmc,
    cores = 4
  )
  
  #compute mcmc diagnostics and report
  summ <- data.frame(summary(fit)$summary)
  summ$variable <- rownames(summ)
  
  cat("\n\n")
  print(head(summ[order(summ$n_eff), c("variable", 
                                       "Rhat", "n_eff", 
                                       "mean", "sd")]))
  
  cat("\n\n")
  print(head(summ[order(summ$Rhat, decreasing = T), c("variable", 
                                                      "Rhat", "n_eff", 
                                                      "mean", "sd")]))
  
  sink()
  
  #save output files
  save(fit, file = paste0(stan_output_dir, file_name,".stanfit"))
  fwrite(summ, file = paste0(stan_output_dir, file_name, ".summ"))
  
}, mc.cores = 32)
