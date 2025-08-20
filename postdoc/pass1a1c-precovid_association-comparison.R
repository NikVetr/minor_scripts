#### packages ####
library(MotrpacHumanPreSuspensionData)
library(MotrpacRatTraining6moData)
library(MotrpacBicQC)
library(parallel)
library(data.table)
library(cmdstanr)
library(posterior)

#### setting the stage ####

#bayesian model paths
model_dir <- "~/scripts/minor_scripts/postdoc/"
model_name <- "bivariate_correlation_uncertainty_fast_marginalize-out-latents.stan"
# model_name <- "bivariate-t_correlation_uncertainty_fast.stan"
model_name <- "bivariate_var-infl-normal_correlation_uncertainty_fast_marginalize-out-latents_vectorized.stan"
model_path <- paste0(model_dir, model_name)
model_string <-  paste0(readLines(model_path), collapse = "\n")
mod <- cmdstan_model(model_path)

#number of bootstrap replicates for weighted correlation
nboot <- 2E3

#gene ortholog matching
rgd_orthologs <- fread("~/data/smontgom/RGD_ORTHOLOGS_20201001.txt", header = T, sep = "\t")
gencode_gene_map <- rdg_mapping <- fread("~/data/smontgom/gencode.v39.RGD.20201001.human.rat.gene.ids.txt")
gencode_gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID <- gsub(gencode_gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID, pattern = "\\..*", replacement = "")
gene_map <- gencode_gene_map

#rat mapping
ome_code_map <- c(
  "metab"              = "metab-merged-all",   # use full merged table
  "prot-ph"            = "prot-ph",            # phosphoproteome
  "prot-pr"            = "prot-pr",            # global proteome
  "transcript-rna-seq" = "transcript-rna-seq"  # transcriptomics
)
tissue_code_map <- list(
  adipose = list(
    default            = "t70-white-adipose",   # choose WAT‑SC by default
    alternate      = "t69-brown-adipose"
  ),
  blood = list(
    default            = "t30-blood-rna",       # whole blood transcriptome
    alternate             = "t31-plasma"           # plasma (metabolomics)
  ),
  muscle = list(
    default            = "t55-gastrocnemius",   # gastrocnemius
    alternate             = "t56-vastus-lateralis"
  )
)
ome_folder_map <- c(
  "metab"              = "metabolomics",
  "prot-ph"            = "proteomics",
  "prot-pr"            = "proteomics",
  "transcript-rna-seq" = "transcriptomics"
)

#point to the authenticated gsutil
good_gsutil <- normalizePath("~/google-cloud-sdk/bin/gsutil")
Sys.setenv(
  PATH = paste(dirname(good_gsutil), Sys.getenv("PATH"), sep = .Platform$path.sep)
)
Sys.setenv(CLOUDSDK_CONFIG = path.expand("~/.config/gcloud"))

#function to retrieve rat data
get_rat_da <- function(tissue, ome,
                       contrast        = c("exercise", "control-only")[1],
                       tissue_variant  = "default") {
  
  tc <- tissue_code_map[[tissue]][[tissue_variant]]
  if (is.null(tc)) stop("Unknown tissue or variant")
  oc <- ome_code_map[[ome]]
  of <- ome_folder_map[[ome]]
  
  # wildcard * swallows the middle chunk (norm & model tags)
  pattern <- sprintf(
    "gs://motrpac-data-hub/analysis/rat-acute-06/%s/da/rat-acute-06_%s_%s_da_*-%s_v1.0.txt",
    of, tc, oc, contrast
  )
  
  # there should be exactly one match
  file <- system(sprintf("gsutil ls %s", shQuote(pattern)), intern = TRUE)
  if (!length(file)) stop("No file matched: ", pattern)
  
  message("Downloading: ", file)
  dl_read_gcp(file, gsutil_path = Sys.which("gsutil"), tmpdir = tempdir())
}

#metadata
tissues <- c("adipose", "blood", "muscle")
omes <- c("metab", "prot-ph", "prot-pr", "transcript-rna-seq")[4]
tissue_variants <- c("default", "alternate", "permuted", "permuted_alternate")

#load all human data to subset later
if(!exists("hda")){
  hda <- load_differential_analysis()  
}

#pre-download Rat effects
cache_dir <- "~/data/smontgom/rat_acute_DA_cache"          # pick any writable folder
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

rat_path <- function(tissue, ome,
                     contrast       = c("exercise", "control-only")[1],
                     tissue_variant = "default") {
  
  tc <- tissue_code_map[[tissue]][[tissue_variant]]
  if (is.null(tc))           stop("unknown tissue / variant")
  
  oc <- ome_code_map[[ome]]
  of <- ome_folder_map[[ome]]
  if (is.null(oc) || is.null(of))
    stop("unknown ome label")
  
  method_tag <- switch(
    ome,
    "metab"              = "timewise-limma-t2f",
    "prot-ph"            = "timewise-limma-bc-t2f",
    "prot-pr"            = "timewise-limma-bc-t2f",
    "transcript-rna-seq" = "timewise-limma-bc-t2f"
  )
  
  sprintf(
    "gs://motrpac-data-hub/analysis/rat-acute-06/%s/da/rat-acute-06_%s_%s_da_%s-%s_v1.0.txt",
    of, tc, oc, method_tag, contrast
  )
}

combos <- expand.grid(tissue         = tissues,
                      tissue_variant = tissue_variants,
                      ome            = omes,
                      stringsAsFactors = FALSE)

download_one <- function(idx) {
  row <- combos[idx, ]
  
  remote <- rat_path(row$tissue, row$ome,
                     contrast       = "exercise",
                     tissue_variant = row$tissue_variant)
  
  # confirm remote exists (1 network roundtrip)
  if (system(sprintf("%s ls %s", shQuote(good_gsutil), shQuote(remote))) != 0){
    warning("Remote file not found: ", remote, ". Skipping...")
    return()  
  }
  
  local <- normalizePath(                 # expands "~", keeps path otherwise
    file.path(cache_dir, basename(remote)),
    mustWork = FALSE
  )
  
  dir.create(dirname(local), showWarnings = FALSE, recursive = TRUE)
  
  if (!file.exists(local) || file.info(local)$size == 0) {
    cmd <- sprintf("%s cp %s %s",
                   shQuote(good_gsutil), shQuote(remote), shQuote(local))
    message("Downloading → ", local)
    status <- system(cmd)
    if (status != 0) stop("gsutil cp failed for: ", remote)
  } else {
    message("✔ found cached file: ", local)
  }
  return(local)
}

cache_paths <- unlist(
  lapply(seq_len(nrow(combos)),
           download_one)
)

load_rat_da_if <- function(tissue, ome,
                           tissue_variant = "default",
                           contrast       = "exercise") {
  
  local_file <- file.path(
    cache_dir,
    basename(rat_path(tissue, ome, contrast, tissue_variant))
  )
  
  if (!file.exists(local_file)) return(NULL)        # ← graceful miss
  data.table::fread(local_file, data.table = FALSE)
}

#### analyze ####
# human + rat DA results comparison
results <- list()
for(tissue_variant in tissue_variants){
  for(tissue in tissues){
    for(ome in omes){
      
      #download rat data
      if(tissue_variant == "permuted"){
        rat_tissue <- "permuted"
        rda_sub <- load_rat_da_if(tissue, ome,
                                  tissue_variant = "default",
                                  contrast       = "exercise")
      } else if(tissue_variant == "permuted_alternate"){
        rat_tissue <- "permuted_alternate"
        rda_sub <- load_rat_da_if(tissue, ome,
                                  tissue_variant = "alternate",
                                  contrast       = "exercise")
      } else {
        rat_tissue <- tissue_code_map[[tissue]][[tissue_variant]]
        rda_sub <- load_rat_da_if(tissue, ome,
                                  tissue_variant = tissue_variant,
                                  contrast       = "exercise")
      }
      
      if (is.null(rda_sub)) {
        message("⏭  no rat file for ",
                paste(tissue, tissue_variant, ome, sep = " / "))
        next
      }
      
      #get df from rat data
      f <- function(d, t, p) 2 * pt(-abs(t), d) - p
      rda_sub$df <- unlist(mclapply(seq_len(nrow(rda_sub)), function(i) {
        t <- rda_sub$t_stat[i]
        p <- rda_sub$p_value[i]
        f1 <- f(1, t, p)
        f2 <- f(1e6, t, p)
        if (is.na(f1) || is.na(f2) || sign(f1) == sign(f2)) {
          return(NA_real_)  # or log the issue
        }
        uniroot(function(d) f(d, t, p), interval = c(1, 1e6))$root
      }, mc.cores = 12L))
      
      #subset human results
      hda_sub <- hda[[tissue]][[ome]]
      
      #check consistency among SEs recovered in two ways for humans
      CI_SEs <- (hda_sub$CI.R - hda_sub$CI.L) / (2 * qt(0.975, df = hda_sub$degrees_of_freedom))
      # t_SEs <- hda_sub$logFC / hda_sub$t
      # n <- nrow(hda_sub)
      # inds <- seq(1, n, length.out = 2E3)
      # plot(CI_SEs[inds], t_SEs[inds]); abline(0,1, col = 2)
      # hist(CI_SEs - t_SEs)
      # quantile(CI_SEs - t_SEs)
      # rel_diff <- (CI_SEs - t_SEs) / CI_SEs
      # summary(rel_diff)
      #close enough
      hda_sub$SE <- CI_SEs
      
      #### comparison ####
      
      #match human and rat genes (or other analytes)
      hda_sub$gene_id <- gsub("\\..*$", "", hda_sub$feature_id)
      rda_sub$human_gene_id <- gene_map$HUMAN_ORTHOLOG_ENSEMBL_ID[match(rda_sub$feature_id, gene_map$RAT_ENSEMBL_ID)]
      rda_sub <- rda_sub[!is.na(rda_sub$human_gene_id),]
      matched_genes <- intersect(rda_sub$human_gene_id, hda_sub$gene_id)
      
      rda_dat <- rda_sub[match(matched_genes, rda_sub$human_gene_id),]
      hda_dat <- hda_sub[match(matched_genes, hda_sub$gene_id),]
      
      #weighted corr
      # weights_gfx <- rowMeans(1/cbind(rda_dat$logFC_se^2, hda_dat$SE^2))
      # weights_fgx <- 1/rowMeans(cbind(rda_dat$logFC_se^2, hda_dat$SE^2))
      # wCorr::weightedCorr(rda_dat$logFC, hda_dat$logFC, weights = weights_gfx, method = "Pearson")
      # wCorr::weightedCorr(rda_dat$logFC, hda_dat$logFC, weights = weights_fgx, method = "Pearson")
      
      #### bayesian disattenuated corr ####
      
      #inverse-variance weighted correlation first?
      x_err <- cbind(rda_dat$logFC, hda_dat$logFC)
      sd_x_err <- cbind(rda_dat$logFC_se, hda_dat$SE)
      df <- cbind(rda_dat$df, hda_dat$degrees_of_freedom)
      
      if(grepl("permuted", tissue_variant)){
        rda_inds <- sample(1:nrow(rda_dat))
        x_err <- cbind(rda_dat$logFC[rda_inds], hda_dat$logFC)
        sd_x_err <- cbind(rda_dat$logFC_se[rda_inds], hda_dat$SE)
        df <- cbind(rda_dat$df[rda_inds], hda_dat$degrees_of_freedom)
      }
      
      iv_weights <- 1 / rowMeans(sd_x_err^2)
      
      # inverse-variance-weighted Pearson 
      cor_weighted <- wCorr::weightedCorr(x_err[,1], x_err[,2],
                                          weights = iv_weights,
                                          method  = "Pearson")
      cor_uw <- cor(x_err)[1,2]
      cor_weighted_boot <- replicate(n = nboot, expr = {
        idx <- sample(1:nrow(x_err), replace = T)
        cor_weighted <- wCorr::weightedCorr(x_err[idx,1], x_err[idx,2],
                                            weights = iv_weights[idx],
                                            method  = "Pearson")
        cor_weighted
      }, simplify = T)
      
      dat <- list(p=length(matched_genes), 
                  x_err=x_err, 
                  sd_x_err=sd_x_err,
                  df=df)
      
      #### fit bayesian model ###
      fit <- mod$sample(chains = 4, iter_sampling = 1E3, iter_warmup = 2E3,
                        data = dat, adapt_delta = 0.9, parallel_chains = 4,
                        refresh = 10, max_treedepth = 12, 
                        thin = 1, init = 0.1)
      
      #### mcmc diagnostics ###
      
      #check convergence
      check_diag <- T
      if(check_diag){
        summ <- fit$summary("r")
        print(summ[order(summ$ess_bulk),])
        print(summ[order(summ$rhat, decreasing = T),])
      }
      
      #extract samples and record
      samps_r <- data.frame(as_draws_df(fit$draws("r")))$r
      out <- list(human_tissue = tissue,
                  rat_tissue = rat_tissue,
                  ome = ome, 
                  sample_cor = cor_uw,
                  ip_wcor = cor_weighted,
                  ip_wcor_boot = cor_weighted_boot,
                  bayes_cor_samps = samps_r,
                  dat = dat,
                  fit = fit,
                  rda_dat = rda_dat,
                  hda_dat = hda_dat,
                  genes = matched_genes)
      
      results[[tissue]][[rat_tissue]][[ome]] <- out
    } 
  }
}

#### plot comparison of results ####
human_tissues <- c("adipose", "blood", "muscle")
rat_tissue_types <- c("default", "alternate", 
                     "permuted", "permuted_alternate")
omes <- c("metab", "prot-ph", "prot-pr", "transcript-rna-seq")[4]

par(mfrow = c(3,2))
for(human_tissue in human_tissues){
  for(rat_tissue_type in rat_tissue_types[4]){
    rat_tissue <- switch(
      rat_tissue_type,
      "permuted" = "permuted",
      "permuted_alternate" = "permuted_alternate",
      tissue_code_map[[human_tissue]][[rat_tissue_type]]
    )
    for(ome in omes){
      
      #recover results
      res <- results[[human_tissue]][[rat_tissue]][[ome]]
      
      if(is.null(res)){
        next()
      }
      
      #sample corr
      scor <- res$sample_cor
      #precision-weighted corr
      pcor <- res$ip_wcor
      #posterior mean
      pmean <- mean(res$bayes_cor_samps)
      #observations
      hobs <- res$hda_dat$logFC
      robs <- res$rda_dat$logFC
      nplot <- min(2E3, res$dat$p)
      sub_idx <- round(seq(1, res$dat$p, length.out = nplot))
      
      title_string <- paste0("human tissue = ", human_tissue,
                             ", rat tissue = ", rat_tissue,
                             ",\n-ome = ", ome)
      plot(hobs[sub_idx], robs[sub_idx], col = adjustcolor(1, 0.25),
           xlab = "human DE", ylab = "rat DE", main = title_string)
      hist(results[[human_tissue]][[rat_tissue]][[ome]]$bayes_cor_samps,
           xlab = "correlation", main = title_string,
           xlim = c(-1,1), breaks = 10, freq = F)
      abline(v = scor, col = 2)
      abline(v = pcor, col = 3)
    }
  }
}

