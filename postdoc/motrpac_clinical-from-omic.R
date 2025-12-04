#NOTE: should have lots of muscle proteomics and phosphoproteomics
#(very little in adipose and blood tho)

#load libraries
library(uniLasso)
library(glmnet)
library(MotrpacHumanPreSuspensionData)
library(data.table)
library(kernlab)
library(parallel)

#set seed for nested cv split across available data for different methods
seed <- 1
set.seed(seed)

#### directories ####
project_dir <- "~/repos/motrpac-predict/"
results_dir <- paste0(project_dir, "results/")
figures_dir <- paste0(project_dir, "figures/")
subfigures_dir <- paste0(figures_dir, "subfigures/")

#### functions ####
source("~/scripts/minor_scripts/postdoc/motrpac_clinical-x-omic_functions.R")
source("~/repos/flexible-sum-to-zero/Nik/R/functions.R")

#### analysis params ####
phenotypic_predictors <- T
molecular_predictors <- T
sex_interactions <- F
residualize_molecular_feats <- F
do_ISIS_hclust <- T
incl_latent_variables <- F
use_pf <- F #use pathfinder instead of full MCMC
log_transform <- F
unnorm_VO2max <- T
linear_recalibration <- F
winsorize <- F
use_median_mad <- F
target_phenotype <- "vo2_peak_mlkg"
# target_phenotype <- "sex"
phenotype_predictors <- c("sex", "age_years", 
                          "body_weight", "body_height_cm", 
                          "bmi", "bp_sys", "bp_dia", "pulse_rest")
phenotype_predictors <- setdiff(phenotype_predictors, target_phenotype)

# c("lean_mass", "fat_mass", "total_bmd", "vat_mass")

n_mol_target <- 1E3 #set to high number to use all of them
n_mol_target_ISIS <- 5E3 #needs to be same size or larger than n_mol_target
n_mol_target_ISIS <- max(n_mol_target, n_mol_target_ISIS)
methods <- c("lasso", "unilasso", "horseshoe")[3]
tissues <- c("blood", "muscle", "adipose")[2]
ome_opts <- c("transcript", "prot", "metab")
ome_sets <- list(
  ome_opts[c(1,3)],
  ome_opts[1],
  ome_opts[3]
)[3]


mcmc_params <- list( #for the Bayesian models  
  n_iter = 1E3,
  adapt_delta = 0.98,
  max_treedepth = 12
)
use_all_for_PCA <- T
standardize_before_PCA <- T
ptrain <- 0.8 #train proportion (1-ptrain is test proportion)
kfold_test <- T #run k-fold test / validation set? 
nfolds <- ifelse(kfold_test, round(1 / (1-ptrain)), 1)

#hs params
hs_pred_type    <- c("pp_mean", "plugin")[1]
hs_estimator    <- c("mean", "median")[1]
hs_beta_thresh  <- 1e-2       # threshold to call a coefficient "non-zero"
use_annotation <- T

# load molecular data

#or get it another way #NO this data is no more or less up-to-date
# data_dir <- "~/data/human-precovid-sed-adu/"
# dir.create(data_dir, showWarnings = FALSE, recursive = T)
# system(paste0("gsutil -m cp -r gs://motrpac-data-hub/analysis/human-precovid-sed-adu/* ", data_dir))
# list.files(data_dir)


#subset molecular data
# molecular_data <- load_qc_local(data_dir, tissue, ome)

# ome <- "all"

fig_paths <- character()
for(method in methods){
  for(tissue in tissues){
    for(omes in ome_sets){

      #### get / prep data ####
      md_omes <- list()
      a2s_keys <- list()
      for(ome in omes){
        molecular_data <- load_qc()
        matching_omes <- names(molecular_data[[tissue]])[grep(ome, names(molecular_data[[tissue]]))]
        analyte2subome <- lapply(molecular_data[[tissue]][matching_omes], 
                                 function(oi) rownames(oi$qc_norm))
        # str_lol(molecular_data, n_names = 20)
        if(!(ome %in% names(molecular_data[[tissue]]))){
          
          #combine all metabolites into one data frane
          mds <- lapply(setNames(matching_omes, matching_omes), function(oi){
            molecular_data[[tissue]][[oi]]$qc_norm
          })
          
          # convert each df to data.table, lift rownames to a column "analyte"
          mds_dt <- lapply(seq_along(mds), function(i) {
            df <- mds[[i]]
            rn <- rownames(df)
            setDT(df, keep.rownames = "analyte")
            df[ , src := i]  # optional: keep source index for diagnostics
            if (!is.null(rn)) df$analyte <- rn
            df
          })
          combined_dt <- rbindlist(mds_dt, use.names = TRUE, fill = TRUE)
          md <- as.data.frame(combined_dt)
          rownames(md) <- md$analyte
          md$analyte <- NULL
          md$src <- NULL
          
          mds_meta <- lapply(setNames(matching_omes, matching_omes), function(oi) 
            cbind(molecular_data[[tissue]][[oi]]$sample_metadata, subome = oi))
          
          md_meta <- unique(do.call(rbind, mds_meta)[,c("Timepoint", "vialLabel", "pid", "Sex", "BMI", "subome")])
          
        } else {
          md <- molecular_data[[tissue]][[ome]]$qc_norm
          md_meta <- molecular_data[[tissue]][[ome]]$sample_metadata
        }
        md_meta <- md_meta[md_meta$Timepoint == "pre_exercise", 
                           c("vialLabel", "pid", "Sex", "BMI", "subome")]
        vial2pid <- setNames(md_meta$pid, md_meta$vialLabel)
        # vial2subome <- setNames(md_meta$subome, md_meta$vialLabel)
        # md$ENSG <-  gsub("\\..*", "", rownames(md))
        md <- md[,colnames(md) %in% md_meta$vialLabel]
        colnames(md) <- vial2pid[colnames(md)]
        md <- compact_by_name(md)
        
        #which sub-ome did each analyte come from?
        # apply(md, 1, function(molrow) 
        #   unique(vial2subome[colnames(md)][which(!is.na(molrow))]))
        # intersect(rownames(molecular_data$muscle$`metab-u-ionpneg`$qc_norm),
        #           rownames(molecular_data$muscle$`metab-u-rpneg`$qc_norm))
        # nsubome <- names(analyte2subome)
        # intersecting_anal_names <- sapply(nsubome, function(xi){
        #   sapply(nsubome, function(yi){
        #     length(intersect(analyte2subome[[xi]], analyte2subome[[yi]]))})  })
        # intersecting_anal_names[upper.tri(intersecting_anal_names)] #ok phew none here!
        analyte2subome_df <- do.call(rbind, lapply(names(analyte2subome), function(a2si) 
          data.frame(mol = analyte2subome[[a2si]], subome = a2si)))
        analyte2subome_key <- setNames(nm = analyte2subome_df$mol, 
                                       object = analyte2subome_df$subome)
        a2s_keys[[ome]] <- analyte2subome_key
        
        #drop NA columns of md
        md <- md[,!apply(apply(md, 2, is.na), 2, any)]
        
        # if(!(ome %in% names(molecular_data[[tissue]]))){
        #   split_md <- split(data.frame(t(md)), rownames(t(md)))
        #   md_list <- lapply(split_md, function(smdi){
        #     smdi_rows <- apply(smdi, 1, identity, simplify = F)
        #     smdi_rows <- sapply(smdi_rows, function(smdir) smdir[!is.na(smdir)])
        #     values <- data.frame(unlist(smdi_rows, use.names = FALSE))
        #     rownames(values) <- unlist(sapply(smdi_rows, names))
        #     return(values)
        #   })
        #   n_indiv_per_n_analyte <- table(sapply(md_list, nrow))
        #   indiv_to_keep <- sapply(md_list, nrow) == as.numeric(names(n_indiv_per_n_analyte)[which.max(n_indiv_per_n_analyte)])
        #   md_list <- md_list[indiv_to_keep]
        #   md <- do.call(cbind, md_list)
        #   colnames(md) <- names(md_list)
        # }
        
        md_omes[[ome]] <- md
      }
      
      #consolidate molecular partitions
      shared_pids <- Reduce(intersect, lapply(lapply(md_omes, colnames), unique))
      md_omes <- lapply(md_omes, function(md_ome) md_ome[,shared_pids])
      
      #load phenotypic data
      pd <- as.data.frame(fread("~/data/motrpac_precovid_sedentary_clinical_obs.txt"))
      pd <- pd[pd$pid %in% shared_pids,]
      md_omes <- lapply(md_omes, function(md_ome) t(md_ome[,as.character(pd$pid)]))
      
      #pre-specify train-test split(s)
      nsamp <- nrow(pd)
      if(nfolds > 1){
        test_breaks <- round(seq(1, nsamp+1, length.out = nfolds + 1))
        dat_random_order <- sample(1:nsamp)
        train_sets <- list()
        test_sets <- list()
        for(fi in 1:nfolds){
          test_sets[[fi]] <- dat_random_order[test_breaks[fi]:(test_breaks[fi+1]-1)]
          train_sets[[fi]] <- setdiff(1:nsamp, test_sets[[fi]])  
        }
      } else {
        ntrain <- ceiling(nsamp * ptrain)
        train_sets <- list()
        train_sets[[1]] <- sample(1:nsamp, ntrain)
        test_sets <- list()
        test_sets[[1]] <- setdiff(1:nsamp, train_sets[[1]])
      }
      
      
      #initialize results objects
      cv_hyperparams <- list()
      nonzero_coefs <- list()
      pred_obs <- list()
      pred_obs_train <- list()
      pred_obs_om <- list()
      summs <- list()
      fits <- list()
      norms <- list()
      X_tests <- list()
      X_trains <- list()
      
      #### run analyses ####
      for(fi in 1:nfolds){
        
        print(paste0("cv-set ", fi))
        
        train_inds <- train_sets[[fi]]
        test_inds <- test_sets[[fi]]
        
        #recover outcome variable
        y_all <- y_all_raw <- pd[,target_phenotype]
        y_train_raw <- y_all_raw[train_inds]
        
        #apply data manipulation of interest to outcome
        if(unnorm_VO2max & target_phenotype == "vo2_peak_mlkg"){
          y_all <- y_all * pd$body_weight
        }
        if(log_transform){
          y_all <- log(y_all)
        }
        
        y_train_unnorm <- y_all[train_inds]
        y_mu <- mean(y_train_unnorm)
        y_sd <- sd(y_train_unnorm)
        y_train <- (y_train_unnorm - y_mu) / y_sd
        
        #identify phenotypic predictors
        basic_feat_train <- pd[train_inds, phenotype_predictors]
        basic_feat_train$sex <- (basic_feat_train$sex == "Male") * 1 - 1/2
        
        #and molecular predictors
        mol_feat_train_sub_list <- list()
        mol_PCs_train_list <- list()
        mol_kPCA_train_list <- list()
        molPCA_mu_list <- list()
        molPCA_sd_list <- list()
        for(md_name in names(md_omes)){
          md <- md_omes[[md_name]]
          mol_feat_train <- md[train_inds,]
          n_mol <- min(n_mol_target, ncol(mol_feat_train))
          n_mol_ISIS <- min(n_mol_target_ISIS, ncol(mol_feat_train))
          
          if(!do_ISIS_hclust){
            mol_feat_corrs <- apply(mol_feat_train, 2, function(xp) cor(xp, y_train))
            mol_inds <- order(abs(mol_feat_corrs), decreasing = T)[1:n_mol]
            mol_feat_train_sub <- mol_feat_train[,mol_inds]  
            
            #check correlations between top correlates of outcome?
            n2check <- 10
            summary(cor(mol_feat_train[,mol_inds[1:n2check]])[upper.tri(diag(n2check))])
            summary(ppcor::pcor(mol_feat_train[,mol_inds[1:n2check]])$estimate[upper.tri(diag(n2check))])
            ppcor::pcor(cbind(y_train, mol_feat_train[,mol_inds[1:n2check]]))$estimate[,1]
          } else {
            
            #ISIS screen
            tmp_out <- capture.output({
              sis_fit <- SIS::SIS(
                x      = mol_feat_train,
                y      = y_train,
                nsis   = n_mol_ISIS,
                iter   = TRUE
              )
            })
            mol_inds_ISIS  <- sis_fit$ix0
            # n_mol <- min(length(mol_inds_SIS), n_mol)
            
            mol_feat_train_sub_ISIS <- mol_feat_train[,mol_inds_ISIS, drop = F]
            ISIS_cor_y_train <- cor(mol_feat_train_sub_ISIS, y_train)
            ISIS_cor_y_train <- setNames(object = as.numeric(ISIS_cor_y_train), rownames(ISIS_cor_y_train))
            ISIS_cor <- cor(mol_feat_train_sub_ISIS)
            ISIS_dist <- 1-abs(ISIS_cor)
            
            #cut hclust and pick best feature in cluster
            hclust_ISIS <- hclust(as.dist(ISIS_dist))
            hclust_id <- cutree(hclust_ISIS, k = n_mol)
            mol_groups <- split(names(hclust_id), hclust_id)
            mol_inds <- sapply(mol_groups, function(mg) mg[which.max(abs(ISIS_cor_y_train[mg]))])
            mol_feat_train_sub <- mol_feat_train[,mol_inds]
            
          }
          
          #compute PC scores on molecular data
          #do we want to do "supervised PCA"?
          if(use_all_for_PCA){
            mol_PCA_input_train <- mol_feat_train
          } else {
            mol_PCA_input_train <- mol_feat_train_sub
          }
          n_kpc <- min(80L, nrow(mol_PCA_input_train) - 1L)
          
          if(standardize_before_PCA){
            molPCA_mu <- colMeans(mol_PCA_input_train)
            molPCA_sd <- apply(mol_PCA_input_train, 2, sd)
            mol_PCA_input_train <- t(t(mol_PCA_input_train) - molPCA_mu)
            mol_PCA_input_train <- t(t(mol_PCA_input_train) / molPCA_sd)
          }
          
          hyper_kPCA <- F
          if(hyper_kPCA){
            med_dist <- as.numeric(stats::median(dist(mol_PCA_input_train)))
            kpar     <- list(sigma = 1 / (2 * med_dist^2))
            mol_kPCA_train <- kernlab::kpca(
              mol_PCA_input_train,
              kernel   = "rbfdot",
              kpar     = kpar,
              features = n_kpc
            )
            # training scores (will have n_kpc columns)
            mol_PCs_train <- kernlab::predict(mol_kPCA_train, mol_PCA_input_train)
          } else {
            mol_kPCA_train <- kpca(mol_PCA_input_train, kernel = "rbfdot") 
            mol_PCs_train <- mol_kPCA_train@rotated[,1:n_kpc]  
          }
          
          #store in list
          molPCA_sd_list[[md_name]] <- molPCA_sd
          molPCA_mu_list[[md_name]] <- molPCA_mu
          mol_kPCA_train_list[[md_name]] <- mol_kPCA_train
          mol_feat_train_sub_list[[md_name]] <- mol_feat_train_sub
          mol_PCs_train_list[[md_name]] <- mol_PCs_train
          
        }
        
        #initialize training set list
        if(phenotypic_predictors){
          x_raw_train_list <- list(
            phenotypes        = basic_feat_train
          )  
        } else {
          x_raw_train_list <- list()
        }
        
        #add in molecular features
        ma_names <- paste0("molecular_analyte.", names(mol_feat_train_sub_list))
        if(molecular_predictors){
          x_raw_train_list <- c(x_raw_train_list, setNames(mol_feat_train_sub_list, ma_names))  
        }

        if(incl_latent_variables){
          x_raw_train_list <- c(x_raw_train_list, 
                                setNames(mol_PCs_train_list,
                                         paste0("molecular_kPCs.", names(mol_PCs_train_list)))
          )
        }
        
        if(residualize_molecular_feats & molecular_predictors){
          resid_mol_train_list <- lapply(mol_feat_train_sub_list, function(mol_anal){
            residualize_feats(output = mol_anal, input = basic_feat_train)
          })
          resid_mol_feat_train_list <- lapply(resid_mol_train_list, function(rml) rml$R)
          x_raw_train_list[ma_names] <- resid_mol_feat_train_list
        }
        
        if(sex_interactions){
          mol_feat_train_sex_inter_list <- lapply(x_raw_train_list[ma_names], function(mol_anal){
            mol_anal * basic_feat_train$sex
          })
          max_names <- gsub("molecular_analyte", "molecular_analyte_x_sex", names(mol_feat_train_sex_inter_list))
          names(mol_feat_train_sex_inter_list) <- max_names
          basic_feat_train_sex_inter <- basic_feat_train[,-which(colnames(basic_feat_train) == "sex")] * 
            basic_feat_train$sex
          x_raw_train_list[c("phenotypes_x_sex", max_names)] <- c(list(
            phenotypes_x_sex        = basic_feat_train_sex_inter),
            ifelse2(molecular_predictors, mol_feat_train_sex_inter_list, NULL))
        }
        
        # 2) compute per-block normalization stats on TRAIN
        x_mu_list <- lapply(x_raw_train_list, colMeans)
        x_sd_list <- lapply(x_raw_train_list, function(M) apply(M, 2, sd))
        
        # 3) z-score on TRAIN
        x_train_list <- mapply(zs, x_raw_train_list, x_mu_list, x_sd_list, SIMPLIFY = FALSE)
        
        # 4) stable PC names (only for that block if present)
        kpc_names_in_list <- names(x_train_list)[grep("^molecular_kPCs\\.", names(x_train_list))]
        if (length(kpc_names_in_list) > 0){
          # Loop and assign column names from the corresponding training block
          for(kpc_name in kpc_names_in_list){
            if(kpc_name %in% names(x_train_list)){
              colnames(x_train_list[[kpc_name]]) <- colnames(x_train_list[[kpc_name]])
            }
          }
        }
        
        # 5) record normalization constants (and transforms) for this fold
        norms[[fi]] <- list(
          x_mu_list = x_mu_list,
          x_sd_list = x_sd_list
        )
        
        #examine marginal associations?
        # apply(x1_train, 2, function(xp) cor(xp, y_train))
        # hist(apply(x2_train, 2, function(xp) cor(xp, y_train)))
        # hist(apply(x3_train, 2, function(xp) cor(xp, y_train)))
        
        #get some high-level metadata
        X_train <- do.call(cbind, x_train_list)
        p_vec   <- vapply(x_train_list, ncol, integer(1L))
        names(p_vec) <- names(x_train_list)
        G <- length(p_vec)
        block <- rep(names(p_vec), times = p_vec)
        block_id <- rep(1:G, times = p_vec)
        n <- nrow(X_train)
        X_trains[[fi]] <- X_train
        
        # prep inputs for test sets
        y_test_unnorm <- y_all[test_inds]
        basic_feat_test <- pd[test_inds, phenotype_predictors]
        basic_feat_test$sex <- (basic_feat_test$sex == "Male") * 1 - 1/2
        
        #molecular features for test set
        mol_feat_test_sub_list <- list()
        mol_PCs_test_list <- list()
        for(md_name in names(md_omes)){
          md <- md_omes[[md_name]]
          mol_feat_test <- md[test_inds,, drop = FALSE]
          
          # Get column indices from the saved training sub-features
          mol_inds_names <- colnames(mol_feat_train_sub_list[[md_name]])
          mol_feat_test_sub <- mol_feat_test[, mol_inds_names, drop = FALSE]
          
          #now project to the earlier PC axes
          if (use_all_for_PCA) {
            mol_PCA_input_test <- mol_feat_test
          } else {
            mol_PCA_input_test <- mol_feat_test_sub
          }
          
          if (standardize_before_PCA) {
            # Use the corresponding mean/sd from the training list
            mol_PCA_input_test <- t( t(mol_PCA_input_test) - molPCA_mu_list[[md_name]] )
            mol_PCA_input_test <- t( t(mol_PCA_input_test) / molPCA_sd_list[[md_name]] )
          }
          
          # embed test using the corresponding kPCA model
          mol_kPCA_train <- mol_kPCA_train_list[[md_name]]
          # Get n_kpc from the *saved training PCs* for this 'ome'
          n_kpc_current <- ncol(mol_PCs_train_list[[md_name]]) 
          
          mol_PCs_test_full <- kernlab::predict(mol_kPCA_train, mol_PCA_input_test)
          
          if (ncol(mol_PCs_test_full) < n_kpc_current) {
            stop(paste0("Test KPCA for ome '", md_name, "' returned fewer PCs than training."))
          }
          
          mol_PCs_test <- mol_PCs_test_full[, seq_len(n_kpc_current), drop = FALSE]
          
          #store in list
          mol_feat_test_sub_list[[md_name]] <- mol_feat_test_sub
          mol_PCs_test_list[[md_name]] <- mol_PCs_test
        }
        
        #initialize test set list
        if(phenotypic_predictors){
          x_raw_test_list <- list(
            phenotypes = basic_feat_test
          )
        } else {
          x_raw_test_list <- list()
        }
        
        #add in molecular features
        ma_names <- paste0("molecular_analyte.", names(mol_feat_test_sub_list))
        if(molecular_predictors){
          x_raw_test_list <- c(x_raw_test_list, setNames(mol_feat_test_sub_list, ma_names))  
        }

        if(incl_latent_variables){
          kpc_names <- paste0("molecular_kPCs.", names(mol_PCs_test_list))
          x_raw_test_list <- c(x_raw_test_list, setNames(mol_PCs_test_list, kpc_names))
        }
        
        if(molecular_predictors & residualize_molecular_feats){
          # Get the names of the molecular blocks (e.g., "molecular_analyte.ome1", etc.)
          ma_names_in_list <- names(x_raw_test_list)[grep("^molecular_analyte\\.", names(x_raw_test_list))]
          
          for(ma_name in ma_names_in_list){
            # Get the original 'ome' name (e.g., "ome1")
            md_name <- gsub("molecular_analyte\\.", "", ma_name)
            
            # Get the model (B) from the corresponding training list item
            B_matrix <- resid_mol_train_list[[md_name]]$B 
            
            # Apply projection: Y_test = Y_test - X_test %*% B
            x_raw_test_list[[ma_name]] <- x_raw_test_list[[ma_name]] - 
              cbind(1, as.matrix(basic_feat_test)) %*% B_matrix
          }
        }
        
        if(sex_interactions){
          ma_names_in_list <- names(x_raw_test_list)[grep("^molecular_analyte\\.", names(x_raw_test_list))]
          mol_feat_test_sex_inter_list <- lapply(x_raw_test_list[ma_names_in_list], function(mol_anal){
            mol_anal * basic_feat_test$sex
          })
          max_names <- gsub("molecular_analyte", "molecular_analyte_x_sex", names(mol_feat_test_sex_inter_list))
          names(mol_feat_test_sex_inter_list) <- max_names
          basic_feat_test_sex_inter <- basic_feat_test[,-which(colnames(basic_feat_test) == "sex")] *
            basic_feat_test$sex
          x_raw_test_list <- c(x_raw_test_list, 
                               list(phenotypes_x_sex = basic_feat_test_sex_inter),
                               ifelse2(molecular_predictors, NULL, mol_feat_test_sex_inter_list))
        }
        
        # use TRAIN means/sds
        x_test_list <- Map(zs, x_raw_test_list, norms[[fi]]$x_mu_list, norms[[fi]]$x_sd_list)
        
        # align colnames for kPCs
        kpc_names_in_list <- names(x_test_list)[grep("^molecular_kPCs\\.", names(x_test_list))]
        if (length(kpc_names_in_list) > 0){
          # Loop and assign column names from the corresponding training block
          for(kpc_name in kpc_names_in_list){
            if(kpc_name %in% names(x_train_list)){
              colnames(x_test_list[[kpc_name]]) <- colnames(x_train_list[[kpc_name]])
            }
          }
        }
        
        #compile total data frame
        X_test <- do.call(cbind, x_test_list)
        X_tests[[fi]] <- X_test
        
        # make X sparse to save mem + speed
        # if (!inherits(X_train, "dgCMatrix")) X <- Matrix(X_train, sparse = TRUE)
        
        #### lasso ####
        if(method == "lasso"){
          
          # options
          family <- "gaussian"
          alpha <- 1
          nfolds_internal <- 10
          maxit <- 1e5
          standardize <- FALSE   # set TRUE if columns aren't already on comparable scales
          
          # anchor block weights by sqrt(size)
          w_anchor <- sqrt(p_vec / p_vec[1])
          
          # coarse multiplicative grid around anchors (you already had a nice grid; both are fine)
          f_grid <- 2^(-6:3)
          grid_list <- setNames(rep(list(f_grid), G), paste0("f", seq_len(G)))
          combos <- do.call(expand.grid, c(grid_list, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
          max_grid <- 1024
          if (nrow(combos) > max_grid) {
            combos <- combos[sample.int(nrow(combos), max_grid), , drop = FALSE]
          }
          make_pf <- function(i) {
            mult <- as.numeric(combos[i, paste0("f", seq_len(G)), drop = TRUE])  # f1..fG
            w <- w_anchor * mult
            pf <- rep_by_group(w, p_vec)
            pf / mean(pf)
          }
          foldid <- sample(rep_len(1:nfolds_internal, n))
          
          # cores (avoid saturating everything; leave 1–2 cores free if you like)
          mc_cores <- min(detectCores(), 8)
          cat("# running", nrow(combos), "weight combos on", mc_cores, "cores...\n")
          
          # run grid in parallel; each job returns a small list
          results <- mclapply(
            X = seq_len(nrow(combos)),
            FUN = function(i) {
              pf <- make_pf(i)
              
              if(i %%100==0) mcprint(paste0("(", i, " / ", nrow(combos), ") "))
              mult <- as.numeric(combos[i, paste0("f", seq_len(G)), drop = TRUE])
              w_used <- w_anchor * mult
              cvfit <- cv.glmnet(X_train, y_train,
                                 family = family,
                                 nfolds = nfolds_internal,
                                 foldid = foldid,
                                 standardize = standardize,
                                 alpha = alpha,
                                 penalty.factor = pf,
                                 intercept = TRUE,
                                 maxit = maxit,
                                 parallel = FALSE)  # do NOT parallelize inside
              
              list(i = i,
                   w = w_used,
                   cvmin = min(cvfit$cvm),
                   lambda_min = cvfit$lambda.min)
            },
            mc.cores = mc_cores,
            mc.set.seed = TRUE   # independent RNG streams per job
          )
          
          # collect and pick the best
          tab <- do.call(rbind, lapply(results, function(z) {
            cbind(i = z$i,
                  {
                    foo <- t(as.matrix(z$w))
                    colnames(foo) <- paste0("w", seq_len(G))
                    foo
                  },
                  cvmin = z$cvmin,
                  lambda_min = z$lambda_min)
          }))
          tab <- as.data.frame(tab)
          best_row <- tab[which.min(tab$cvmin), ]
          best_lambda <- best_row$lambda_min
          cat("# best weights = ",
              paste(paste0(names(p_vec), "=", 
                           signif(as.numeric(best_row[paste0("w", seq_len(G))]), 3)),
                    collapse = ", "),
              "  lambda.min = ", signif(best_lambda, 4), "\n", sep = "")
          
          # refit final model at best weights and lambda.min
          pf_best <- (function() {
            w_best <- as.numeric(best_row[paste0("w", seq_len(G))])
            pf <- rep_by_group(w_best, p_vec)
            pf / mean(pf)
          })()
          
          final_fit <- glmnet(X_train, y_train,
                              family = family,
                              standardize = standardize,
                              alpha = alpha,
                              penalty.factor = pf_best,
                              intercept = TRUE,
                              lambda = best_lambda,
                              maxit = maxit)
          cm       <- coef(final_fit, s = best_lambda)
          beta_hat <- as.numeric(cm)[-1]
          a_hat    <- as.numeric(cm)[1]
          nz_mask  <- as.numeric(cm != 0)[-1] == 1
          names(beta_hat) <- colnames(X_train)
          fits[[fi]] <- final_fit
          
          #### uni-lasso ####
        } else if (method == "unilasso") {
          
          # options
          family <- "gaussian"
          nfolds_internal <- 10
          maxit <- 1e5
          standardize <- FALSE        # you already standardized X upstream
          foldid <- sample(rep_len(1:nfolds_internal, n))
          
          # === Optional block reweighting via feature scaling ===
          use_block_rescale <- FALSE  # set TRUE if you want differential penalties
          ul_block_mult <- setNames(rep(1, length(p_vec)), names(p_vec))  # named by block
          
          if (use_block_rescale) {
            # choose per-block multipliers (smaller -> less penalty effectively)
            pf_target <- rep_by_group(ul_block_mult, p_vec)
            # emulate penalty.factor = pf_target by scaling columns by s_k = 1/pf_k
            s_vec <- 1 / pf_target
            X_train_ul <- sweep(X_train, 2, s_vec, `*`)
            X_test_ul  <- sweep(X_test,  2, s_vec, `*`)
          } else {
            X_train_ul <- X_train
            X_test_ul  <- X_test
            s_vec      <- rep(1, ncol(X_train))  # identity back-map
          }
          
          # === Package uniLasso with CV (no penalty.factor here) ===
          final_fit <- uniLasso::cv.uniLasso(
            x = X_train_ul, y = y_train,
            family = family,
            standardize = standardize,
            nfolds = nfolds_internal,
            foldid = foldid,
            lower.limits = 0,   # stage-2 nonnegative (default)
            maxit = maxit
          )
          fits[[fi]] <- final_fit
          
          best_lambda <- final_fit$lambda.min
          
          # Coefs (on scaled design if use_block_rescale == TRUE)
          cm <- coef(final_fit, s = best_lambda)   # intercept + gamma
          gamma_hat_scaled <- as.numeric(cm)[-1]
          a_hat            <- as.numeric(cm)[1]
          
          # Back-map to original X scale if you scaled the design:
          beta_hat <- gamma_hat_scaled * s_vec
          names(beta_hat) <- colnames(X_train)
          
          # exact nonzeros mask for reporting
          nz_mask <- as.numeric(cm != 0)[-1] == 1

          #### horseshoe ####
        } else if(method == "horseshoe"){
          
          # prepare inputs
          blk_names <- names(p_vec)
          
          # Heuristic default for expected nonzeros per block (edit if you have priors):
          p0_guess_base <- c(phenotypes = 5, molecular_analyte = 50, molecular_kPCs = 20)
          p0_guess_base_matches <- sapply(names(p0_guess_base), function(p0gb) which(grepl(p0gb, blk_names)), simplify = F)
          p0_guess <- rep(NA, length(blk_names))
          p0_guess[unlist(p0_guess_base_matches)] <- p0_guess_base[rep(names(p0_guess_base_matches), sapply(p0_guess_base_matches, length))]
          p0_guess[is.na(p0_guess)] <- pmax(1L, round(0.05 * p_vec[is.na(p0_guess)]))
          names(p0_guess) <- blk_names
          
          # slab scales: keep your old “phenotypes=4, others=2” by name heuristic, else 2
          slab_scale <- if (length(blk_names)) {
            sapply(blk_names, function(nm) if (grepl("pheno", nm, ignore.case = TRUE)) 4 else 2)
          } else {
            rep(2, G)
          }
          slab_df <- 5
          
          # compute tau0 by block
          tau0_block <- function(p0, d, n, sigma = 1) (p0 / (d - p0)) * (sigma / sqrt(n))
          hs_tau0 <- mapply(function(p0, d) tau0_block(p0, d, n),
                            p0 = p0_guess, d = p_vec, SIMPLIFY = TRUE)
          hs <- list(
            hs_tau0    = hs_tau0,
            slab_scale = slab_scale,
            slab_df    = slab_df
          )
          
          #compute importance score annotation for HS local scales
          #first simulate fake annotation
          n_annot <- 3
          mol_annots <- which(block_id == which(grepl("molecular_analyte", blk_names)))
          annot <- matrix(0, ncol(X_train), n_annot) 
          annot[mol_annots,] <- rnorm(length(mol_annots) * n_annot)
          
          # pack data for Stan
          dat <- list(
            n = n,
            d = ncol(X_train),
            x = X_train,
            y = as.vector(y_train),
            G = as.integer(G),
            group = as.integer(block_id),
            hs_tau0 = as.vector(hs_tau0),
            slab_scale = as.vector(slab_scale),
            slab_df = slab_df,
            K = n_annot,
            annot = annot
          )
          
          # fit model
          if(use_annotation){
            hs_mod <- cmdstanr::cmdstan_model("~/repos/motrpac-predict/src/Stan/multi-horseshoe-informative.stan")  
          } else {
            hs_mod <- cmdstanr::cmdstan_model("~/scripts/minor_scripts/postdoc/multi-horseshoe.stan")  
          }
          
          
          if(use_pf){
            fit <- hs_mod$pathfinder(
              data = dat, draws = 2000, num_paths = 8,
              refresh = ceiling(mcmc_params[["n_iter"]] / 20)
            )
          } else {
            fit <- hs_mod$sample(
              data = dat,
              chains = 4,
              iter_warmup = mcmc_params[["n_iter"]],
              iter_sampling = mcmc_params[["n_iter"]],
              thin = 1,
              adapt_delta = mcmc_params[["adapt_delta"]],
              max_treedepth = mcmc_params[["max_treedepth"]],
              parallel_chains = 4,
              refresh = ceiling(mcmc_params[["n_iter"]] / 20)
            )
            
            # check diagnostics
            summ <- fit$summary()
            summs[[fi]] <- summ
            summ_cols <- c("variable", "mean", "rhat", "ess_bulk", "ess_tail")
            print(summ[order(summ$rhat, decreasing = TRUE), summ_cols])
            print(summ[order(summ$ess_bulk), summ_cols])
            
          }
          fits[[fi]] <- fit
          
          # extract parameters
          beta_draws <- posterior::as_draws_matrix(fit$draws("beta"))
          a_draws    <- as.numeric(posterior::as_draws_matrix(fit$draws("a")))
          sigma_draws<- as.numeric(posterior::as_draws_matrix(fit$draws("sigma")))
          
          if (hs_estimator == "median") {
            beta_hat <- apply(beta_draws, 2, median)
            a_hat    <- median(a_draws)
          } else {
            beta_hat <- colMeans(beta_draws)
            a_hat    <- mean(a_draws)
          }
          
          if (hs_pred_type == "plugin") {
            yhat_test_z <- as.numeric(a_hat + X_test %*% beta_hat)
            yhat_test   <- yhat_test_z * y_sd + y_mu
          } else { # posterior predictive mean
            m <- min(1000L, nrow(beta_draws))
            idx <- sample(seq_len(nrow(beta_draws)), m)
            B   <- beta_draws[idx, , drop = FALSE]          # m x d
            A   <- a_draws[idx]
            # n_test x m matrix of z-scale predictive means
            Yhat_z <- tcrossprod(X_test, B)
            Yhat_z <- sweep(Yhat_z, 2L, A, "+")
            yhat_test <- rowMeans(Yhat_z) * y_sd + y_mu     # mean-only back-transform
          }
          
        }
        
        #### make predictions ####
        if (method == "horseshoe") {
          # Thresholding path for HS (since it’s not sparse at exactly zero)
          nz <- which(abs(beta_hat) > hs_beta_thresh)
        } else if (method == "unilasso") {
          # exact nonzeros from coef() we saved as nz_mask
          nz <- which(nz_mask)
        } else { # lasso
          nz_mask <- as.numeric(cm != 0)[-1] == 1
          nz <- which(nz_mask)
        }
        
        #print summary of # of nz coefs across blocks
        cat(paste(c(paste0("# nonzeros total = ", 
                           (if (is.logical(nz)) sum(nz) else length(nz))), 
                    paste(u <- unique(x <- block[if (is.logical(nz)) which(nz) else nz]), 
                          tabulate(match(x, u), length(u)), sep = " = ")), collapse = " | "), "\n")
        
        nz_coefs <- split(setNames(beta_hat[nz], colnames(X_train)[nz]), block[nz])
        
        #try doing what that nature med paper did?
        only_molecular <- c(T, F)
        X_test_orig <- X_test
        y_test_unnorm_orig  <- y_test_unnorm
        y_train_unnorm_orig <- y_train_unnorm
        X_train_orig <- X_train
        y_train_unnorm_orig  <- y_train_unnorm
        y_train_unnorm_orig <- y_train_unnorm
        for(om in only_molecular){
          
          #we modify these so they need to be reset
          X_test <- X_test_orig
          y_test_unnorm  <- y_test_unnorm_orig
          y_train_unnorm <- y_train_unnorm_orig
          X_train <- X_train_orig
          y_train_unnorm  <- y_train_unnorm_orig
          y_train_unnorm <- y_train_unnorm_orig
          
          if(om){
            set_to_test_mean <- F
            if(set_to_test_mean){
              #set phenotypic variables to their mean value
              X_test[,grepl("phenotype", block)] <-
                t(t(X_test[,grepl("phenotype", block)]) -
                    t(X_test[,grepl("phenotype", block)]) + 
                    colMeans(X_test[,grepl("phenotype", block)]))
              X_train[,grepl("phenotype", block)] <-
                t(t(X_train[,grepl("phenotype", block)]) -
                    t(X_train[,grepl("phenotype", block)]) + 
                    colMeans(X_train[,grepl("phenotype", block)])) 
            } else {
              #set phenotypic variables to their mean value
              X_test[,grepl("phenotype", block)] <- 0
              X_train[,grepl("phenotype", block)] <- 0
            } 
          }
          
          # 3D) predictions on test (z-scale)
          if (method == "horseshoe") {
            yhat_test_z <- as.numeric(a_hat + X_test %*% beta_hat)
            yhat_train_z <- as.numeric(a_hat + X_train %*% beta_hat)
          } else if (method == "lasso") {
            yhat_test_z  <- as.numeric(predict(final_fit, newx = X_test,  s = best_lambda))
            yhat_train_z <- as.numeric(predict(final_fit, newx = X_train, s = best_lambda))
          } else if (method == "unilasso"){
            yhat_test_z  <- as.numeric(predict(final_fit,
                                               newx = if (use_block_rescale) X_test_ul else X_test,
                                               s = best_lambda))
            yhat_train_z <- as.numeric(predict(final_fit,
                                               newx = if (use_block_rescale) X_train_ul else X_train,
                                               s = best_lambda))
          }
          
          #debugging check... we are getting positive train r2 on the z-scale right
          sst_z <- sum( (y_train - mean(y_train))^2 )
          sse_z <- sum( (y_train - yhat_train_z)^2 )
          r2_train_fold_z <- 1 - sse_z/sst_z
          cat(paste0(ifelse(om, "om ", ""), "R2 train z = ", r2_train_fold_z, "\n"))
          
          #post-hoc recalibration on train to adjust for mean-regression
          if(linear_recalibration){
            cal_out <- calibrate_linear(y_train_z = y_train, 
                                        yhat_train_z = yhat_train_z, 
                                        yhat_test_z = yhat_test_z)
            yhat_train_z <- cal_out$yhat_train_cal_z
            yhat_test_z <- cal_out$yhat_test_cal_z
          }
          
          #undo centering and scaling operations
          yhat_test  <- yhat_test_z  * y_sd + y_mu
          yhat_train <- yhat_train_z * y_sd + y_mu
          
          #exponentiate yhat if we log-transformed earlier
          if(log_transform){ #need to use smearing backtransform
            resid_log <- y_train_unnorm - yhat_train
            smear_factor <- mean(exp(resid_log), na.rm = TRUE)
            
            yhat_test <- exp(yhat_test) * smear_factor
            y_test_unnorm <- exp(y_test_unnorm)
            
            yhat_train <- exp(yhat_train) * smear_factor
            y_train_unnorm <- exp(y_train_unnorm)
          }
          
          #get vo2max back on the per kg scale
          if(unnorm_VO2max & target_phenotype == "vo2_peak_mlkg"){
            yhat_test <- yhat_test / pd$body_weight[test_inds]
            y_test_unnorm <- y_test_unnorm / pd$body_weight[test_inds]
            
            yhat_train <- yhat_train / pd$body_weight[train_inds]
            y_train_unnorm <- y_train_unnorm / pd$body_weight[train_inds]
          }
          
          # compute test R^2
          sst <- sum( (y_test_unnorm - mean(y_test_unnorm))^2 )
          sse <- sum( (y_test_unnorm - yhat_test)^2 )
          r2_test <- 1 - sse/sst
          cat(paste0(ifelse(om, "om ", ""), sprintf("test R^2 = %.4f\n", r2_test)))
          
          # cv_hyperparams row: keep structure; fill NAs for unused fields
          if (method == "lasso") {
            cv_hyperparams[[fi]] <- best_row
          } else if (method == "unilasso") {
            cv_hyperparams[[fi]] <- data.frame(
              i = NA,
              lambda_min = best_lambda,
              method = "unilasso",
              stringsAsFactors = FALSE
            )
          } else { # horseshoe
            # create a row with the same columns as best_row if possible; else a superset is fine
            hs_row <- data.frame(
              i = NA, w1 = NA, w2 = NA, w3 = NA,
              cvmin = NA, lambda_min = NA,
              hs_tau01 = hs$hs_tau0[1], hs_tau02 = hs$hs_tau0[2], hs_tau03 = hs$hs_tau0[3],
              slab1 = hs$slab_scale[1], slab2 = hs$slab_scale[2], slab3 = hs$slab_scale[3],
              slab_df = hs$slab_df,
              method = "horseshoe",
              stringsAsFactors = FALSE
            )
            cv_hyperparams[[fi]] <- hs_row
          }
          
          #save the full predictions to the standard lists
          #and the molecular only predictions to a special list
          if(!om){
            nonzero_coefs[[fi]] <- nz_coefs
            pred_obs[[fi]] <- data.frame(obs = y_test_unnorm, pred = yhat_test,
                                         pid = pd$pid[test_inds])
            pred_obs_train[[fi]] <- data.frame(obs = y_train_unnorm, pred = yhat_train,
                                               pid = pd$pid[train_inds])  
          } else {
            pred_obs_om[[fi]] <- data.frame(obs = y_test_unnorm, pred = yhat_test,
                                         pid = pd$pid[test_inds]) 
          }
        }
        X_test <- X_test_orig #put the original X_test back
        X_train <- X_train_orig #put the original X_train back
        
      }
      
      #### aggregate results ####
      all_pred_obs <- do.call(rbind, pred_obs)[,c("obs", "pred")]
      all_cv_hyperparams <- do.call(rbind, cv_hyperparams)
      if(is.null(unlist(nonzero_coefs))){
        all_nonzero_coefs <- NULL
      } else {
        all_nonzero_coefs <- build_sparse_mats(nonzero_coefs)
      }
      rss <- sum((all_pred_obs$obs - all_pred_obs$pred)^2)
      tss <- sum((all_pred_obs$obs - mean(all_pred_obs$obs))^2)
      r2 <- 1 - rss/tss
      
      #check predictions on training data too
      all_pred_obs_train_list <- do.call(rbind, pred_obs_train)
      all_pred_obs_train_list <- split(all_pred_obs_train_list, all_pred_obs_train_list$pid)
      all_pred_obs_train <- do.call(rbind, lapply(all_pred_obs_train_list, colSums))[,c("obs", "pred")] / (nfolds-1)
      all_pred_obs_train <- as.data.frame(all_pred_obs_train)
      rss_train <- sum((all_pred_obs_train$obs - all_pred_obs_train$pred)^2)
      tss_train <- sum((all_pred_obs_train$obs - mean(all_pred_obs_train$obs))^2)
      r2_train <- 1 - rss_train/tss_train
      
      #check predictions on training data too
      all_pred_obs_om <- do.call(rbind, pred_obs_om)[,c("obs", "pred")]
      rss_om <- sum((all_pred_obs_om$obs - all_pred_obs_om$pred)^2)
      tss_om <- sum((all_pred_obs_om$obs - mean(all_pred_obs_om$obs))^2)
      r2_om <- 1 - rss_om/tss_om
      
      #### save results ####
      results_subdir <- paste0(results_dir,
                               target_phenotype, "/", 
                               tissue, "/", 
                               ome, "/", 
                               method, "/")
      if(!dir.exists(results_subdir)){dir.create(results_subdir, recursive = T)}
      res_path <- paste0(results_subdir, "res",
                         ifelse(sex_interactions, "_sex-interaction", ""),
                         ".RData")
      res <- list(
        fits = fits,
        all_pred_obs = all_pred_obs,
        r2 = r2,
        cv_hyperparams = cv_hyperparams,
        nonzero_coefs = nonzero_coefs,
        pred_obs = pred_obs,
        all_pred_obs_train = all_pred_obs_train,
        summs = summs,
        fits = fits,
        target_phenotype = target_phenotype,
        tissue = tissue,
        ome = ome,
        method = method,
        test_sets = test_sets,
        norms = norms,
        phenotype_predictors = phenotype_predictors
      )
      save(res, file = res_path)
      
      #need to save predictions in both the train set 
      #(for eg meta-modeling) and test sets
      #and might as well save fitted models too
      
      
      #### meta-model ####
      
      #### plot ####
      npanel <- G+4
      fig_path <- paste0(subfigures_dir, "clinical-x-omics-prediction_", 
                         paste0(unlist(omes), collapse = ";"), "-", tissue, "-", target_phenotype, "-", method,
                         ifelse(sex_interactions, "_sex-interaction", ""), 
                         ifelse(residualize_molecular_feats, "_ortho-mol", ""), 
                         ".pdf")
      fig_paths <- c(fig_paths, fig_path)
      cairo_pdf(fig_path, width = 600/72, height = 500/72 * npanel / 4)
      par(mfrow = c(ceiling(npanel / 2),2), mar = c(7,5,6,5))
      if(G == 3){
        par(oma = c(1,1,1,1))
      } else {
        par(oma = c(3,1.5,1.5,1.5))
      }
      
      #plot validation set predictions
      if(length(pred_obs) <= 8){
        base_cols <- RColorBrewer::brewer.pal(length(pred_obs), "Dark2")  
      } else {
        base_cols <- rainbow(length(pred_obs))
      }
      
      #plot validation set predictions
      dotcols <- rep(adjustcolor(base_cols, 0.8), sapply(pred_obs, nrow))
      plot(all_pred_obs$obs, all_pred_obs$pred,
           main = paste0(ome, ", ", tissue, ", ", target_phenotype, ", ", method),
           xlab = paste0("observed"),
           ylab = paste0("predicted"), pch = 19,
           col = dotcols)
      title(main = bquote(R^2 %~~% .(paste0(round(r2, 3), " (over ", nfolds, " folds)"))), line = 0.75, 
            cex.main = 1)
      abline(0, 1, col = 2, lty = 2)
      legend("bottomright",
             legend = c("dot cols indicate fold", "1-to-1 line"),
             lty    = c(NA, 2),
             pch    = c(19, NA),
             col    = c("black", 2),
             bty    = "n")
      
      #plot training set prediction averages
      plot(all_pred_obs_train$obs, all_pred_obs_train$pred,
           main = paste0(ome, ", ", tissue, ", ", target_phenotype, ", ", method),
           xlab = paste0("observed"),
           ylab = paste0("predicted"), pch = 19,
           col = dotcols)
      title(main = bquote(R^2 %~~% .(paste0(round(r2_train, 3), " (over ", nfolds, " folds), TRAINING SET MEAN"))), line = 0.75, 
            cex.main = 1)
      abline(0, 1, col = 2, lty = 2)
      legend("bottomright",
             legend = c("dot cols indicate fold", "1-to-1 line"),
             lty    = c(NA, 2),
             pch    = c(19, NA),
             col    = c("black", 2),
             bty    = "n")
      
      
      #plot training set prediction averages
      plot(all_pred_obs_om$obs, all_pred_obs_om$pred,
           main = paste0(ome, ", ", tissue, ", ", target_phenotype, ", ", method),
           xlab = paste0("observed"),
           ylab = paste0("predicted"), pch = 19,
           col = dotcols)
      title(main = bquote(R^2 %~~% .(paste0(round(r2_om, 3), " (over ", nfolds, " folds), MOLECULAR ANALYTES ONLY"))), line = 0.75, 
            cex.main = 1)
      abline(0, 1, col = 2, lty = 2)
      legend("bottomright",
             legend = c("dot cols indicate fold", "1-to-1 line"),
             lty    = c(NA, 2),
             pch    = c(19, NA),
             col    = c("black", 2),
             bty    = "n")
      
      #plot estimated coefficients
      order_by_nsig <- T
      for(group in unique(block)){
        
        if(is.null(all_nonzero_coefs[[group]])){
          plot.new()
          plot.window(0:1, 0:1)
          text(0.5, 0.5, labels = paste0("no ", group, "\nwith non-zero coefs"))
          next()
        }
        
        nzc_sub <- as.matrix(all_nonzero_coefs[[group]])
        nzc_sub_ncoef <- apply(nzc_sub, 2, function(nzc) sum(abs(nzc) > 1E-6))
        nzc_sub_emag <- apply(nzc_sub, 2, function(nzc) mean(abs(nzc)[abs(nzc) > 1E-6]))
        
        min_ncoef <- 2
        nzc_sub <- nzc_sub[,nzc_sub_ncoef >= min_ncoef,drop=FALSE]
        nzc_sub_emag <- nzc_sub_emag[nzc_sub_ncoef >= min_ncoef]
        nzc_sub_ncoef <- nzc_sub_ncoef[nzc_sub_ncoef >= min_ncoef]
        if(order_by_nsig){
          nzcs_ord <- order(nzc_sub_ncoef, decreasing = T)  
        } else {
          nzcs_ord <- order(nzc_sub_emag, decreasing = T)
        }
        
        #reorder subset
        nzc_sub <- nzc_sub[,nzcs_ord,drop=FALSE]
        nzc_sub_emag <- nzc_sub_emag[nzcs_ord]
        nzc_sub_ncoef <- nzc_sub_ncoef[nzcs_ord]
        
        #start plotting
        if(length(nzcs_ord) == 0){
          plot.new()
          plot.window(0:1, 0:1)
          text(0.5, 0.5, labels = paste0("no ", group, "\nwith non-zero coefs"))
          next()
        }
        
        plot(nzc_sub_ncoef, type = "l", 
             ylab =  paste0("# non-zero estimates (/ ", nfolds, " folds)"),
             xaxt = "n", yaxt = "n", xlab = "", lwd = 3,
             main = paste0(tissue, ", ", ome, ", ", group))
        axis(2, at = min(nzc_sub_ncoef):max(nzc_sub_ncoef), 
             labels = min(nzc_sub_ncoef):max(nzc_sub_ncoef), las = 1)
        n_feats <- ncol(nzc_sub)
        abline(v = 1:n_feats, lty = 3, lwd = 0.75)
        
        pusr <- par("usr")
        segments(x0 = 1:n_feats, x1 = 1:n_feats,
                 y0 = pusr[3], y1 = pusr[3] - diff(pusr[3:4])/20, xpd = NA)
        labcex <- min(1, (diff(pusr[1:2]) / length(nzc_sub_ncoef)) / (strwidth("a") * 2.25))
        text(x = 1:n_feats + strwidth("a"), 
             y = pusr[3] - diff(pusr[3:4])/15, 
             labels = names(nzc_sub_ncoef), 
             xpd = NA,
             srt = 45, pos = 2, cex = labcex)
        if(grepl("pcs", tolower(group))){
          cols <- adjustcolor(c("grey10", "grey10"), 0.75)
        } else {
          cols <- adjustcolor(c("blue", "red"), 0.75)  
        }
        
        effects_by_analyte <- do.call(rbind, lapply(1:n_feats, function(ai){
          nz_vals <- nzc_sub[abs(nzc_sub[,ai]) > 1E-6, ai]
          data.frame(abs_val = abs(nz_vals), col = cols[(sign(nz_vals) == 1) + 1],
                     index = ai)
        }))
        range_effects <- c(0, max(effects_by_analyte$abs_val))
        range_yax <- pusr[3:4] + c(0.1, -0.1)
        effects_by_analyte$vloc <- effects_by_analyte$abs_val / 
          diff(range_effects) *  diff(range_yax) + range_yax[1]
        ytick_vals <- pretty(range_effects)
        ytick_vals <- ytick_vals[ytick_vals >= range_effects[1] &
                                   ytick_vals <= range_effects[2]
        ]
        ytick_locs <-  ytick_vals / diff(range_effects) *  diff(range_yax) + range_yax[1]
        points(x = effects_by_analyte$index, 
               y = effects_by_analyte$vloc, 
               col = effects_by_analyte$col,
               pch = 19)
        axis(4, at = ytick_locs, labels = ytick_vals, xpd = NA, las = 1)
        text(x = pusr[2] + diff(pusr[1:2]) / 4, y = mean(pusr[3:4]), 
             labels = "magnitude of relative effect", srt = 270, xpd = NA)
        legend("topright",
               legend = c("# non-zero folds", "Positive effect", "Negative effect"),
               col = c("black", cols),
               lty = c(1, NA, NA),
               lwd = c(3, NA, NA),
               pch = c(NA, 19, 19),
               pt.cex = 1, box.lty = 2, cex = 0.5)
      }
      
      if(length(nz_coefs) > 0){
        multihist(nz_coefs, stacked_hist = T, legend_location = "above", 
                  horiz_leg = F, xpd = NA, leg_loc_y_extra = 0.2, 
                  xlab = "standardized coefficient", freq = T, 
                  plot_cols = RColorBrewer::brewer.pal(G, name = "Dark2"))  
      } else {
        plot.new()
        plot.window(0:1, 0:1)
        text(0.5, 0.5, labels = paste0("no non-zero coefs at all"))
      }
      
      dev.off()
      
      #convert to png too
      magick::image_write(magick::image_read_pdf(fig_path, density = 300), 
                          path = gsub("\\.pdf$", ".png", fig_path), format = "png")
      
    }
  }
}

#combine all pdfs from run into one
all_figures_path <- paste0(figures_dir, "clinical-x-omics-prediction_", 
                           paste0(omes, collapse = "|"), "-", 
                           paste0(tissues, collapse = "|"), "-", 
                           target_phenotype, "-", 
                           paste0(methods, collapse = "|"),
                           ifelse(sex_interactions, "_sex-interaction", ""), 
                           ifelse(residualize_molecular_feats, "_ortho-mol", ""), 
                           ".pdf")
pdftools::pdf_combine(input = fig_paths, output = all_figures_path)
