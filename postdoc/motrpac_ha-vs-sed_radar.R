library(data.table)
library(plotrix)
library(smacof)
library(ape)

#### functions ####
source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")
source("~/scripts/minor_scripts/postdoc/circular_mds.R")
source("~/scripts/minor_scripts/postdoc/1Drepel_v2.R")
source("~/repos/flexible-sum-to-zero/Nik/R/functions.R")

#### high-level parameters ####
fig_dir <- "~/motrpac_ha-vs-sed/"
pres_fig_dir <- paste0(fig_dir, "presentation/")
do_signflip <- F
focal_sex_opts <- c("Male", "Female", "Both", "MF")[3]
use_combined_sample_for_cmat <- T
use_combined_within_sex_for_cmat <- T

#### download data ####
data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_sets/"
# data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/raw/data_sets/"
dir.create(data_dir, showWarnings = FALSE, recursive = T)
path_to_file <- "gs://motrpac-data-hub/human-main/phenotype/human-main-ha-adu/curated/data_sets/*"
# path_to_file <- "gs://motrpac-data-hub/human-main/phenotype/human-main-ha-adu/raw/data_sets/*"
# system(paste0("gsutil -m cp -r ", path_to_file, " ", data_dir))

#read in data
sd_ph <- as.data.frame(fread("~/data/motrpac_precovid_sedentary_clinical_obs.txt"))
sd_ph$avg_vo2cart <- (sd_ph$vo2cart1 + sd_ph$vo2cart2) / 2
ha_path <- paste0(data_dir, "human-main-ha-adu_clinical_screening_ds_crf-curated_v1.0.txt")
# ha_ph <- as.data.frame(fread(ha_path))
# ha_ph$lowerbody_volume
# cat(paste0("\"", colnames(sd_ph), "\","))

#### Sam's data file ####
#read in Sam's integrated data
sam_path <- "~/data/HA_SED_screening_draft1.txt"
ph <- as.data.frame(fread(sam_path))

#plot missing values?
cairo_pdf("~/HA_SED_screening_draft1-missing.pdf",
          width = 300/72, height = 500/72, pointsize = 12)
par(mar = c(2,5,2,2))
missing_prop <- apply(is.na(ph), 2, mean)
plot(missing_prop*0, missing_prop, xlim = c(0,0.55),
     xaxt = "n", xlab = "", yaxt = "n", ylab = "missing-ness proportion",  bty="n",
     pch = 10, cex = 0.25)
axis(2)

#label
labcex <- 0.45
minsep <- strheight("W") * 1.25 * labcex
xdisp <- strwidth("w") * labcex
laby <- minsep_fit(missing_prop, w = minsep, bounds = c(-0.05,1.05))
text(0.4, y = laby, labels = names(missing_prop), xpd = NA, cex = labcex, pos = 4)
elbow_l <- 0.05  
for(xi in 1:length(laby)){
  bezier_curve(y0 = missing_prop[xi], x0 = 0, y1 = laby[xi], x1 = 0.4 + xdisp * 1.5, xpd = NA, 
               ends = "flat")
}

dev.off()

#### consolidate data ####

#subset sam's data to baseline
ph <- ph[,missing_prop < 1]
ph$visit_code
sd_ph$age_years
sd_ph$bmi
intersect(sd_ph$pid, ph$pid)
intersect(colnames(ph), colnames(sd_ph))

#my old file
str(unique(sd_ph$pid))
#sam's file
str(unique(ph$pid))
intersect(ph$pid, sd_ph$pid)

sd_ph_sam <- ph[ph$cohort == "SED",]
best_idx <- lapply(seq_len(nrow(sd_ph)), function(i) {
  d2 <- (sd_ph_sam$wtkg - sd_ph$body_weight[i])^2 +
    (sd_ph_sam$htcm - sd_ph$body_height_cm[i])^2 +
    (sd_ph_sam$age_psca - sd_ph$age_years[i])^2
  d2[sd_ph_sam$sex != sd_ph$sex[i]] <- Inf
  which(d2 == min(d2))
})
if(all(sapply(best_idx, length) == 1)){
  best_idx <- unlist(best_idx)
} else {
  stop("multiple indices matched")
}
all(cbind(sd_ph[,"body_height_cm"] - sd_ph_sam[best_idx, "htcm"],
          sd_ph[,"body_weight"] - sd_ph_sam[best_idx, "wtkg"],
          sd_ph[,"age_years"] - sd_ph_sam[best_idx, "age_psca"]) == 0)

sd_ph_sam_sub <- sd_ph_sam[best_idx,]
sd_ph$leg_ext_1rm_kg - sd_ph_sam$le_rest_rmat[best_idx]

# function to find best-correlated column in sd_ph_sam_sub for each sd_ph column
best_matches <- sapply(names(sd_ph), function(col_ph) {
  x <- sd_ph[[col_ph]]
  if (is.numeric(x)) {
    cors <- sapply(sd_ph_sam_sub, function(y) {
      if (is.numeric(y) & !all(is.na(y))) cor(x, y, use = "complete.obs") else NA
    })
    names(which.max(abs(cors)))
  } else {
    # for non-numeric, measure concordance
    agrees <- sapply(sd_ph_sam_sub, function(y) {
      if (!all(is.na(y)) & length(unique(y)) > 1 & length(unique(x)) >1)
        fisher.test(x = as.character(x), 
                    y = as.character(y), simulate.p.value = T)$p.value
      else NA
    })
    names(which.min(agrees))
  }
})
unlist(best_matches)


#can just try using Sam's data alone?
sam_colnames <- setNames(colnames(ph), colnames(ph))
data_totals <- do.call(rbind, lapply(sam_colnames, function(scn){
  counts <- sapply(split(ph[,scn], ph$cohort), function(x) sum(!is.na(x)))
  out <- data.frame(feature = scn,
                    total = sum(counts),
                    ngroups = sum(counts > 0),
                    t(counts))
  out
}))
rownames(data_totals) <- NULL
split(data_totals, data_totals$ngroups)
split(paste0(data_totals$feature, " (SED = ", data_totals$SED,
             ", HARE = ", data_totals$HA.RE,
             ", HAEE = ", data_totals$HA.EE, ")"
), data_totals$ngroups)

#nope
#ok
#let's just try doing the matching manually
col_map <- c(
  pid                 = "pid",                 # numeric/char ok; ngroups=2
  study_label         = "study",               # integer code; recode to label later if needed
  sex                 = "sex",                 # character; ngroups=3
  age_years           = "age_psca",            # integer; ngroups=3
  visit_code          = "visit_code",          # character; ngroups=2
  body_weight         = "wtkg",                # numeric; ngroups=3
  body_height_cm      = "htcm",                # numeric; ngroups=3
  bmi                 = "bmi",                 # numeric; ngroups=3
  waist_circum_cm     = "wccm",                # numeric; ngroups=2
  bp_sys              = "sys_bp",              # numeric; ngroups=2
  bp_dia              = "dias_bp",             # numeric; ngroups=2
  pulse_rest          = "hr_min",              # integer; ngroups=2
  avg_vo2cart        = "avg_vo2cart_cpet",    # numeric; ngroups=2
  watts_max           = "wattend_cpet",        # integer; ngroups=2
  wtkg_cpet           = "wtkg_cpet",           # numeric; ngroups=3
  hr_end_cpet         = "hrend_cpet",          # integer; ngroups=2
  knee_iso_torque     = "peaktorq_iske",       # numeric; ngroups=3
  leg_press_1rm_lbs   = "lp_rest_corrected_rmat", # numeric; ngroups=1 fallback (avoid logical)
  leg_ext_1rm_lbs     = "le_rest_corrected_rmat", # numeric; ngroups=1 fallback
  chest_press_1rm_lbs = "cp_rest_corrected_rmat", # numeric; ngroups=1 fallback
  vo2_peak_L          = "VO2max_L",            # numeric; ngroups=3
  vo2_peak_mlkg       = "VO2max",              # numeric; ngroups=3
  o2_pulse            = "O2_pulse",            # numeric; ngroups=2
  handgrip_max        = "peak_handgrip"        # integer; ngroups=3
)

sd_ph <- sd_ph[, names(col_map)]
re_ph <- ph[ph$cohort == "HA-RE", col_map]
ee_ph <- ph[ph$cohort == "HA-EE", col_map]
# ee_ph[,colMeans(is.na(ee_ph)) == 1] <- 0
colnames(re_ph) <- colnames(ee_ph) <- names(col_map)

#quickly check sample means
rbind(
  SED = round(apply(sd_ph, 2, function(x) mean(as.numeric(x), na.rm = T)), 2),
  EE = round(apply(ee_ph, 2, function(x) mean(as.numeric(x), na.rm = T)), 2),
  RE = round(apply(re_ph, 2, function(x) mean(as.numeric(x), na.rm = T)), 2)
)

#### filter to focal variables ####

evs <- list()
mean_vals_list <- list()
if(focal_sex_opts == "MF"){
  focal_sexes <- c("Male", "Female")
} else {
  focal_sexes <- focal_sex_opts
}

#specify plot parameters
variables_to_inspect <- c("bmi", "waist_circum_cm", "bp_sys", "bp_dia", "pulse_rest", 
                          "vo2cart1", "vo2cart2", "watts_max", "wtkg_cpet", "hr_end_cpet", 
                          "grip1", "grip2", "grip3", "knee_iso_torque", "leg_press_1rm_lbs", 
                          "leg_ext_1rm_lbs", "chest_press_1rm_lbs", "lean_mass", "fat_mass", 
                          "total_bmd", "vat_mass", "vo2_peak_L", "vo2_peak_mlkg", "o2_pulse", 
                          "leg_press_1rm_kg", "leg_ext_1rm_kg", "chest_press_1rm_kg", "handgrip_max")
variables_to_inspect <- intersect(variables_to_inspect, colnames(sd_ph))
unit_names <- c(
  chest_press_1rm_lbs = "lbs",
  pulse_rest          = "bpm",
  hr_end_cpet         = "bpm",
  bp_dia              = "mmHg",
  bp_sys              = "mmHg",
  leg_press_1rm_lbs   = "lbs",
  waist_circum_cm     = "cm",
  bmi                 = "kg/m²",
  wtkg_cpet           = "kg",
  leg_ext_1rm_lbs     = "lbs",
  knee_iso_torque     = "N·m",
  handgrip_max        = "kg",
  o2_pulse            = "mL/beat",
  vo2_peak_L          = "L/min",
  watts_max           = "W",
  vo2_peak_mlkg       = "mL/kg/min"
)

axis_labels <- c(
  bmi                 = "BMI",
  waist_circum_cm     = "Waist Circumference",
  bp_sys              = "Systolic BP",
  bp_dia              = "Diastolic BP",
  pulse_rest          = "Resting HR",
  watts_max           = "Peak Power",
  wtkg_cpet           = "Body Mass (CPET)",
  hr_end_cpet         = "Peak Exercise HR",
  knee_iso_torque     = "Knee Ext. Torque",
  leg_press_1rm_lbs   = "Leg Press 1RM",
  leg_ext_1rm_lbs     = "Leg Ext. 1RM",
  chest_press_1rm_lbs = "Chest Press 1RM",
  vo2_peak_L          = "V̇O₂peak (absolute)",
  vo2_peak_mlkg       = "V̇O₂peak (relative)",
  o2_pulse            = "O₂ Pulse",
  handgrip_max        = "Max Handgrip"
)

#preprocess data

#fix vo2 numbers from ml to L
# re_ph$vo2_peak_L <- re_ph$vo2_peak_L * 1E3
# ee_ph$vo2_peak_L <- ee_ph$vo2_peak_L * 1E3
# re_ph$vo2_peak_mlkg <- re_ph$vo2_peak_mlkg * 1E3
# ee_ph$vo2_peak_mlkg <- ee_ph$vo2_peak_mlkg * 1E3
sd_ph$vo2_peak_L <- sd_ph$vo2_peak_L / 1E3
sd_ph$vo2_peak_mlkg <- sd_ph$vo2_peak_mlkg / 1E3

#compare means to ensure they are all on roughly the same order of magnitude
data.frame(units = unit_names[variables_to_inspect],
           SED = apply(sd_ph[,variables_to_inspect], 2, function(x) mean(as.numeric(x), na.rm = T)),
           HARE = apply(re_ph[,variables_to_inspect], 2, function(x) mean(as.numeric(x), na.rm = T)),
           HAEE = apply(ee_ph[,variables_to_inspect], 2, function(x) mean(as.numeric(x), na.rm = T)))

#save complete data frames before filtration
sd_ph_all <- sd_ph; re_ph_all <- re_ph; ee_ph_all <- ee_ph

for(focal_sex in focal_sexes){
  
  if(focal_sex == "Both"){
    variables_to_condition_on <- c("sex", "age_years", "body_weight", "body_height_cm")
  } else {
    variables_to_condition_on <- c("age_years", "body_weight", "body_height_cm")
  }
  
  #filter by sex (or make sex a numeric variable)
  sd_ph_all -> sd_ph; re_ph_all -> re_ph; ee_ph_all -> ee_ph
  if(focal_sex == "Both"){
    sd_ph$sex <- (sd_ph$sex == "Female") * 1
    re_ph$sex <- (re_ph$sex == "Female") * 1
    ee_ph$sex <- (ee_ph$sex == "Female") * 1
  } else {
    sd_ph <- sd_ph[sd_ph$sex == focal_sex,]
    re_ph <- re_ph[re_ph$sex == focal_sex,]
    ee_ph <- ee_ph[ee_ph$sex == focal_sex,]
  }
  
  #get correlation matrix
  cmat <- cor(sd_ph[,variables_to_inspect], use = "com")
  #or pooled version
  if(use_combined_sample_for_cmat){
    if(use_combined_within_sex_for_cmat){
      cv_list <- lapply(list(SED_m = sd_ph_all[sd_ph_all$sex == "Male",variables_to_inspect], 
                             SED_f = sd_ph_all[sd_ph_all$sex == "Female",variables_to_inspect], 
                             HAEE_m = ee_ph_all[ee_ph_all$sex == "Male",variables_to_inspect], 
                             HAEE_f = ee_ph_all[ee_ph_all$sex == "Female",variables_to_inspect], 
                             HARE_m = re_ph_all[re_ph_all$sex == "Male",variables_to_inspect], 
                             HARE_f = re_ph_all[re_ph_all$sex == "Female",variables_to_inspect]),
                        cov, use = "pairwise.complete.obs") 
    } else {
      cv_list <- lapply(list(SED = sd_ph_all[,variables_to_inspect], 
                             HAEE = ee_ph_all[,variables_to_inspect], 
                             HARE = re_ph_all[,variables_to_inspect]), 
                        cov, use = "pairwise.complete.obs")
    }
  } else {
    cv_list <- lapply(list(SED = sd_ph[,variables_to_inspect], 
                           HAEE = ee_ph[,variables_to_inspect], 
                           HARE = re_ph[,variables_to_inspect]), 
                      cov, use = "pairwise.complete.obs")
  }
  if(use_combined_within_sex_for_cmat){
    cv_list[["HAEE_m"]][is.na(cv_list[["HAEE_m"]])] <- 0
    cv_list[["HAEE_f"]][is.na(cv_list[["HAEE_f"]])] <- 0
    cv_pooled <- Reduce("+", cv_list) / ((cv_list[["HAEE_m"]] != 0) + 
                                           (cv_list[["HAEE_f"]] != 0) + 4)
  } else {
    cv_list[["HAEE"]][is.na(cv_list[["HAEE"]])] <- 0
    cv_pooled <- Reduce("+", cv_list) / ((cv_list[["HAEE"]] != 0) + 2)
  }
  
  # convert to correlation matrix
  isds <- diag(1 / sqrt(diag(cv_pooled)))
  cmat <- isds %*% cv_pooled %*% isds
  
  #strongest correlations?
  top_elements <- arrayInd(order(cmat, decreasing = T), 
                           .dim = dim(cmat))[-c(1:nrow(cmat)),]
  top_elements <- top_elements[(1:(nrow(top_elements)/2))*2,]
  data.frame(var1 = variables_to_inspect[top_elements[,1]], 
             var2 = variables_to_inspect[top_elements[,2]], 
             r = round(cmat[top_elements], 3))
  
  # cvmat <- cov(sd_ph[,variables_to_inspect], use = "com")
  # L_cvmat <- chol(Matrix::nearPD(x = cvmat)$mat)
  
  #flip signs to positivity
  if(do_signflip){
    signflips <- flip_signs(cmat, objective = "count")$signs
    if(mean(signflips == -1) > 0.5){signflips <- -signflips}
    sd_ph[,variables_to_inspect] <- t(t(as.matrix(sd_ph[,variables_to_inspect])) * signflips)
    re_ph[, variables_to_inspect] <- t(t(as.matrix(re_ph[, variables_to_inspect])) * signflips)
    ee_ph[, variables_to_inspect] <- t(t(as.matrix(ee_ph[, variables_to_inspect])) * signflips)
    cmat <- diag(signflips) %*% cmat %*% t(diag(signflips))
  } else {
    signflips <- rep(1, ncol(cmat))
  }
  
  #get names back in there
  rownames(cmat) <- colnames(cmat) <- rownames(cv_pooled)
  
  #compute distances
  dists <- 1-cmat^2
  dists <- dists / max(dists) * pi
  
  #### fit models ####
  
  # compute mean values of conditioning variables
  # mean_vals <- cbind(sapply(sd_ph[ , variables_to_condition_on], 
  #                     function(x) mean(x, na.rm = TRUE)),
  #                    sapply(re_ph[ , variables_to_condition_on], 
  #                           function(x) mean(x, na.rm = TRUE)),
  #                    sapply(ee_ph[ , variables_to_condition_on], 
  #                           function(x) mean(x, na.rm = TRUE)))
  mean_vals <- cbind(SED = sapply(variables_to_condition_on, 
                                  function(v2c) mean(
                                    sapply(split(sd_ph[,v2c], sd_ph$sex), mean),
                                    na.rm = TRUE)),
                     HARE = sapply(variables_to_condition_on, 
                                   function(v2c) mean(
                                     sapply(split(re_ph[,v2c], re_ph$sex), mean),
                                     na.rm = TRUE)),
                     HAEE = sapply(variables_to_condition_on, 
                                   function(v2c) mean(
                                     sapply(split(ee_ph[,v2c], ee_ph$sex), mean),
                                     na.rm = TRUE)))
  
  mean_vals <- rowMeans(mean_vals)
  mean_vals <- as.data.frame(as.list(mean_vals))
  mean_vals_list[[focal_sex]] <- mean_vals
  
  backend <- c("dglm", "homoskedastic")[2]
  
  # ensure ev exists with rows in the same order as variables_to_inspect
  if (!exists("ev") || !all(ev$variable == variables_to_inspect)) {
    ev <- data.frame(variable = variables_to_inspect, stringsAsFactors = FALSE)
  }
  
  # store means/sds (optional, for inspection)
  ev$SED_mu <- ev$HARE_mu <-ev$HAEE_mu <-ev$SED_sd <-ev$HARE_sd <-ev$HAEE_sd <- 
    ev$SED_rad <- ev$HARE_rad <- ev$HAEE_rad <- NA_real_
  n_missing <- c(SED = 0L, HARE = 0L, HAEE = 0L)
  pair_names <- c("HARE_in_SED","HAEE_in_SED","SED_in_HARE","SED_in_HAEE","HARE_in_HAEE","HAEE_in_HARE")
  for (nm in pair_names) {
    qnm <- paste0("q_", nm)
    znm <- paste0("z_", nm)
    if (!qnm %in% names(ev)) ev[[qnm]] <- NA_real_
    if (!znm %in% names(ev)) ev[[znm]] <- NA_real_
  }
  
  for (j in seq_along(variables_to_inspect)) {
    feat <- variables_to_inspect[j]
    
    s_sed  <- fit_one(sd_ph,  feat, variables_to_condition_on, mean_vals, backend = backend, 
                      interact_with = ifelse(focal_sex == "Both", "sex", ""))
    s_hare <- fit_one(re_ph,  feat, variables_to_condition_on, mean_vals, backend = backend, 
                      interact_with = ifelse(focal_sex == "Both", "sex", ""))
    s_haee <- fit_one(ee_ph,  feat, variables_to_condition_on, mean_vals, backend = backend, 
                      interact_with = ifelse(focal_sex == "Both", "sex", ""))
    
    if (s_sed$n == 0) n_missing["SED"] <- n_missing["SED"] + 1L
    if (s_hare$n == 0) n_missing["HARE"] <- n_missing["HARE"] + 1L
    if (s_haee$n == 0) n_missing["HAEE"] <- n_missing["HAEE"] + 1L
    
    F_sed  <- if (s_sed$n  > 0) make_norm_cdf(s_sed$sd)  else NULL
    F_hare <- if (s_hare$n > 0) make_norm_cdf(s_hare$sd) else NULL
    F_haee <- if (s_haee$n > 0) make_norm_cdf(s_haee$sd) else NULL
    
    ev$SED_mu[j]  <- s_sed$mu
    ev$HARE_mu[j] <- s_hare$mu
    ev$HAEE_mu[j] <- s_haee$mu
    ev$SED_sd[j]  <- s_sed$sd
    ev$HARE_sd[j] <- s_hare$sd
    ev$HAEE_sd[j] <- s_haee$sd
    
    # HARE vs SED
    if (!is.null(F_sed) && is.finite(s_hare$mu) && is.finite(s_sed$mu)) {
      ev$HARE_rad[j] <- F_sed(s_hare$mu - s_sed$mu)
    } else {
      ev$HARE_rad[j] <- NA_real_
    }
    
    # HAEE vs SED
    if (!is.null(F_sed) && is.finite(s_haee$mu) && is.finite(s_sed$mu)) {
      ev$HAEE_rad[j] <- F_sed(s_haee$mu - s_sed$mu)
    } else {
      ev$HAEE_rad[j] <- NA_real_
    }
    
    # SED vs HA midpoint: average of available parts
    q_parts <- numeric(0)
    if (!is.null(F_hare) && is.finite(s_sed$mu) && is.finite(s_hare$mu)) {
      q_parts <- c(q_parts, F_hare(s_sed$mu - s_hare$mu))
    }
    if (!is.null(F_haee) && is.finite(s_sed$mu) && is.finite(s_haee$mu)) {
      q_parts <- c(q_parts, F_haee(s_sed$mu - s_haee$mu))
    }
    ev$SED_rad[j] <- if (length(q_parts)) mean(q_parts) else NA_real_
    
    #all pairwise comparisons
    
    # HARE vs SED
    ev$q_HARE_in_SED[j]  <- .safe_q(s_hare$mu, s_sed$mu,  F_sed)
    ev$z_HARE_in_SED[j]  <- .safe_z(s_hare$mu, s_sed$mu,  s_sed$sd)
    
    # HAEE vs SED
    ev$q_HAEE_in_SED[j]  <- .safe_q(s_haee$mu, s_sed$mu,  F_sed)
    ev$z_HAEE_in_SED[j]  <- .safe_z(s_haee$mu, s_sed$mu,  s_sed$sd)
    
    # SED vs HARE
    ev$q_SED_in_HARE[j]  <- .safe_q(s_sed$mu,  s_hare$mu, F_hare)
    ev$z_SED_in_HARE[j]  <- .safe_z(s_sed$mu,  s_hare$mu, s_hare$sd)
    
    # SED vs HAEE
    ev$q_SED_in_HAEE[j]  <- .safe_q(s_sed$mu,  s_haee$mu, F_haee)
    ev$z_SED_in_HAEE[j]  <- .safe_z(s_sed$mu,  s_haee$mu, s_haee$sd)
    
    # HARE vs HAEE
    ev$q_HARE_in_HAEE[j] <- .safe_q(s_hare$mu, s_haee$mu, F_haee)
    ev$z_HARE_in_HAEE[j] <- .safe_z(s_hare$mu, s_haee$mu, s_haee$sd)
    
    # HAEE vs HARE
    ev$q_HAEE_in_HARE[j] <- .safe_q(s_haee$mu, s_hare$mu, F_hare)
    ev$z_HAEE_in_HARE[j] <- .safe_z(s_haee$mu, s_hare$mu, s_hare$sd)
    
  }
  
  # clamp to [0,1]
  ev$SED_rad  <- pmin(1, pmax(0, ev$SED_rad))
  ev$HARE_rad <- pmin(1, pmax(0, ev$HARE_rad))
  ev$HAEE_rad <- pmin(1, pmax(0, ev$HAEE_rad))
  
  
  # # initialize output storage
  # ev <- data.frame(variable = variables_to_inspect, expected = NA_real_)
  # 
  # # loop over each outcome variable to find expected values at mean
  # for (i in seq_along(variables_to_inspect)) {
  #   outcome <- variables_to_inspect[i]
  #   fmla <- as.formula(paste(outcome, "~", paste(variables_to_condition_on, collapse = " + ")))
  #   fit <- lm(fmla, data = sd_ph)
  #   
  #   #expected val and residual sd
  #   pred_val <- predict(fit, newdata = mean_vals)
  #   ev$SED_exp[i] <- pred_val
  #   ev$SED_sd[i] <- sigma(fit)
  # }
  # 
  # #placeholder variables for now
  # ev$HARE_sd <- ev$HAEE_sd <- ev$SED_sd
  # # ev$HARE_exp <- ev$SED_exp + as.numeric(mvtnorm::rmvnorm(n = 1, sigma = cvmat))
  # # ev$HAEE_exp <- ev$SED_exp + as.numeric(mvtnorm::rmvnorm(n = 1, sigma = cvmat))
  # ev$HARE_exp <- ev$SED_exp + rnorm(length(ev$SED_exp)) * ev$HARE_sd / 2
  # ev$HAEE_exp <- ev$SED_exp + rnorm(length(ev$SED_exp)) * ev$HAEE_sd / 2
  # 
  # #get glass delta-style values
  # exp_names <- c("SED_exp", "HARE_exp", "HARE_exp")
  # exp_mins <- apply(ev[,exp_names], 1, min)
  # exp_maxs <- apply(ev[,exp_names], 1, max)
  # 
  # ev$axis_centers <- rowMeans(ev[,exp_names])
  # ev$asix_sds <- rowMeans(ev[,gsub("exp", "sd", exp_names)])
  # ev$axis_lb <- ev$axis_centers - 2 * ev$asix_sds
  # ev$axis_lb[ev$axis_lb > exp_mins] <- exp_mins[ev$axis_lb > exp_mins]
  # ev$axis_ub <- ev$axis_centers + 2 * ev$asix_sds
  # ev$axis_ub[ev$axis_ub < exp_maxs] <- exp_mins[ev$axis_ub < exp_maxs]
  # 
  # ev$SED_rad <- (ev$SED_exp - ev$axis_lb) / (ev$axis_ub - ev$axis_lb)
  # ev$HARE_rad <- (ev$HARE_exp - ev$axis_lb) / (ev$axis_ub - ev$axis_lb)
  # ev$HAEE_rad <- (ev$HAEE_exp - ev$axis_lb) / (ev$axis_ub - ev$axis_lb)
  
  #get coordinates on circle
  cmds <- circ_mds(dists, maxit = 1E5, n_starts = 10)
  ev$angle <- atan2(y = cmds$coords$y, x = cmds$coords$x)
  
  ev$SED_x  <- ev$SED_rad  * cos(ev$angle)
  ev$SED_y  <- ev$SED_rad  * sin(ev$angle)
  ev$HARE_x <- ev$HARE_rad * cos(ev$angle)
  ev$HARE_y <- ev$HARE_rad * sin(ev$angle)
  ev$HAEE_x <- ev$HAEE_rad * cos(ev$angle)
  ev$HAEE_y <- ev$HAEE_rad * sin(ev$angle)
  
  #reorder ev by circle order
  ev <- ev[order(ev$angle),]
  evs[[focal_sex]] <- ev
  
}

#### radar plot ####

x_base <- list(li = 1.5, lo = 0.35,
               ri = 3.5, ro = 4.65)

cairo_pdf("~/HAEE-HARE-SED_Exp-Comparison.pdf", width = 800/72, height = 800/72, 
          pointsize = 20)

# https://vetr.dev/layout-gui/#7.8|2x_base$ri;1133;1537;4x_base$ri|,,,
mat <- matrix(c(
  c(0, 0, 0, 0, 0, 0, 0, 0),
  c(0, 2, 2, 0, 0, 3, 3, 0),
  c(0, 2, 2, 1, 1, 3, 3, 0),
  c(0, 0, 0, 1, 1, 0, 0, 0),
  c(0, 0, 0, 4, 4, 0, 0, 0),
  c(0, 0, 0, 4, 4, 0, 0, 0),
  c(0, 0, 0, 0, 0, 0, 0, 0)
), nrow = 7, byrow = TRUE)

layout(mat)
par(mar = c(0,0,0,0))
plot(NULL, xlim = c(-2,2), ylim = c(-2,2), asp = 1, 
     xaxt = "n", xlab = "", yaxt = "n", ylab = "",  bty="n",
     cex = 0)
feat_labels <- setNames(paste0(c("","-")[(signflips==-1) + 1], rownames(cmat)),
                        rownames(cmat))
# feat_labels <- axis_labels[names(feat_labels)]

xyt <- c(0,0)
label_sunburst(cmds$coords, labels = feat_labels, 
               rotate = c("none","radial","tangent")[2], 
               cex = 0.75, xyt)

#draw lines from origin to pts
annot_lwd <- 0.75
segments(x0 = 0, y0 = 0, x1 = cmds$coords$x, y1 = cmds$coords$y, 
         lty = 3, lwd = annot_lwd, col = adjustcolor(1, 0.8))

#draw circle boundary
t <- seq(0, 2*pi, length.out = 128)
circle_coords <- data.frame(x = sin(t), y = cos(t))
lines(circle_coords, lwd = annot_lwd, col = adjustcolor(1, 0.3))

#draw the polygons
cols_border <- c(SED = "darkblue", 
                 HAEE = "darkred", 
                 HARE = "darkorange")
cols_fill <- sapply(cols_border, adjustcolor, 0.4)
panel_titles <- c(SED = "SEDENTARY", 
                  HAEE = "HIGHLY-ACTIVE (ENDURANCE)", 
                  HARE = "HIGHLY-ACTIVE (RESISTANCE)")

polygon(ev$SED_x, ev$SED_y, 
        col = cols_fill["SED"], border = cols_border["SED"])
polygon(ev$HARE_x, ev$HARE_y, 
        col = cols_fill["HARE"], border = cols_border["HARE"])
polygon(ev$HAEE_x, ev$HAEE_y, 
        col = cols_fill["HAEE"], border = cols_border["HAEE"])

#draw individual distributions
plot_base <- function(xyt = c(0,0)){
  cmds_tcoords <- data.frame(t(t(cmds$coords) + xyt))
  plot(NULL, xlim = c(-2,2), ylim = c(-2,2), asp = 1, 
       xaxt = "n", xlab = "", yaxt = "n", ylab = "",  bty="n",
       cex = 0)
  label_sunburst(cmds$coords, labels = feat_labels, 
                 rotate = c("none","radial","tangent")[2], 
                 cex = 0.75, xyt = xyt)
  segments(x0 = xyt[1], y0 = xyt[2], x1 = cmds_tcoords$x, y1 = cmds_tcoords$y, 
           lty = 3, lwd = annot_lwd, col = adjustcolor(1, 0.8))
  lines(t(t(circle_coords) + xyt), lwd = annot_lwd, col = adjustcolor(1, 0.3))
}

xyt1 <- c(0,1)
plot_base(xyt1)
polygon(ev$SED_x + xyt1[1], ev$SED_y + xyt1[2], 
        col = cols_fill["SED"], border = cols_border["SED"])
arctext(
  x = panel_titles["SED"],
  center = xyt1,
  radius = x_base$ri,
  middle = 3/4*pi,
  cex   = 0.9, 
  col = cols_fill["SED"],
  font = 2,
  xpd = NA
)

xyt2 <- c(0,1)
plot_base(xyt2)
polygon(ev$HAEE_x + xyt2[1], ev$HAEE_y + xyt2[2], 
        col = cols_fill["HAEE"], border = cols_border["HAEE"])
arctext(
  x = panel_titles["HAEE"],
  center = xyt2,
  radius = 3,
  middle = 1/4*pi + 1/24*pi,
  cex   = 0.9, 
  col = cols_fill["HAEE"],
  font = 2,
  xpd = NA
)

xyt3 <- c(0,-0.5)
plot_base(xyt = xyt3)
polygon(ev$HARE_x + xyt3[1], ev$HARE_y + xyt3[2], 
        col = cols_fill["HARE"], border = cols_border["HARE"])
arctext(
  x = panel_titles["HARE"],
  center = xyt3,
  radius = 3.25,
  cex   = 0.9, 
  col = cols_fill["HARE"],
  font = 2,
  xpd = NA, 
  middle = 3/2*pi, clockwise = F
)


dev.off()

#### linear MDS (viewpoint plot) ####
extra_x_poly <- 0.02
mds_method <- c("cailliez", "isoMDS", "smacof", "linear")[3]
dists <- dists / max(dists) #renorm to max dist
if(mds_method == "cailliez") {
  fit <- cmdscale(as.dist(dists), k = 1, add = T)
  if("list" %in% class(fit)) fit <- fit$points
  x1d <- setNames(as.numeric(fit), rownames(fit))
} else if(mds_method == "isoMDS") {
  fit <- MASS::isoMDS(as.dist(dists), k = 1)$points
  x1d <- setNames(as.numeric(fit), rownames(fit))
} else if(mds_method == "smacof") {
  fit <- smacof::smacofSym(as.dist(dists), ndim = 1, type = "ordinal")$conf
  x1d <- setNames(as.numeric(fit), rownames(fit))
} else if(mds_method == "linear") {
  x1d <- setNames(1:nrow(dists)/nrow(dists), rownames(dists))
}

ord_lin <- order(x1d)
x1d <- x1d[ord_lin]
x_plot  <- (x1d - min(x1d)) / diff(range(x1d))
labs_ord <- feat_labels[names(x_plot)]

#recover parameters for plotting
yvals <- list(
  SED = setNames(ev$SED_rad, ev$variable)[names(labs_ord)],
  HARE = setNames(ev$HARE_rad, ev$variable)[names(labs_ord)],
  HAEE = setNames(ev$HAEE_rad, ev$variable)[names(labs_ord)]
)
yvals_raw <- list(
  SED = setNames(ev$SED_mu, ev$variable)[names(labs_ord)],
  HARE = setNames(ev$HARE_mu, ev$variable)[names(labs_ord)],
  HAEE = setNames(ev$HAEE_mu, ev$variable)[names(labs_ord)]
)
yvals_raw <- lapply(yvals_raw, round, 2)

#plotting
panel_titles <- c(SED = "SEDENTARY", 
                  HAEE = "HIGHLY-ACTIVE\n(ENDURANCE)", 
                  HARE = "HIGHLY-ACTIVE\n(RESISTANCE)")
fig_path <- paste0(fig_dir, "HAEE-HARE-SED_Exp-Comparison_Linear-", mds_method,".pdf")
cairo_pdf(fig_path,
          width = 1000/72, height = 300/72, pointsize = 12)
par(mfrow = c(1, 4), oma = c(0,4,0,0))

####
pops <- c("SED", "HARE", "HAEE")[c(2,1,3)]
links <- c("bezier", "sharp")[1]
par(mar = c(11, 1, 11, 1))

for(pop in pops){
  
  #initialize plotting window
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "n")
  
  #draw lines to peaks corresponding to traits
  missing_vals <- is.na(yvals[[pop]])
  segments(x0 = x_plot[!missing_vals], 
           y0 = 0, 
           x1 = x_plot[!missing_vals], 
           y1 = yvals[[pop]][!missing_vals], 
           lty = 3, col = adjustcolor(1, 0.5))
  
  #draw mountrain range
  xvals_input <- c(-extra_x_poly, x_plot[!missing_vals], 1 + extra_x_poly)
  yvals_input <- c(0, yvals[[pop]][!missing_vals], 0)
  polyx <- c(xvals_input, rev(xvals_input))
  polyy <- c(yvals_input, rep(0, length(yvals_input)))
  polygon(polyx, polyy,
          col = cols_fill[pop], border = cols_border[pop])
  
  # --- inside-only horizontal midline at y = 0.5 (set to 0 if desired) ---
  y_mid <- 0.5  # change to 0 to draw along the baseline
  n <- length(polyx)
  if(n >= 3){
    xints <- numeric(0)
    for(k in seq_len(n)){
      i <- k
      j <- if (k < n) k + 1 else 1
      xi <- polyx[i]; xj <- polyx[j]
      yi <- polyy[i]; yj <- polyy[j]
      # does edge [i->j] cross the horizontal y = y_mid?
      if( (yi <= y_mid && y_mid < yj) || (yj <= y_mid && y_mid < yi) ){
        t <- (y_mid - yi) / (yj - yi)
        xints <- c(xints, xi + t * (xj - xi))
      }
    }
    if(length(xints) >= 2){
      xints <- sort(xints)
      for(m in seq(1, length(xints) - 1, by = 2)){
        segments(x0 = xints[m], y0 = y_mid,
                 x1 = xints[m + 1], y1 = y_mid,
                 lwd = 1, lty = 1, col = 1, xpd = NA)
      }
    }
  }
  
  
  #label panel title
  text(x = -0.05, y = c(SED = 1.6, HARE = 1.8, HAEE = 1.8)[pop], labels = panel_titles[pop], pos = 4,
       cex = 2.5, font = 2, col = cols_fill[pop], xpd = NA)
  
  if(pop == pops[1]){
    
    #get horiz label locations
    xyr <- xyrat()
    minsep_scale <- 1.25
    minsep <- strwidth("w") * minsep_scale
    ydisp <- strheight("w")
    labx <- seq(0, 1, length = length(x_plot)) #uniformly distributed
    labx <- minsep_fit(x_plot, w = minsep, bounds = c(-0.05,1.05))
    labx_peaks <- minsep_fit(x_plot, w = minsep, bounds = c(-0.02,1.00))
    elbow_l <- 0.05
    
    #draw vertical axis
    yax_xloc <- -0.04
    yax_locs <- 0:4/4
    yax_labs <- c("0", paste0("0.", 1:3*25), "1")
    segments(yax_xloc,0,yax_xloc,1, xpd = NA,lwd = 1.5)
    segments(yax_xloc,yax_locs,yax_xloc-0.02,yax_locs, xpd = NA, lwd = 1.5)
    text(x = yax_xloc-0.01, y = yax_locs, labels = yax_labs, pos = 2, xpd = NA,
         cex = 1.25)
    text(yax_xloc-0.17, y = 0.5, labels = "reciprocal quantile", 
         srt = 90, xpd = NA, cex = 1.25)
    
    #draw vertical axis guiding lines
    segments(-0.02,yax_locs,4.55,yax_locs, xpd = NA, 
             lwd = 0.75, col = adjustcolor(1, 0.15), lty = 2)
    # segments(-0.02,0.5,4.55,0.5, xpd = NA, 
    #          lwd = 1, col = adjustcolor(1, 0.25), lty = 1) #middle line
    
  }
  
  #draw horizontal labels
  if(links == "sharp"){
    segments(x_plot, 0, x_plot, -elbow_l, xpd = NA)
    segments(x_plot, -elbow_l, labx, -elbow_l*3, xpd = NA)
    segments(labx, -elbow_l*3, labx, -elbow_l*4, xpd = NA)  
  } else {
    for(xi in 1:length(labx)){
      if(!missing_vals[xi]){
        bezier_curve(x_plot[xi], y0 = 0, labx[xi], -elbow_l*4, xpd = NA, 
                     ends = "steep")
      }
      
    }
  }
  text(x = labx[!missing_vals] + minsep / minsep_scale, 
       y = -elbow_l*4 - ydisp/4, 
       labels = labs_ord[!missing_vals], 
       xpd = NA, srt = 90, pos = 2)
  
  #label peaks with true heights
  for(xi in 1:length(labx_peaks)){
    if(!missing_vals[xi]){
      bezier_curve(x_plot[xi], y0 = yvals[[pop]][xi], 
                   labx_peaks[xi], yvals[[pop]][xi] + elbow_l*4, xpd = NA, 
                   ends = "steep", xpd = NA)
      peak_lab <- paste0(pretty_large_number(yvals_raw[[pop]][names(labs_ord)][xi]),
                         " ", unit_names[names(labs_ord)[xi]])
      text(x = labx_peaks[xi] - minsep / minsep_scale / 2, 
           y =  yvals[[pop]][xi] + elbow_l*4 + ydisp/4, 
           labels = peak_lab, xpd = NA, srt = 90, pos = 4)
    }
  }
  
}

#get combined in there too
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n")
max_quantiles <- apply(do.call(cbind, yvals), 1, max, na.rm = T)
top_pop <- unlist(apply(do.call(cbind, yvals), 1, which.max))
segments(x0 = x_plot, y0 = 0, x1 = x_plot, 
         y1 = max_quantiles, lty = 3, 
         col = adjustcolor(1, 0.5))
for(pop in pops){
  missing_vals <- is.na(yvals[[pop]])
  xvals_input <- c(-extra_x_poly, x_plot[!missing_vals], 1 + extra_x_poly)
  yvals_input <- c(0, yvals[[pop]][!missing_vals], 0)
  polyx <- c(xvals_input, rev(xvals_input))
  polyy <- c(yvals_input, rep(0, length(yvals_input)))
  polygon(polyx, polyy,
          col = cols_fill[pop], border = cols_border[pop])
}
if(links == "sharp"){
  segments(x_plot, 0, x_plot, -elbow_l, xpd = NA)
  segments(x_plot, -elbow_l, labx, -elbow_l*3, xpd = NA)
  segments(labx, -elbow_l*3, labx, -elbow_l*4, xpd = NA)  
} else {
  for(xi in 1:length(labx)){
    bezier_curve(x_plot[xi], y0 = 0, labx[xi], -elbow_l*4, xpd = NA, 
                 ends = "steep")
  }
}
text(x = labx + minsep / minsep_scale, y = -elbow_l*4 - ydisp/4, 
     labels = labs_ord, 
     xpd = NA, srt = 90, pos = 2)

text(x = -0.05, y = 1.45, labels = "COMBINED", pos = 4,
     cex = 2.5, font = 2, col = adjustcolor(1, 0.4), xpd = NA)

#label top pop
segments(x0 = x_plot, y0 = max_quantiles, x1 = x_plot, 
         y1 = 1.0, lty = 1, cols_border[pops[top_pop]])
for(xi in 1:length(labx)){
  bezier_curve(x_plot[xi], y0 = 1, labx[xi], 1.05, xpd = NA, 
               ends = "steep", col = cols_border[pops[top_pop[xi]]])
}
text(x = labx - minsep / minsep_scale / 2, y = 1.06, 
     labels = pops[top_pop], col = cols_border[pops[top_pop]],
     xpd = NA, srt = 90, pos = 4)


dev.off()

magick::image_write(magick::image_read_pdf(fig_path, density = 300), 
                    path = gsub("\\.pdf$", ".png", fig_path), format = "png")

#### linear opposing views ####

#prep data objects
y_plot <- x_plot
pops <- list(li = c("HAEE", "HARE"),
             lo = c("HAEE", "SED"),
             ri = c("HARE", "HAEE"),
             ro = c("HARE", "SED")
)
width_peaks <- 1
links <- c("bezier", "sharp")[1]
values_to_use <- c("q", "z")[2]
line_lwd <- 1
panel_titles <- c(SED = "SEDENTARY", 
                  HAEE = "HIGHLY-ACTIVE (ENDURANCE)", 
                  HARE = "HIGHLY-ACTIVE (RESISTANCE)")
MFplot <- focal_sex_opts == "MF"
sex_cols <- c("Male" = "darkolivegreen", Female = "darkmagenta")
sex_cols_fill <- sapply(sex_cols, adjustcolor, 0.4)
col_miss_fill <- adjustcolor("grey", 0.3)
col_miss_border <- "grey40"
shrink_gray_poly_y_by <- 0.015
lty_miss <- 2

# --- z-mode setup (robust symmetric scaling to [0,1]) ---
z_cols <- grep("^z_.*_in_.*$", names(ev), value = TRUE)
z_all  <- unlist(ev[z_cols], use.names = FALSE)
z_all  <- z_all[is.finite(z_all)]
zr     <- if(length(z_all)) as.numeric(quantile(abs(z_all), 0.98, na.rm = TRUE)) else 3
zr <- ceiling(zr)
fig_path <- paste0(fig_dir, "HAEE-HARE-SED_Exp-Comparison_Linear_compact-", 
                   values_to_use, 
                   ifelse(MFplot, "_sex-stratified", ""),".pdf")
cairo_pdf(fig_path,
          width = 750/72, height = 400/72, pointsize = 12)
par(mfrow = c(1, 1), oma = c(0,0,0,0), mar = c(4,0,4,0))

#initialize plotting window
plot(NA, xlim = c(-2, 7), ylim = c(0, 1), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", bty = "n")

#first get vertical label locations
xyr <- xyrat()
minsep_scale <- 1.5
minsep <- strheight("W") * minsep_scale
laby <- seq(0, 1, length = length(y_plot)) #uniformly distributed
laby <- minsep_fit(y_plot, w = minsep, bounds = c(-0.025,1.025))
if(MFplot){
  peak_lab_cex <- 0.8
  laby_peaks_info <- minsep_fit(y_plot, w = minsep * peak_lab_cex, 
                                bounds = c(-0.01,c(q=1.02,z=1.0)[values_to_use]), 
                                resize = T)
  if(class(laby_peaks_info) == "numeric"){
    laby_peaks <- laby_peaks_info
  } else {
    laby_peaks <- laby_peaks_info$y
    peak_lab_cex <- laby_peaks_info$scale_factor * 1.25  
  }
  elbow_l <- 0.05
  
} else {
  laby_peaks <- minsep_fit(y_plot, w = minsep, 
                           bounds = c(-0.01,c(q=1.02,z=1.0)[values_to_use]))
  peak_lab_cex <- 1
  elbow_l <- 0.05
}
xdisp <- strwidth("w") * peak_lab_cex

#plot vertical axis labels
midp_x <- 2.5
text(y = laby, 
     x = midp_x, 
     labels = axis_labels[labs_ord], 
     xpd = NA)

#draw line out to maximum width label
nfeat <- length(y_plot)
lab_w <- strwidth(axis_labels[labs_ord]) + xdisp
max_width_label <- max(lab_w) + xdisp
bezier_insertions <- midp_x + c(1,-1) * max(lab_w) / 2
missing_left <- is.na(ev[match(names(y_plot), ev$variable), paste0(pops$li[1], "_mu")]) 
missing_right <- is.na(ev[match(names(y_plot), ev$variable), paste0(pops$ri[1], "_mu")])
segments(x0 = c(midp_x + lab_w[!missing_right] / 2, midp_x - lab_w[!missing_left] / 2),
         x1 = c(rep(bezier_insertions[1], sum(!missing_right)), 
                rep(bezier_insertions[2], sum(!missing_left))),
         y0 = c(laby[!missing_right], laby[!missing_left]), 
         y1 = c(laby[!missing_right], laby[!missing_left]), xpd = NA
)

#draw lines out to vertical axes 
for(xi in 1:length(laby)){
  
  #right curve
  if(!missing_right[xi]){
    bezier_curve(y1 = y_plot[xi], x0 = bezier_insertions[1], 
                 y0 = laby[xi], x1 = x_base$ri, 
                 xpd = NA, 
                 ends = "flat")  
  }
  
  #left curve
  if(!missing_left[xi]){
    bezier_curve(y1 = y_plot[xi], x0 = bezier_insertions[2], 
                 y0 = laby[xi], x1 = x_base$li, 
                 xpd = NA, 
                 ends = "flat")
  }
}

#draw horizontal axes
xax_xloc <- 1.04
xax_locs <- 0:4/4
if (values_to_use == "q") {
  xax_labs <- c("0", ".25", ".5", ".75", "1")
} else {
  z_ticks <- c(-zr, -zr/2, 0, zr/2, zr)
  xax_labs <- as.character(round(z_ticks, 2))
}

#inner
segments(x0 = x_base$li, x1 = x_base$li - width_peaks, 
         y0 = xax_xloc, y1 = xax_xloc, col = 1, xpd = NA, lwd = line_lwd)
segments(x0 = x_base$ri, x1 = x_base$ri + width_peaks, 
         y0 = xax_xloc, y1 = xax_xloc, col = 1, xpd = NA, lwd = line_lwd)
segments(x0 = xax_locs + x_base$ri, x1 = xax_locs + x_base$ri,
         y0 = xax_xloc, y1 = xax_xloc + 0.02, xpd = NA, lwd = line_lwd)
segments(x0 = x_base$li - xax_locs, x1 = x_base$li - xax_locs, 
         y0 = xax_xloc, y1 = xax_xloc + 0.02, xpd = NA, lwd = line_lwd)
text(x = xax_locs + x_base$ri, y = xax_xloc, labels = xax_labs, pos = 3, xpd = NA,
     cex = 0.75)
text(x = x_base$li - xax_locs, y = xax_xloc, labels = xax_labs, pos = 3, xpd = NA,
     cex = 0.75)

#outer
segments(x0 = x_base$lo, x1 = x_base$lo - width_peaks, 
         y0 = xax_xloc, y1 = xax_xloc, col = 1, xpd = NA, lwd = line_lwd)
segments(x0 = x_base$ro, x1 = x_base$ro + width_peaks, 
         y0 = xax_xloc, y1 = xax_xloc, col = 1, xpd = NA, lwd = line_lwd)
segments(x0 = xax_locs + x_base$ro, x1 = xax_locs + x_base$ro, 
         y0 = xax_xloc, y1 = xax_xloc + 0.02, xpd = NA, lwd = line_lwd)
segments(x0 = x_base$lo - xax_locs, x1 = x_base$lo - xax_locs, 
         y0 = xax_xloc, y1 = xax_xloc + 0.02, xpd = NA, lwd = line_lwd)
text(x = xax_locs + x_base$ro, y = xax_xloc, labels = xax_labs, pos = 3, xpd = NA,
     cex = 0.75)
text(x = x_base$lo - xax_locs, y = xax_xloc, labels = xax_labs, pos = 3, xpd = NA,
     cex = 0.75)

#label horizontal axis
hax_label <- c(q = "x-pop quantiles", z = "z-scores")[values_to_use]
hax_label_width <- strwidth(hax_label) + xdisp * ifelse(MFplot, 2, 1)
hax_label_y <- 1.1
text(x = midp_x, y = hax_label_y, labels = hax_label, xpd = NA, font = 2)
bezier_curve(y1 = xax_xloc, x0 = midp_x - hax_label_width / 2, 
             y0 = hax_label_y, x1 = x_base$li, 
             xpd = NA, 
             ends = "flat")
bezier_curve(y1 = xax_xloc, x0 = midp_x + hax_label_width / 2, 
             y0 = hax_label_y, x1 = x_base$ri, 
             xpd = NA, 
             ends = "flat")

# text(xax_xloc-0.17, y = x_base$li - 1, labels = "reciprocal quantile",
#      srt = 90, xpd = NA, cex = 1.25)

#draw vertical guiding lines
y_bot   <- -0.05
segments(x0 = c(xax_locs + x_base$ri, x_base$li - xax_locs),
         x1 = c(xax_locs + x_base$ri, x_base$li - xax_locs),
         y0 = y_bot, y1 = 1.04, xpd = NA, 
         lwd = 0.75, col = adjustcolor(1, 0.15), lty = 2)

segments(x0 = c(xax_locs + x_base$ro, x_base$lo - xax_locs),
         x1 = c(xax_locs + x_base$ro, x_base$lo - xax_locs),
         y0 = y_bot, y1 = 1.04, xpd = NA, 
         lwd = 0.75, col = adjustcolor(1, 0.15), lty = 2)

#draw lines and peaks corresponding to traits

# compact polygon mountains for li, lo, ri, ro using quantiles
reg_order <- c("li","lo","ri","ro")
extra_y_poly <- 0.01

#get maximum height for ev df
# evs: list of equal-shaped data frames with numeric columns
ev_max <- as.data.frame(do.call(pmax,
                             c(lapply(evs, function(x) as.matrix(x[,-1])), list(na.rm = TRUE))))
ev_max <- cbind(variable = evs[[1]]$variable, ev_max)

for(focal_sex in focal_sexes){
  
  #subset data to the appropriate subset
  ev <- evs[[focal_sex]]
  
  for(reg in reg_order){
    
    # pick population pair and baseline
    pair <- pops[[reg]]                 # e.g., c("HAEE","HARE")
    x0   <- x_base[[reg]]               # baseline x for this region
    dir  <- if(substr(reg,1,1) == "l") -1 else 1  # left = -1, right = +1
    
    # choose column from ev for quantiles (z-score path left as a stub)
    if(values_to_use == "q"){
      colnm <- paste0("q_", pair[1], "_in_", pair[2])
      xvals <- setNames(ev[[colnm]], ev$variable)[names(labs_ord)]
      xvals_max <- setNames(ev_max[[colnm]], ev_max$variable)[names(labs_ord)]
    } else {
      colnm <- paste0("z_", pair[1], "_in_", pair[2])
      zraw  <- setNames(ev[[colnm]], ev$variable)[names(labs_ord)]
      xvals <- z_to_unit(zraw)
      zraw_max  <- setNames(ev_max[[colnm]], ev_max$variable)[names(labs_ord)]
      xvals_max <- z_to_unit(zraw_max)
    }
    
    # mask missing, compute curve x at each vertical y
    miss <- is.na(xvals)
    yy <- y_plot[!miss]
    xx <- x0 + dir * width_peaks * xvals[!miss]
    xx_max <- x0 + dir * width_peaks * xvals_max[!miss]
    
    # dashed guides
    segments(x0 = x0, x1 = xx, y0 = yy, y1 = yy, lty = 3, col = adjustcolor(1, 0.5))
    
    # polygon path: go up along the curve in y, then back down along the baseline
    o <- order(yy); yy <- yy[o]; xx <- xx[o]
    ycap <- c(min(yy) - extra_y_poly, yy, max(yy) + extra_y_poly)
    xcap <- c(x0, xx, x0)
    
    #connect inner to outer polys
    # if(reg %in% c("li","ri")){
    #   outer_reg <- if(reg == "li") "lo" else "ro"
    #   segments(x0 = xx, x1 = x_base[[outer_reg]], y0 = yy, y1 = yy, lwd = 1, lty = 1, col = 1)
    # }
    
    # connect inner to outer polys (even if inner is NA, when outer exists)
    if(reg %in% c("li","ri")){
      outer_reg <- if(reg == "li") "lo" else "ro"
      
      # OUTER values (aligned with current scaling)
      pair_out <- pops[[outer_reg]]
      if(values_to_use == "q"){
        col_out   <- paste0("q_", pair_out[1], "_in_", pair_out[2])
        xvals_out <- setNames(ev[[col_out]], ev$variable)[names(labs_ord)]
      } else {
        col_out   <- paste0("z_", pair_out[1], "_in_", pair_out[2])
        zraw_out  <- setNames(ev[[col_out]], ev$variable)[names(labs_ord)]
        xvals_out <- z_to_unit(zraw_out)  # from your z-setup
      }
      miss_out <- is.na(xvals_out)
      
      # target y's are those present in OUTER layer
      yy_tar <- y_plot[!miss_out]
      
      # we already have inner curve sorted as yy (y) and xx (x) above
      if(length(yy) >= 2 && length(yy_tar)){
        # interpolate x_on_inner at each yy_tar along the inner curve
        idx <- findInterval(yy_tar, yy)             # yy[idx] <= yy_tar < yy[idx+1]
        x_on <- rep(NA_real_, length(yy_tar))
        
        # exact vertex matches
        eq <- match(yy_tar, yy, nomatch = 0)
        if(any(eq > 0)) x_on[eq > 0] <- xx_max[eq[eq > 0]]
        
        # linear interpolation for interior points
        ok <- idx > 0 & idx < length(yy) & is.na(x_on)
        if(any(ok)){
          t <- (yy_tar[ok] - yy[idx[ok]]) / (yy[idx[ok] + 1] - yy[idx[ok]])
          x_on[ok] <- xx_max[idx[ok]] + t * (xx_max[idx[ok] + 1] - xx_max[idx[ok]])
        }
        
        # draw connectors for those with a valid intersection
        ok2 <- !is.na(x_on)
        if(any(ok2)){
          segments(x0 = x_on[ok2], x1 = x_base[[outer_reg]],
                   y0 = yy_tar[ok2], y1 = yy_tar[ok2],
                   lwd = 1, lty = 1, col = 1)
        }
      }
    }
    
    
    #label peaks with true heights
    if(MFplot){
      if(reg %in% c("lo","ro") & focal_sex == focal_sexes[[1]]){
        
        split_bezier <- F
        
        if(split_bezier){
          vals_txt_list <- lapply(focal_sexes, function(fs){
            raw_vals <- setNames(evs[[fs]][[paste0(pair[1], "_mu")]], 
                                 evs[[fs]]$variable)[names(labs_ord)]
            raw_vals <- round(raw_vals, 2)
            paste0(c("Male" = "", "Female" = "\n")[fs], 
                   sapply(raw_vals[!miss], pretty_large_number, nsf = 1), " ",
                   unit_names[names(labs_ord)][!miss])  
          })
          vals_txt <- Reduce(paste0, vals_txt_list)
          
          # place labels slightly "above" in y, and a bit farther out in x
          yy_lab <- laby_peaks[!miss]
          xx_curve <- x0 + dir * width_peaks * xvals_max[!miss]
          xx_lab <- x0 + dir * (width_peaks * pmax(xvals_max[!miss], 0.55) + 
                                  elbow_l*4)
          
          for(i in seq_along(yy)){
            extra_x <- elbow_l * 4
            # curved leader from curve to label midpoint
            bezier_curve(x0 = xx_curve[i], y0 = y_plot[!miss][i],
                         x1 = xx_lab[i],   y1 = yy_lab[i],
                         xpd = NA, ends = "flat")
            
            #further lines to sex-specific labels
            bezier_curve(x1 = xx_lab[i] - (elbow_l - extra_x) * dir, 
                         y1 = yy_lab[i] + minsep / 2,
                         x0 = xx_lab[i], y0 = yy_lab[i],
                         xpd = NA, ends = "flat")
            bezier_curve(x1 = xx_lab[i] - (elbow_l - extra_x) * dir, 
                         y1 = yy_lab[i] - minsep / 4,
                         x0 = xx_lab[i], y0 = yy_lab[i],
                         xpd = NA, ends = "flat")
            
            # text label
            xlab <- xx_lab[i] - (elbow_l - extra_x + xdisp) * dir
            ylab <- yy_lab[i]
            dim_lab <- c(w = strwidth(vals_txt[i], cex = peak_lab_cex),
                         h = strheight(vals_txt[i], cex = peak_lab_cex))
            #main text label
            text(x = xlab, y = ylab,
                 labels = vals_txt[i], xpd = NA,
                 pos = if(dir == -1) 2 else 4, cex = peak_lab_cex)
            
            #sex identifiers
            text(x = xlab + dim_lab["w"] * dir, y = ylab + minsep / 2,
                 labels = " ( ♂ ) ", xpd = NA,
                 pos = if(dir == -1) 2 else 4, cex = peak_lab_cex,
                 col = sex_cols["Male"])
            text(x = xlab + dim_lab["w"] * dir, y = ylab - minsep / 4,
                 labels = " ( ♀ ) ", xpd = NA,
                 pos = if(dir == -1) 2 else 4, cex = peak_lab_cex,
                 col = sex_cols["Female"]) 
          }
          
        } else {
          # get male and female raw values for this population (pair[1]) on the same variables
          raw_vals_sex <- lapply(c("Male", "Female"), function(fs) {
            setNames(evs[[fs]][[paste0(pair[1], "_mu")]],
                     evs[[fs]]$variable)[names(labs_ord)]
          })
          names(raw_vals_sex) <- c("Male", "Female")
          raw_vals_sex <- lapply(raw_vals_sex, function(v) round(v, 2))
          
          # pretty-print male & female values, masking by current miss pattern
          vals_M    <- sapply(raw_vals_sex[["Male"]][!miss],   pretty_large_number, nsf = 1)
          vals_F    <- sapply(raw_vals_sex[["Female"]][!miss], pretty_large_number, nsf = 1)
          units_txt <- unit_names[names(labs_ord)][!miss]
          
          # label positions (parallel to your original peak-label code)
          yy_lab   <- laby_peaks[!miss]
          xx_curve <- x0 + dir * width_peaks * xvals_max[!miss]
          xx_lab   <- x0 + dir * (width_peaks * pmax(xvals_max[!miss], 0.55) + elbow_l * 4)
          
          for (i in seq_along(yy_lab)) {
            
            # single curved leader from curve to label
            bezier_curve(x0 = xx_curve[i], y0 = y_plot[!miss][i],
                         x1 = xx_lab[i],   y1 = yy_lab[i],
                         xpd = NA, ends = "flat")
            
            ylab    <- yy_lab[i]
            cex_lab <- peak_lab_cex
            
            # pieces of the label: "M, F unit"
            lab_M     <- vals_M[i]
            lab_F     <- vals_F[i]
            lab_comma <- ", "
            lab_unit  <- paste0(" ", units_txt[i])
            
            # widths at this cex
            w_space     <- strwidth("  ",     cex = cex_lab)
            w_M     <- strwidth(lab_M,     cex = cex_lab)
            w_comma <- strwidth(lab_comma, cex = cex_lab)
            w_F     <- strwidth(lab_F,     cex = cex_lab)
            w_unit  <- strwidth(lab_unit,  cex = cex_lab)
            w_total <- w_M + w_comma + w_F + w_unit
            
            # this is the SAME anchor as your non-sex-split version:
            #   text(x = xx_lab[i] - elbow_l * dir, ..., pos = if(dir == -1) 2 else 4)
            #
            # for dir ==  1 (right side), pos = 4 means "to the right" of x_anchor
            # for dir == -1 (left side), pos = 2 means "to the left"  of x_anchor
            x_anchor <- xx_lab[i] - elbow_l * dir
            
            if (dir == 1) {
              # right side: label runs from x_anchor (left edge) to x_anchor + w_total (right edge)
              x_M     <- x_anchor + w_space
              x_comma <- x_M     + w_M
              x_F     <- x_comma + w_comma
              x_unit  <- x_F     + w_F
              
              text(x = x_M, y = ylab, labels = lab_M,
                   xpd = NA, adj = c(0, 0.5), cex = cex_lab,
                   col = sex_cols["Male"])
              text(x = x_comma, y = ylab, labels = lab_comma,
                   xpd = NA, adj = c(0, 0.5), cex = cex_lab)
              text(x = x_F, y = ylab, labels = lab_F,
                   xpd = NA, adj = c(0, 0.5), cex = cex_lab,
                   col = sex_cols["Female"])
              text(x = x_unit, y = ylab, labels = lab_unit,
                   xpd = NA, adj = c(0, 0.5), cex = cex_lab)
              
            } else {
              # left side (dir == -1): label runs from x_anchor - w_total (left edge) to x_anchor (right edge)
              x_unit  <- x_anchor - w_space
              x_F     <- x_unit  - w_unit
              x_comma <- x_F     - w_F
              x_M     <- x_comma - w_comma
              
              text(x = x_M, y = ylab, labels = lab_M,
                   xpd = NA, adj = c(1, 0.5), cex = cex_lab,
                   col = sex_cols["Male"])
              text(x = x_comma, y = ylab, labels = lab_comma,
                   xpd = NA, adj = c(1, 0.5), cex = cex_lab)
              text(x = x_F, y = ylab, labels = lab_F,
                   xpd = NA, adj = c(1, 0.5), cex = cex_lab,
                   col = sex_cols["Female"])
              text(x = x_unit, y = ylab, labels = lab_unit,
                   xpd = NA, adj = c(1, 0.5), cex = cex_lab)
            }
          } 
        }
      }
    } else {
      if(reg %in% c("lo","ro")){
        
        raw_vals <- setNames(ev[[paste0(pair[1], "_mu")]], ev$variable)[names(labs_ord)]
        raw_vals <- pretty_large_number(raw_vals, 2)
        vals_txt <- paste0(raw_vals[!miss], " ",
                           unit_names[names(labs_ord)][!miss])
        
        # place labels slightly "above" in y, and a bit farther out in x
        yy_lab <- laby_peaks[!miss]
        xx_curve <- x0 + dir * width_peaks * xvals[!miss]
        xx_lab <- x0 + dir * (width_peaks * pmax(xvals[!miss], 0.55) + elbow_l*4)
        
        for(i in seq_along(yy)){
          # curved leader from curve to label
          bezier_curve(x0 = xx_curve[i], y0 = y_plot[!miss][i],
                       x1 = xx_lab[i],   y1 = yy_lab[i],
                       xpd = NA, ends = "flat")
          # text label (vertical like before)
          text(x = xx_lab[i] - elbow_l * dir, y = yy_lab[i],
               labels = vals_txt[i], xpd = NA,
               pos = if(dir == -1) 2 else 4)
        }
      }
    }
    
    # 1. Create a new xvals vector for the grey background
    xvals_grey <- xvals
    # 'miss' is the original logical mask of missing values
    idx_miss <- which(miss)
    idx_present <- which(!miss)
    
    # Interpolate the missing values, if possible
    # We must have at least 2 non-missing points to interpolate between
    if (length(idx_miss) > 0 && length(idx_present) >= 2) {
      # We use y_plot as the 'x' for approx and xvals as the 'y'
      interp <- approx(x = y_plot[idx_present], 
                       y = xvals[idx_present], 
                       xout = y_plot[idx_miss], 
                       rule = 2) # rule=2 extrapolates to nearest value
      xvals_grey[idx_miss] <- interp$y
      
    } else if (length(idx_miss) > 0) {
      # We have NAs, but not enough data to interpolate.
      # Fall back to baseline (0) for missing values.
      xvals_grey[idx_miss] <- 0
    }
    
    # Define the FULL polygon shape using the interpolated 'xvals_grey'
    yy_full <- y_plot # Use all y-values
    # Note: xvals_grey now has no NAs
    xx_full <- x0 + dir * width_peaks * xvals_grey
    
    # Create the polygon path for the *entire* y-range
    ycap_full <- c(min(yy_full) - extra_y_poly, yy_full, max(yy_full) + extra_y_poly)
    xcap_full <- c(x0, xx_full, x0)
    xpoly_full <- c(xcap_full, rep(x0, length(xcap_full)))
    ypoly_full <- c(ycap_full, rev(ycap_full))
    
    # Draw this complete shape as the "missing" (grey/dashed) background
    # polygon(xpoly_full, ypoly_full, 
    #         col = col_miss_fill, 
    #         border = NA)
    
    # keep these in case other code below still expects them
    xpoly <- c(xcap, rep(x0, length(xcap)))
    ypoly <- c(ycap, rev(ycap))
    
    # full grids for this region
    yy_full <- y_plot
    xx_full <- x0 + dir * width_peaks * xvals   # may include NA for missing
    nfeat   <- length(yy_full)
    
    # colours & linetypes
    if (MFplot) {
      col_present_fill   <- sex_cols_fill[focal_sex]
      col_present_border <- sex_cols[focal_sex]
    } else {
      col_present_fill   <- cols_fill[[pair[1]]]
      col_present_border <- cols_border[[pair[1]]]
    }
    col_miss_fill   <- adjustcolor("grey70", alpha.f = 0.35)
    col_miss_border <- adjustcolor("grey40", alpha.f = 0.9)
    lty_miss  <- 2
    lwd_border <- 1
    
    # helper: contiguous runs of TRUE/FALSE
    .get_runs <- function(flag_vec){
      r <- rle(flag_vec)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1
      list(values = r$values, starts = starts, ends = ends)
    }
    
    ## 1) PRESENT blocks: coloured polys
    runs_p <- .get_runs(!miss)
    for (k in which(runs_p$values)) {
      
      idx <- runs_p$starts[k]:runs_p$ends[k]
      
      yy_run <- yy_full[idx]
      xx_run <- xx_full[idx]  # non-NA here
      o <- order(yy_run); yy_run <- yy_run[o]; xx_run <- xx_run[o]
      
      y_min <- min(yy_run)
      y_max <- max(yy_run)
      if (runs_p$starts[k] == 1)     y_min <- y_min - extra_y_poly
      if (runs_p$ends[k]   == nfeat) y_max <- y_max + extra_y_poly
      
      ycap_run <- c(y_min, yy_run, y_max)
      xcap_run <- c(x0,    xx_run,  x0)
      
      xpoly_run <- c(xcap_run, rep(x0, length(xcap_run)))
      ypoly_run <- c(ycap_run, rev(ycap_run))
      
      # if shrink_gray_poly_y_by > 0, draw full polygon border; otherwise border handled manually
      if (shrink_gray_poly_y_by > 0) {
        polygon(xpoly_run, ypoly_run,
                col = col_present_fill,
                border = col_present_border)
      } else {
        polygon(xpoly_run, ypoly_run,
                col = col_present_fill,
                border = NA)
      }
    }
    
    ## 2) MISSING blocks: grey polys, shrunk vertically by `shrink_gray_poly_y_by`
    runs_m <- .get_runs(miss)
    for (k in which(runs_m$values)) {
      
      idx_miss  <- runs_m$starts[k]:runs_m$ends[k]
      i_first   <- min(idx_miss)
      i_last    <- max(idx_miss)
      
      prev_i <- i_first - 1
      next_i <- i_last  + 1
      
      has_prev <- prev_i >= 1     && !miss[prev_i]
      has_next <- next_i <= nfeat && !miss[next_i]
      
      # baseline caps: start from neighbours (if present) or extend outwards
      if (has_prev) {
        y_low_cap <- yy_full[prev_i]
      } else {
        y_low_cap <- yy_full[i_first] - extra_y_poly
      }
      if (has_next) {
        y_high_cap <- yy_full[next_i]
      } else {
        y_high_cap <- yy_full[i_last] + extra_y_poly
      }
      
      # apply vertical shrink_gray_poly_y_by to grey region (top and bottom)
      if (shrink_gray_poly_y_by > 0) {
        y_low_cap  <- y_low_cap  + shrink_gray_poly_y_by
        y_high_cap <- y_high_cap - shrink_gray_poly_y_by
        if (y_low_cap >= y_high_cap) {
          next
        }
      }
      
      nb_idx <- integer(0)
      if (has_prev) nb_idx <- c(nb_idx, prev_i)
      if (has_next) nb_idx <- c(nb_idx, next_i)
      
      if (length(nb_idx) == 2) {
        # two neighbours: build a quadrilateral, ordered to avoid crossing
        ord <- order(yy_full[nb_idx])
        nb_low  <- nb_idx[ord[1]]
        nb_high <- nb_idx[ord[2]]
        
        # path (clockwise, no X):
        # baseline low -> baseline high -> curve high -> curve low
        # original neighbour coordinates (sorted so nb_low has smaller y)
        x_low  <- xx_full[nb_low]
        y_low  <- yy_full[nb_low]
        x_high <- xx_full[nb_high]
        y_high <- yy_full[nb_high]
        
        # target y-values on the curve after shrinking
        y_curve_low  <- y_low  + shrink_gray_poly_y_by
        y_curve_high <- y_high - shrink_gray_poly_y_by
        
        # if shrink is extreme, skip drawing this gap
        if (y_curve_low >= y_curve_high) {
          next
        }
        
        # interpolate the corresponding x-values along the line between neighbours
        x_curve_low  <- interp_x_for_y(y_curve_low,  x_low, y_low, x_high, y_high)
        x_curve_high <- interp_x_for_y(y_curve_high, x_low, y_low, x_high, y_high)
        
        # polygon: baseline low -> baseline high -> curve high -> curve low
        xpoly_m <- c(x0,
                     x0,
                     x_curve_high,
                     x_curve_low)
        ypoly_m <- c(y_low_cap,
                     y_high_cap,
                     y_curve_high,
                     y_curve_low)
        
        if (shrink_gray_poly_y_by > 0) {
          polygon(xpoly_m, ypoly_m,
                  col = col_miss_fill,
                  border = col_miss_border,
                  lty = lty_miss)
        } else {
          polygon(xpoly_m, ypoly_m,
                  col = col_miss_fill,
                  border = NA)
        }
        
      } else if (length(nb_idx) == 1) {
        # only one neighbour (tail region): triangle anchored on that neighbour
        nb <- nb_idx[1]
        mid_y <- (y_low_cap + y_high_cap) / 2
        
        if (yy_full[nb] >= mid_y) {
          # neighbour above midpoint
          xpoly_m <- c(x0, x0, xx_full[nb])
          ypoly_m <- c(y_low_cap, y_high_cap, yy_full[nb])
        } else {
          # neighbour below midpoint
          xpoly_m <- c(x0, xx_full[nb], x0)
          ypoly_m <- c(y_low_cap, yy_full[nb], y_high_cap)
        }
        
        if (shrink_gray_poly_y_by > 0) {
          polygon(xpoly_m, ypoly_m,
                  col = col_miss_fill,
                  border = col_miss_border,
                  lty = lty_miss)
        } else {
          polygon(xpoly_m, ypoly_m,
                  col = col_miss_fill,
                  border = NA)
        }
        
      } else {
        # no neighbours (everything missing for this region): nothing to fill
      }
    }
    
    ## 3) BORDERS on top: only when shrink_gray_poly_y_by == 0 (outer borders only)
    if (shrink_gray_poly_y_by == 0) {
      
      # present borders
      for (k in which(runs_p$values)) {
        idx <- runs_p$starts[k]:runs_p$ends[k]
        
        # curve edge
        lines(xx_full[idx], yy_full[idx],
              col = col_present_border,
              lty = 1, lwd = lwd_border)
        
        # baseline segment for this block
        segments(x0 = x0,
                 y0 = min(yy_full[idx]),
                 x1 = x0,
                 y1 = max(yy_full[idx]),
                 col = col_present_border,
                 lty = 1, lwd = lwd_border)
        
        # caps at global extremes
        if (runs_p$starts[k] == 1) {
          segments(x0 = x0,
                   y0 = min(yy_full[idx]) - extra_y_poly,
                   x1 = xx_full[idx][1],
                   y1 = yy_full[idx][1],
                   col = col_present_border,
                   lty = 1, lwd = lwd_border)
        }
        if (runs_p$ends[k] == nfeat) {
          segments(x0 = x0,
                   y0 = max(yy_full[idx]) + extra_y_poly,
                   x1 = xx_full[idx][length(idx)],
                   y1 = yy_full[idx][length(idx)],
                   col = col_present_border,
                   lty = 1, lwd = lwd_border)
        }
      }
      
      # missing borders: dashed, using neighbour indices
      for (k in which(runs_m$values)) {
        
        idx_miss  <- runs_m$starts[k]:runs_m$ends[k]
        i_first   <- min(idx_miss)
        i_last    <- max(idx_miss)
        
        prev_i <- i_first - 1
        next_i <- i_last  + 1
        
        # extend by -1,+1 to go back to observed traits
        idx_b <- c(prev_i, idx_miss, next_i)
        idx_b <- idx_b[idx_b >= 1 & idx_b <= nfeat]
        
        # curve side: only use positions where xx_full is finite
        ok_curve <- is.finite(xx_full[idx_b]) & is.finite(yy_full[idx_b])
        if (sum(ok_curve) >= 2) {
          lines(xx_full[idx_b][ok_curve],
                yy_full[idx_b][ok_curve],
                col = col_miss_border,
                lty = lty_miss,
                lwd = lwd_border)
        }
        
        # baseline side over the same y-span
        y_span <- range(yy_full[idx_b])
        segments(x0 = x0,
                 y0 = y_span[1],
                 x1 = x0,
                 y1 = y_span[2],
                 col = col_miss_border,
                 lty = lty_miss,
                 lwd = lwd_border)
      }
    }
    
    ## 4) MIDLINE: only through present blocks (not grey gaps)
    mid_val <- 0.5
    x_mid <- x0 + dir * width_peaks * mid_val
    if (!MFplot) {
      
      # loop over present runs and intersect midline with each present polygon
      for (k in which(runs_p$values)) {
        
        idx <- runs_p$starts[k]:runs_p$ends[k]
        yy_run <- yy_full[idx]
        xx_run <- xx_full[idx]
        o <- order(yy_run); yy_run <- yy_run[o]; xx_run <- xx_run[o]
        
        y_min <- min(yy_run)
        y_max <- max(yy_run)
        if (runs_p$starts[k] == 1)     y_min <- y_min - extra_y_poly
        if (runs_p$ends[k]   == nfeat) y_max <- y_max + extra_y_poly
        
        ycap_run <- c(y_min, yy_run, y_max)
        xcap_run <- c(x0,    xx_run,  x0)
        
        xpoly_run <- c(xcap_run, rep(x0, length(xcap_run)))
        ypoly_run <- c(ycap_run, rev(ycap_run))
        
        n_run <- length(xpoly_run)
        if (n_run >= 3) {
          yints <- numeric(0)
          for (j in seq_len(n_run)) {
            i1 <- j
            i2 <- if (j < n_run) j + 1 else 1
            xi <- xpoly_run[i1]; xj <- xpoly_run[i2]
            yi <- ypoly_run[i1]; yj <- ypoly_run[i2]
            
            if ((xi <= x_mid && x_mid < xj) ||
                (xj <= x_mid && x_mid < xi)) {
              t <- (x_mid - xi) / (xj - xi)
              yints <- c(yints, yi + t * (yj - yi))
            }
          }
          
          if (length(yints) >= 2) {
            yints <- sort(yints)
            for (m in seq(1, length(yints) - 1, by = 2)) {
              segments(x0 = x_mid, y0 = yints[m],
                       x1 = x_mid, y1 = yints[m + 1],
                       lwd = 1, lty = 1, col = 1, xpd = NA)
            }
          }
        }
      }
    }
    
    # guiding arrows + labels at bottom for this region
    if(focal_sex == focal_sexes[[1]]){
      # arrow start offset from midline and arrow half-span
      gap   <- 0.15 * width_peaks              # how far from midline arrows begin
      span  <- (width_peaks / 2) - gap         # arrow length from its start to tip
      
      # left arrow: starts left of midline, points further left
      xL0 <- x_mid - gap
      xL1 <- x_mid - (gap + span)
      
      # right arrow: starts right of midline, points further right
      xR0 <- x_mid + gap
      xR1 <- x_mid + (gap + span)
      
      # place labels centered under each arrow, colored by group
      lab_left  <- if (dir == -1) pair[1] else pair[2]
      lab_right <- if (dir == +1) pair[1] else pair[2]
      xL_mid <- (xL0 + xL1) / 2
      xR_mid <- (xR0 + xR1) / 2
      text(x = xL_mid, y = y_bot - ydisp * 0.8, labels = paste0(lab_left, "\nlarger"),
           xpd = NA, cex = 0.8, col = cols_fill[[lab_left]], font = 2)
      text(x = xR_mid, y = y_bot - ydisp * 0.8, labels = paste0(lab_right, "\nlarger"),
           xpd = NA, cex = 0.8, col = cols_fill[[lab_right]], font = 2)  
      
      # draw arrows (unconnected)
      arrows(x0 = xL0, y0 = y_bot, x1 = xL1, y1 = y_bot, lwd = 2, lty = 1, 
             col = cols_border[lab_left], xpd = NA, length = 0.05)
      arrows(x0 = xR0, y0 = y_bot, x1 = xR1, y1 = y_bot, lwd = 2, lty = 1,
             col = cols_border[lab_right], xpd = NA, length = 0.05)
      
      #connect arrows with bezier curve
      arrow_rise <- 0.02
      bezier_curve(x0 = xL0, y0 = y_bot, p=0.5, k=1,
                   x1 = x_mid,   y1 = y_bot + arrow_rise,
                   xpd = NA, ends = "flat",
                   lwd = 2, lty = 1,
                   col = cols_border[lab_left])
      bezier_curve(x0 = xR0, y0 = y_bot,  p=0.5, k=1,
                   x1 = x_mid,   y1 = y_bot + arrow_rise,
                   xpd = NA, ends = "flat",
                   lwd = 2, lty = 1,
                   col = cols_border[lab_right]) 
    }
    
  }
  
}

#draw in titles
text(x = mean(x_base$li - width_peaks, x_base$lo), 
     y = 1.1375, labels = panel_titles[pops$li[1]], pos = 3,
     cex = 1.25, font = 2, col = cols_fill[pops$li[1]], xpd = NA)
text(x = c(x_base$li - width_peaks / 2, x_base$lo - width_peaks / 2), 
     y = 1.0825, labels = paste0("vs. ", c(pops$li[2], pops$lo[2])), pos = 3,
     cex = 1, font = 2, col = cols_fill[c(pops$li[2], pops$lo[2])], xpd = NA)

text(x = mean(x_base$ri + width_peaks, x_base$ro), 
     y = 1.1375, labels = panel_titles[pops$ri[1]], pos = 3,
     cex = 1.25, font = 2, col = cols_fill[pops$ri[1]], xpd = NA)
text(x = c(x_base$ri + width_peaks / 2, x_base$ro + width_peaks / 2), 
     y = 1.0825, labels = paste0("vs. ", c(pops$ri[2], pops$ro[2])), pos = 3,
     cex = 1, font = 2, col = cols_fill[c(pops$ri[2], pops$ro[2])], xpd = NA)

#draw in note for the mean values at the bottom
note_txt <- lapply(setNames(focal_sexes, focal_sexes), function(focal_sex){
  mv <- mean_vals_list[[focal_sex]][1, ]
  paste0(
    ifelse(focal_sex == focal_sexes[1], "note: expectations computed ", "and "), 
    ifelse(focal_sex == "Both", "at the sex midpoint ", paste0("in ", tolower(focal_sex), "s ")),
    "for an individual aged ", round(mv$age_years, 0), " y, weighing ",
    round(mv$body_weight, 0), " kg, and standing ", round(mv$body_height_cm, 0), " cm tall"
  )
})
note_txt <- paste0(note_txt, collapse = ifelse(MFplot, ",\n", ""))

# convert device corners (ndc = 0..1) to user coords, then add small device-relative margins
x_dev_left   <- grconvertX(0, from = "ndc", to = "user")
y_dev_bottom <- grconvertY(0, from = "ndc", to = "user")
dx <- grconvertX(0.005, from = "ndc", to = "user") - grconvertX(0, from = "ndc", to = "user")  # 2% of device width
dy <- grconvertY(0.01, from = "ndc", to = "user") - grconvertY(0, from = "ndc", to = "user")  # 2% of device height

x_note <- x_dev_left + dx
y_note <- y_dev_bottom + dy

text(x = x_note, y = y_note, labels = note_txt,
     adj = c(0, 0), xpd = NA, cex = 0.7, col = adjustcolor(1, 0.65))

#add a small legend for sex
if(MFplot){
  legend(
    x = 2.4, y = -0.1,
    legend = names(sex_cols),
    horiz = TRUE,            # horizontal items
    pch = 22,                # filled square point
    col = sex_cols,          # symbol border + text color
    pt.bg = sex_cols_fill,   # symbol fill color
    cex = 1,
    xpd = NA,
    bty = "n",
    xjust = 0.5, 
    yjust = 0.5,
    pt.cex = 2,
    text.width = 0.35
    
  )
  
}

dev.off()

#convert to png too
magick::image_write(magick::image_read_pdf(fig_path, density = 300), 
                    path = gsub("\\.pdf$", ".png", fig_path), format = "png")


#### quick viz for talk ####


#### 1 trait hists ####
variables_to_inspect
traits <- c("chest_press_1rm_lbs", "vo2_peak_mlkg")

for(trait in traits){
  
  fig_path <- paste0(pres_fig_dir, trait, "-comparison.pdf")
  
  #histogram of interesting variable
  cairo_pdf(fig_path,
            width = 500/72, height = 400/72, pointsize = 16)
  par(mfrow = c(1, 1), oma = c(0,0,0,0), mar = c(4,4,4,4))
  
  trait_dat <- list("HA" = c(re_ph[,trait], ee_ph[,trait]),
                    "SED" = sd_ph[,trait])
  trait_dat$SED <- trait_dat$SED[!is.na(trait_dat$SED)]
  trait_dat$HA <- trait_dat$HA[!is.na(trait_dat$HA)]
  multihist(trait_dat, main = trait, nbreaks = 20, xlab = unit_names[trait], 
            legend_location = "topright", leg_loc_x_extra = -0.1, horiz_leg = F)
  
  dev.off()
  
  magick::image_write(magick::image_read_pdf(fig_path, density = 300), 
                      path = gsub("\\.pdf$", ".png", fig_path), format = "png")
  
}


#### 2 conditional seq ####
traits <- c("chest_press_1rm_lbs", "vo2_peak_mlkg")
trait <- traits[2]
plot_cols <- setNames(adjustcolor(rainbow(n = 2)[c(1,2)], alpha.f = 0.5), 
                      c("HA", "SED"))

# parameters
ha_mu <- ev$HARE_mu[ev$variable == trait]
ha_sd <- ev$HARE_sd[ev$variable == trait]
sed_mu <- ev$SED_mu[ev$variable == trait]
sed_sd <- ev$SED_sd[ev$variable == trait]

n <- 512
x_min <- min(ha_mu - 4 * ha_sd, sed_mu - 4 * sed_sd)
x_max <- max(ha_mu + 4 * ha_sd, sed_mu + 4 * sed_sd)
x <- seq(x_min, x_max, length.out = n)
y_ha <- dnorm(x, mean = ha_mu, sd = ha_sd)
y_sed <- dnorm(x, mean = sed_mu, sd = sed_sd)
ylim <- c(0, max(y_ha, y_sed) * 1.05)
fig_paths <- c()
for (i in 1:5) {
  
  fig_path <- paste0(pres_fig_dir, trait, "-cc_", i, ".pdf")
  fig_paths <- c(fig_paths, fig_path)
  
  cairo_pdf(fig_path, width = 650/72, height = 300/72, pointsize = 16)
  par(mar = c(4.5, 6, 1, 1))
  plot(x, y_ha, type = "n", 
       xlab = "", 
       ylab = "density", bty = "l", ylim = ylim)
  mtext(paste0(trait, " | ", 
               paste0(c("activity_group", variables_to_condition_on), 
                      collapse = ", ")), side = 1, line = 2.5, cex = 0.75, font = 2)
  
  # always draw SED base polygon
  polygon(c(x, rev(x)), c(y_sed, rep(0, length(y_sed))),
          col = plot_cols["SED"], border = 1)
  
  # draw HA polygon in all panels except 1 and 4
  if (i %in% 2:3) {
    polygon(c(x, rev(x)), c(y_ha, rep(0, length(y_ha))),
            col = plot_cols["HA"], border = 1)
  }
  
  # add vertical line and label at HA mu for panels 3, 4, 5
  if (i >= 3) {
    segments(x0 = ha_mu, x1 = ha_mu, y0 = par("usr")[3], y1 = ylim[2] * 0.95, lty = 2, lwd = 2)
    text(x = ha_mu, y = ylim[2] * 0.93, pos = 3,
         labels = expression(mu[HA] * " | " * Predictors), xpd = NA)
  }
  
  # split-shade SED relative to HA mu:
  # - panel 4: only show the right (HA-colored) side; left is fully transparent
  # - panel 5: show both sides, then annotate % left/right of HA mu under SED
  if (i == 4 || i == 5) {
    ix_left <- x <= ha_mu
    left_col  <- adjustcolor(1, alpha.f = 0.0)
    right_col <- plot_cols["HA"]
    
    if (any(ix_left)) {
      polygon(c(x[ix_left], rev(x[ix_left])),
              c(y_sed[ix_left], rep(0, sum(ix_left))),
              col = left_col, border = NA)
    }
    if (any(!ix_left)) {
      polygon(c(x[!ix_left], rev(x[!ix_left])),
              c(y_sed[!ix_left], rep(0, sum(!ix_left))),
              col = right_col, border = NA)
    }
    
    # panel 5: annotate percentages over midpoints of left/right SED polys
    if (i == 5) {
      p_left  <- pnorm(ha_mu, mean = sed_mu, sd = sed_sd)
      p_right <- 1 - p_left
      
      # midpoint x-positions for labels
      x_left_mid  <- mean(range(x[ix_left]))
      x_right_mid <- mean(range(x[!ix_left]))
      # y-positions slightly above the curve at those x's
      y_left_mid  <- approx(x, y_sed, xout = x_left_mid)$y
      y_right_mid <- approx(x, y_sed, xout = x_right_mid)$y
      
      text(x_left_mid,  min(ylim[2], y_left_mid  + 0.08 * diff(ylim)),
           labels = paste0(round(100 * p_left, 1), "%"), cex = 2, pos = 2, col = adjustcolor(1, 0.85))
      text(x_right_mid, min(ylim[2], y_right_mid + 0.08 * diff(ylim)),
           labels = paste0(round(100 * p_right, 1), "%"), cex = 2, col = adjustcolor(1, 0.85))
    }
  }
  
  # legend:
  # - panel 1: only SED
  # - panels 2–5: both HA and SED
  if (i == 1) {
    legend("topright", legend = "SED",
           fill = plot_cols["SED"], border = NA, bty = "n")
  } else {
    legend("topright", legend = c("SED", "HA"),
           fill = plot_cols[c("SED", "HA")], border = NA, bty = "n")
  }
  
  dev.off()
  
  magick::image_write(
    magick::image_read_pdf(fig_path, density = 300),
    path = gsub("\\.pdf$", ".png", fig_path),
    format = "png"
  )
}


all_fig_path <- paste0(pres_fig_dir, trait, "-cc_all.pdf")
pdftools::pdf_combine(input = fig_paths, output = all_fig_path)

#### 3 heatmap ####
source("~/scripts/minor_scripts/postdoc/my_heatmap_ridgeline.R")

# mat_cols <- list(cmat1 = cmat, cmat2 = cmat, cmat3 = cmat) #for debugging multi-heatmap function
# nmat <- length(mat_cols)
# cairo_pdf(fig_path,
#           width = 400/72 * nmat, height = 400/72, pointsize = 16)
# par(mfrow = c(1, 1), oma = c(0,0,0,0), mar = c(3,12,5,2))

mat_cols <- cmat
for(i in 1:1){
  fig_path <- paste0(pres_fig_dir, "all-traits_correlation-heatmap_", i,".pdf")
  cairo_pdf(fig_path,
            width = 500/72, height = 400/72, pointsize = 16)
  par(mfrow = c(1, 1), oma = c(0,0,0,0), mar = c(2,10,5,1))
  
  #set index parameters
  hb <- i!=1
  b2h <- ifelse(i %in% 2:5, i-1, "all")
  b2h <- ifelse2(i == 6, c(1,4), b2h)
  b2h <- ifelse2(i == 7, c(1,3), b2h)
  
  my_heatmap(mat_cols = mat_cols, mat_size_rule = "abs", 
             legend_title = "pooled within-\ngroup correlation", leg_loc_scale = c(0.0, 0.1),
             cex.main = 0.7, plot_guiding_lines = T, use_bezier = F, 
             diag_matters = F, plot_diagonal_labels = T, 
             mds_method = c("olo", "smacof", "hclust")[3], 
             blockmodel = T, highlight_blocks = hb, blocks_to_highlight = b2h, 
             asp = 1, space_scale_ylabs = 1.3,
             plot_numbers = F)
  
  dev.off()
  
  magick::image_write(
    magick::image_read_pdf(fig_path, density = 300),
    path = gsub("\\.pdf$", ".png", fig_path),
    format = "png"
  )  
}


#### END ####

#can maybe also do a bargraph version of this with a 1D MDS?
#sorta like a mountain ridgeline (a viewpoint plot!)
#ha data is here: gs://motrpac-data-hub/human-main/phenotype/human-main-ha-adu
#can maybe also do the quantile of the residual distribution of HA <-> SED?