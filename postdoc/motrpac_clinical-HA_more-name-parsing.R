library(data.table)

#load sed data first
sd_ph <- as.data.frame(fread("~/data/motrpac_precovid_sedentary_clinical_obs.txt"))
cat(paste0("\"", colnames(sd_ph), "\","))

#look at curated data first
data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_sets/"
ha_paths <- list.files(data_dir)
colnms <- list()
for(i in seq_along(ha_paths)){
  ha_path <- paste0(data_dir, ha_paths[i])
  ha_ph <- as.data.frame(fread(ha_path))
  cat(paste0("\n\n", ha_path, ":\n"))
  cat(paste0("\"", colnames(ha_ph), "\","))
  colnms[[ha_path]] <- colnames(ha_ph)
}

data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/raw/data_sets/"
ha_paths <- list.files(data_dir)
colnms <- list()
for(i in seq_along(ha_paths)){
  ha_path <- paste0(data_dir, ha_paths[i])
  ha_ph <- as.data.frame(fread(ha_path))
  colnms[[ha_path]] <- colnames(ha_ph)
}
table(table(unlist(colnms)))
head(sort(table(unlist(colnms)), T), 12)

#look for filenames individually
setNames(unlist(colnms)[unlist(colnms) == "lean_mass"], NULL)
setNames(unlist(colnms)[grepl("lean", unlist(colnms))], NULL)

#### LLM guided filtration ####

# editable synonyms/keywords for your sd_ph columns
# note: plain strings are matched as substrings; anything with ".*" is a regex
manual_keywords <- list(
  pid                  = c("pid", "participant_id", "subject_id", "rec_id", "record_id", "id"),
  Timepoint            = c("timepoint", "time_point", "visit", "visit_code", "tp"),
  study_label          = c("study_label", "study", "arm", "form_name"),
  sex                  = c("sex", "gender", "biological_sex"),
  age_years            = c("age_years", "age_yrs", "age_year", "age", "age_at_visit"),
  visit_code           = c("visit_code", "visit", "visitnum", "visit_number"),
  body_weight          = c("body_weight", "weight", "wt", "wtkg", "weight_kg", "body_weight_kg"),
  body_height_cm       = c("body_height_cm", "height_cm", "height", "stature"),
  bmi                  = c("bmi", "body_mass_index", "index"),
  waist_circum_cm      = c("waist_circum_cm", "waist_circumference", "waist_cm", "waist", "wc"),
  bp_sys               = c("bp_sys", "sbp", "systolic_bp", "systolic", "dia"),
  bp_dia               = c("bp_dia", "dbp", "diastolic_bp", "diastolic", "sys"),
  pulse_rest           = c("pulse_rest", "resting_hr", "rest_hr", "hr_rest", "pulse"),
  vo2cart1             = c("vo2cart1", "vo2_cart_1", "vo2_stage1", "vo2"),
  vo2cart2             = c("vo2cart2", "vo2_cart_2", "vo2_stage2", "vo2"),
  watts_max            = c("watts_max", "wmax", "max_watts", "watts_peak", "peak_watts"),
  wtkg_cpet            = c("wtkg_cpet", "weight_cpet", "body_weight_cpet"),
  hr_end_cpet          = c("hr_end_cpet", "hr_peak_cpet", "heart_rate_end_cpet", "hr_end"),
  grip1                = c("grip", "handgrip1"),  # will be auto-expanded with trial patterns
  grip2                = c("grip", "handgrip2"),
  grip3                = c("grip", "handgrip3"),
  knee_iso_torque      = c("knee_iso_torque", "knee_isometric_torque", "knee_torque", "knee_ext_torque", "torque"),
  leg_press_1rm_lbs    = c("leg_press_1rm_lbs", "leg_press_1rm", "leg_press_1rm_kg", "lpmax"),
  leg_ext_1rm_lbs      = c("leg_ext_1rm_lbs", "leg_ext_1rm", "leg_extension_1rm", "leg_ext_1rm_kg", "lemax"),
  chest_press_1rm_lbs  = c("chest_press_1rm_lbs", "chest_press_1rm", "bench_press_1rm", "chest_press_1rm_kg", "cpmax"),
  lean_mass            = c("lean_mass", "total_lean_mass", "wbtot_lean", "lean"),
  fat_mass             = c("fat_mass", "total_fat_mass", "wbtot_fat", "fat"),
  total_bmd            = c("total_bmd", "whole_body_bmd", "wb_bmd", "wbtot_bmd", "bmd_total"),
  vat_mass             = c("vat_mass", "visceral_fat_mass", "vat"),
  vo2_peak_L           = c("vo2_peak_l", "vo2peak_l", "vo2_l", "vo2_absolute", "vo2max_l"),
  vo2_peak_mlkg        = c("vo2_peak_mlkg", "vo2peak_mlkg", "vo2_rel", "vo2_ml_kg_min", "vo2max_mlkg"),
  o2_pulse             = c("o2_pulse", "oxygen_pulse", "o2pulse"),
  leg_press_1rm_kg     = c("leg_press_1rm_kg", "leg_press_1rm"),
  leg_ext_1rm_kg       = c("leg_ext_1rm_kg", "leg_ext_1rm", "leg_extension_1rm"),
  chest_press_1rm_kg   = c("chest_press_1rm_kg", "chest_press_1rm", "bench_press_1rm"),
  handgrip_max         = c("handgrip_max", "grip_max", "max_grip", "handgrip_peak",
                           "grip", "handgrip", "grip_strength", "handgrip_strength")
)

# unresolved fields only
manual_keywords <- list(
  # body height (likely in hwwt)
  body_height_cm = c(
    "body_height_cm", "height_cm", "height", "heightcm",
    "stature", "stature_cm", "htcm", "ht_cm", "ht",
    "standing_height", "meas_height", "height.*cm", "ht.*cm"
  ),
  
  # cpet watts peak
  watts_max = c(
    "watts_max", "watt_max", "wmax", "wattpeak", "watt_peak",
    "peak_watts", "peak_power", "power_peak", "power_max",
    "workrate_peak", "workrate_max", "wrpeak", "pwrpeak",
    "wattend_cpet", "watt.*peak", "peak.*watt", "power.*peak", "work.*rate.*peak",
    "wrpeak", "wr_peak", "workrate_peak", "pwr_peak", "peak_workrate", "peak_wr"
  ),
  
  # vo2 absolute peak (L·min^-1)
  vo2_peak_L = c(
    "vo2_peak_l", "vo2peak_l", "peak_vo2", "vo2peak",
    "vo2max", "vo2_max", "vo2_abs", "vo2_l", "vo2_l_min", "vo2_lmin",
    "vo2peak_cpet", "vo2.*peak.*l", "peak.*vo2", "vo2.*max.*l",
    "vo2peak", "vo2_peak", "vo2max", "vo2_max", "vo2peak_l", "vo2_abs", "vo2_l", "vo2_l_min", "peak_vo2", "vo2peak_cpet"
  ),
  
  # vo2 relative peak (mL·kg^-1·min^-1)
  vo2_peak_mlkg = c(
    "vo2_peak_mlkg", "vo2peak_mlkg", "vo2_rel", "relative_vo2",
    "vo2_ml_kg_min", "vo2_mlkg_min", "vo2_mlkg", "vo2mlkg",
    "vo2max_mlkg", "peak_vo2_mlkg", "vo2.*ml.*kg", "ml.*kg.*vo2", "vo2_kg"
  ),
  
  # oxygen pulse (often derived)
  o2_pulse = c(
    "o2_pulse", "o2pulse", "oxygen_pulse", "o2_pulse_peak",
    "oxygen_pulse_peak", "o2pulse_cpet", "oxygen.*pulse"
  ),
  
  # 1RM in kg (search broadly; we can still derive from lbs if absent)
  leg_press_1rm_kg = c(
    "leg_press_1rm_kg", "leg_press_1rm", "lp_1rm_kg", "lp1rm_kg", "lp1rm",
    "leg_press_max_kg", "leg_press.*1rm.*kg", "kg.*1rm.*press.*leg",
    "lpmaxrest_fara", "lpmax.*"
  ),
  leg_ext_1rm_kg = c(
    "leg_ext_1rm_kg", "leg_ext_1rm", "leg_extension_1rm", "le_1rm_kg", "le1rm_kg", "le1rm",
    "leg_extension_max_kg", "leg.*ext.*1rm.*kg", "kg.*1rm.*ext.*leg",
    "leg_extension.*kg", "lemaxrest_fara", "lemax.*"
  ),
  chest_press_1rm_kg = c(
    "chest_press_1rm_kg", "chest_press_1rm", "bench_press_1rm",
    "cp1rm_kg", "bench_1rm_kg", "bench1rm", "chest.*press.*1rm.*kg",
    "kg.*1rm.*press.*chest", "bench.*1rm.*kg", "cpmaxrest_fara", "cpmax.*"
  ),
  
  # knee torque (iso / isokinetic)
  knee_iso_torque = c(
    "knee_iso_torque", "knee_isometric_torque", "knee_torque",
    "knee_ext_torque", "knee_flex_torque", "ext_torque", "flex_torque",
    "ext_peak_torque", "flex_peak_torque", "peak_torque", "ptq", "pt_peak",
    "nm", "newton_meter", "isometric", "iso", "isokinetic", "isok", "iske", "mvc",
    "knee.*torque", "torque.*knee", "knee.*iso", "iso.*knee", "isok.*knee", "ptq.*knee"
  )
)


# helper: make auto patterns like 'lean.*mass' to catch 'total_lean_mass'
make_auto_patterns <- function(x) {
  x <- tolower(x)
  tokens <- unlist(strsplit(gsub("[^a-z0-9]+", "_", x), "_"))
  tokens <- tokens[tokens != ""]
  pats <- character(0)
  if (length(tokens) >= 2) {
    # forward and reverse token-join patterns
    pats <- c(
      paste(tokens, collapse=".*"),
      paste(rev(tokens), collapse=".*")
    )
  }
  # allow simple stem patterns when helpful
  stems <- gsub("(s|ies)$", "", tokens, perl = TRUE)
  stems <- unique(stems[nchar(stems) > 2])
  if (length(stems) >= 2) pats <- unique(c(pats, paste(stems, collapse=".*")))
  pats
}

# score matched names to rank results
score_matches <- function(names_vec, keywords) {
  if (length(names_vec) == 0) return(numeric(0))
  ln <- tolower(names_vec)
  lk <- tolower(keywords)
  # string distance component (smaller is better)
  dmin <- vapply(ln, function(s) min(adist(s, lk)), numeric(1))
  # exact-substring bonus if any keyword fully contained
  bonus_contains <- vapply(ln, function(s) as.numeric(any(grepl(paste(lk, collapse="|"), s, fixed = TRUE))), numeric(1)) * 3
  # boundary-ish match bonus (keyword as whole word)
  boundary_re <- paste0("\\b(", paste(gsub("([.^$|()*+?{}\\[\\]\\\\])","\\\\\\1", lk), collapse="|"), ")\\b")
  bonus_boundary <- vapply(ln, function(s) as.numeric(grepl(boundary_re, s, perl=TRUE)), numeric(1)) * 2
  # starts/ends with bonus
  bonus_start <- vapply(ln, function(s) as.numeric(any(startsWith(s, lk))), numeric(1)) * 2
  bonus_end   <- vapply(ln, function(s) as.numeric(any(endsWith(s, lk))), numeric(1)) * 1
  # final score: higher is better
  score <- (-dmin) + bonus_contains + bonus_boundary + bonus_start + bonus_end
  names(score) <- names_vec
  score
}

# generate search patterns for a target name
patterns_for_target <- function(target_name) {
  mk <- manual_keywords[[target_name]]
  if (is.null(mk)) mk <- character(0)
  auto <- make_auto_patterns(target_name)
  unique(c(mk, auto))
}

# core search
# colnms: named list(file_path -> character vector of column names)
# returns invisibly a named list of data.frames with matches and counts
scan_target <- function(target_name, colnms, top_n = 20, show_examples_per_col = 3, verbose = TRUE) {
  pats <- patterns_for_target(target_name)
  if (length(pats) == 0) pats <- make_auto_patterns(target_name)
  if (length(pats) == 0) pats <- tolower(target_name)
  
  # build a single regex OR of all patterns; treat entries containing '.*' as regex, others as escaped substrings
  esc <- function(s) gsub("([.^$|()*+?{}\\[\\]\\\\])","\\\\\\1", s, perl = TRUE)
  re_parts <- ifelse(grepl("\\.\\*", pats), pats, esc(pats))
  big_re <- paste0("(", paste(re_parts, collapse = "|"), ")")
  
  # per-file matches
  matches_by_file <- lapply(colnms, function(v) {
    v[grepl(big_re, v, ignore.case = TRUE, perl = TRUE)]
  })
  
  # frequency across files (unique per file to avoid double-counting)
  per_file_unique <- unlist(lapply(matches_by_file, unique), use.names = FALSE)
  freq <- sort(table(per_file_unique), decreasing = TRUE)
  matched_cols <- names(freq)
  
  if (verbose) {
    cat("\n== ", target_name, " ==\n", sep = "")
    cat("# keywords: ", paste(pats, collapse = ", "), "\n", sep = "")
    if (length(matched_cols) == 0) {
      cat("no matches found.\n")
      return(invisible(list()))
    }
  }
  
  # rank by score + frequency
  scores <- score_matches(matched_cols, pats)
  # align frequency vector to score names
  freq_vec <- as.numeric(freq[matched_cols])
  names(freq_vec) <- matched_cols
  # weighted rank: prioritize frequency, then score
  rank_val <- freq_vec * 5 + scores
  ord <- order(rank_val, decreasing = TRUE, na.last = TRUE)
  ranked <- matched_cols[ord]
  
  # assemble small example of which files each column appears in
  examples <- lapply(ranked, function(col) {
    hits <- names(Filter(function(vec) any(tolower(vec) == tolower(col)), matches_by_file))
    base <- basename(hits)
    head(base, show_examples_per_col)
  })
  
  top <- head(ranked, top_n)
  for (nm in top) {
    cat("- ", nm, " [files:", freq_vec[[nm]], "]", sep = "")
    ex <- examples[[nm]]
    if (length(ex) > 0) cat(" e.g. {", paste(ex, collapse = ", "), "}", sep = "")
    cat("\n")
  }
  
  # best guess suggestion
  best_guess <- if (length(top)) top[1] else NA_character_
  cat("> suggested:", best_guess, "\n")
  
  # return structured info invisibly (could be used later to extract columns)
  res <- data.frame(
    match = ranked,
    files = freq_vec[ranked],
    score = as.numeric(scores[ranked]),
    stringsAsFactors = FALSE
  )
  invisible(list(
    target = target_name,
    patterns = pats,
    table = res,
    examples = examples
  ))
}

# iterate across all sd_ph columns (you can subset this vector if you prefer)
scan_all_targets <- function(sd_colnames, colnms, top_n = 20, verbose = TRUE) {
  out <- vector("list", length(sd_colnames))
  names(out) <- sd_colnames
  for (nm in sd_colnames) {
    if (verbose) cat("\n", strrep("=", 60), "\n", sep = "")
    out[[nm]] <- scan_target(nm, colnms, top_n = top_n, verbose = verbose)
  }
  invisible(out)
}

scan_all_targets(names(manual_keywords)[2:3], colnms, top_n = 20)


#### merge results ####

# candidate columns in priority order for each target variable
best_map_candidates <- list(
  pid                 = c("pid", "rec_id"),
  Timepoint           = c("visit_code", "d_visit", "days_visit"),
  study_label         = c("form_name"),
  sex                 = c("sex_psca"),
  age_years           = c("age_psca"),
  visit_code          = c("visit_code"),
  
  # anthropometrics
  body_weight         = c("wtkgavg_hwwt", "wtkg_hwwt", "wtkg3_hwwt", "wtkg2_hwwt", "wtkg1_hwwt",
                          "wtkg_dxas", "wtkg_cpet", "wtkg_acee", "wtkg_acre", "wtkg_pcaa", "wtkg_pcab"),
  body_height_cm      = c("htcmavg_hwwt", "htcm3_hwwt", "htcm2_hwwt", "htcm1_hwwt",
                          "hht_cpet", "sht_cpet", "htin_psca"),  # in_psca needs conversion
  bmi                 = c("calcbmi_hwwt", "bmi_psca"),
  waist_circum_cm     = c("wccmavg_hwwt", "wccm3_hwwt", "wccm2_hwwt", "wccm1_hwwt"),
  
  # resting BP/HR
  bp_sys              = c("sysavg_bphr", "sys3_bphr", "sys2_bphr", "sys1_bphr",
                          "bpsys_acee", "bpsys_acre", "sysbp_iske"),  # prefer BPHR
  bp_dia              = c("diasavg_bphr", "dias3_bphr", "dias2_bphr", "dias1_bphr",
                          "bpdia_acee", "bpdia_acre", "diabp_iske"),
  pulse_rest          = c("pulseavg_bphr", "pulse3_bphr", "pulse2_bphr", "pulse1_bphr"),
  
  # CPET
  vo2cart1            = c("vo2cart1_cpet"),
  vo2cart2            = c("vo2cart2_cpet"),
  wtkg_cpet           = c("wtkg_cpet"),
  hr_end_cpet         = c("hrend_cpet"),
  
  # peak outputs you said you already have
  watts_max           = c("wattpeak_cpet", "watts_max", "wrpeak_cpet", "workrate_peak", "wattend_cpet"),
  vo2_peak_L          = c("vo2peak_l", "vo2peak_cpet", "vo2_abs", "vo2_l"),
  vo2_peak_mlkg       = c("vo2_ml_kg_min", "vo2_rel", "vo2peak_mlkg"),
  
  # grip
  grip1               = c("trial1stat_grip", "trial1kg_grip"),
  grip2               = c("trial2stat_grip", "trial2kg_grip"),
  grip3               = c("trial3stat_grip", "trial3kg_grip"),
  handgrip_max        = c("handgrip_max"),  # usually derived; see below
  
  # strength 1RM (lbs exist; kg derived)
  leg_press_1rm_lbs   = c("lpmaxrest_fara"),
  leg_ext_1rm_lbs     = c("lemaxrest_fara"),
  chest_press_1rm_lbs = c("cpmaxrest_fara"),
  leg_press_1rm_kg    = c("leg_press_1rm_kg"),   # if native kg exists, else derive
  leg_ext_1rm_kg      = c("leg_ext_1rm_kg"),
  chest_press_1rm_kg  = c("chest_press_1rm_kg"),
  
  # DEXA
  lean_mass           = c("total_lean_mass", "wbtot_lean"),
  fat_mass            = c("total_fat_mass", "wbtot_fat"),  # avoid pfat which is percent
  total_bmd           = c("total_bmd", "wbtot_bmd"),
  vat_mass            = c("vat_mass", "vat_mass_unob"),
  
  # derived if VO2peak + HR available
  o2_pulse            = c("o2_pulse"),
  
  # torque (proxy if only isokinetic available)
  knee_iso_torque     = c("knee_iso_torque", "peaktorqavg_iske", "peaktorq_iske")
)

# pick first column in candidates that exists in df; returns column name or NA
first_present <- function(df, candidates) {
  candidates[ candidates %in% names(df) ][1]
}

# set or derive a column safely; modifies df in-place by reference when used inside a function
ensure_height_cm <- function(df) {
  if (!"body_height_cm" %in% names(df)) {
    nm <- first_present(df, best_map_candidates$body_height_cm)
    if (!is.na(nm)) {
      if (nm == "htin_psca") {
        df$body_height_cm <- df$htin_psca * 2.54
      } else {
        df$body_height_cm <- df[[nm]]
      }
    } else {
      df$body_height_cm <- NA_real_
    }
  }
  df
}

# grip max from trials if not already present
ensure_handgrip_max <- function(df) {
  if (!"handgrip_max" %in% names(df)) {
    t1 <- first_present(df, c("trial1stat_grip","trial1kg_grip"))
    t2 <- first_present(df, c("trial2stat_grip","trial2kg_grip"))
    t3 <- first_present(df, c("trial3stat_grip","trial3kg_grip"))
    if (!is.na(t1) || !is.na(t2) || !is.na(t3)) {
      a <- if (!is.na(t1)) df[[t1]] else NA
      b <- if (!is.na(t2)) df[[t2]] else NA
      c3 <- if (!is.na(t3)) df[[t3]] else NA
      df$handgrip_max <- pmax(a, b, c3, na.rm = TRUE)
    }
  }
  df
}

# 1RM kg derivations from lbs if kg missing
to_kg <- function(x) x * 0.453592
ensure_1rm_kg <- function(df) {
  if (!"leg_press_1rm_kg" %in% names(df))  if ("lpmaxrest_fara" %in% names(df))  df$leg_press_1rm_kg   <- to_kg(df$lpmaxrest_fara)
  if (!"leg_ext_1rm_kg"  %in% names(df))  if ("lemaxrest_fara" %in% names(df))  df$leg_ext_1rm_kg    <- to_kg(df$lemaxrest_fara)
  if (!"chest_press_1rm_kg"%in% names(df)) if ("cpmaxrest_fara" %in% names(df)) df$chest_press_1rm_kg<- to_kg(df$cpmaxrest_fara)
  df
}

# O2 pulse derived when possible (mL/beat)
ensure_o2_pulse <- function(df) {
  if (!"o2_pulse" %in% names(df)) {
    vnm <- first_present(df, best_map_candidates$vo2_peak_L)
    hnm <- first_present(df, c("hrend_cpet", "hr_peak_cpet"))
    if (!is.na(vnm) && !is.na(hnm)) {
      df$o2_pulse <- 1000 * df[[vnm]] / df[[hnm]]
    }
  }
  df
}

# choose a single column per target and copy to canonical name
apply_best_map <- function(df, map_candidates) {
  out <- df
  # resolve straightforward mappings
  for (target in names(map_candidates)) {
    if (target %in% c("body_height_cm","handgrip_max","leg_press_1rm_kg","leg_ext_1rm_kg","chest_press_1rm_kg","o2_pulse")) next
    src <- first_present(out, map_candidates[[target]])
    if (!is.na(src)) out[[target]] <- out[[src]]
  }
  # special derivations
  out <- ensure_height_cm(out)
  out <- ensure_handgrip_max(out)
  out <- ensure_1rm_kg(out)
  out <- ensure_o2_pulse(out)
  out
}

# ---- example usage on a merged row-level df ----
# df <- merge_all_your_sources_somehow(...)
# df_std <- apply_best_map(df, best_map_candidates)
