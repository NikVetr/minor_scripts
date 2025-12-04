
#load sed data first
sd_ph <- as.data.frame(fread("~/data/motrpac_precovid_sedentary_clinical_obs.txt"))

#print out column names for language model to parse
data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_sets/"
ha_paths <- list.files(data_dir)
for(i in seq_along(ha_paths)){
  ha_path <- paste0(data_dir, ha_paths[i])
  ha_ph <- as.data.frame(fread(ha_path))
  cat(paste0("\n\n", ha_path, ":\n"))
  cat(paste0("\"", colnames(ha_ph), "\","))
}

cat(paste0("\"", colnames(sd_ph), "\","))

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

#can try combining like strings?
all_strings <- unique(unlist(colnms))
str(all_strings)
all_strings[grepl("hrexer", all_strings)]

#language model output
expected <- colnames(sd_ph)

# simple normalizer: lower, remove non-alnum, collapse sequences
normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("%", "pct", x)
  x <- gsub("kg|_kg|lbs|_lbs|cm|_cm|mmhg", "", x)   # strip units tokens commonly embedded
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# seed synonym map: left = normalized expected, right = normalized observed(s)
# multiple observed synonyms allowed per expected (vector)
synonyms <- list(
  pid                  = c("pid"),
  timepoint            = c("timepoint","timepointorder","timepointversion"), # may not truly match semantically
  study_label          = character(0),
  sex                  = c("sex"),
  age_years            = c("calculatedage","age","age_years"),
  visit_code           = c("visit_code","visitcode"),
  body_weight          = c("wtkg","body_weight"),
  body_height_cm       = c("htcm","body_height","body_height_cm"),
  bmi                  = c("bmi"),
  waist_circum_cm      = c("wccm","waist","waist_circum"),
  bp_sys               = c("sys_bp","avg_ecgsys_cpet","ecgsys"), # last two are imperfect fallbacks
  bp_dia               = c("dias_bp","avg_ecgdia_cpet","ecgdia"),
  pulse_rest           = c("hr_min","hr_rest","rest_hr"),
  vo2cart1             = c("avg_vo2cart_cpet","vo2_avg"),
  vo2cart2             = c("avg_vco2cart_cpet","vco2_avg"),
  watts_max            = c("wattend_cpet","watts_max","watts_peak"),
  wtkg_cpet            = c("wtkg_cpet"),
  hr_end_cpet          = c("hrend_cpet","hr_end_cpet"),
  grip1                = c("grip1","peak_handgrip"),   # note: not same; peak used as proxy if needed
  grip2                = c("grip2","avg_handgrip"),
  grip3                = c("grip3","sd_handgrip"),
  knee_iso_torque      = c("peaktorq_iske","peaktorqavg_iske"),
  leg_press_1rm_lbs    = c("lp_1rm_lbs"),
  leg_ext_1rm_lbs      = c("le_1rm_lbs"),
  chest_press_1rm_lbs  = c("cp_1rm_lbs"),
  lean_mass            = c("lean_mass"),
  fat_mass             = c("fat_mass"),
  total_bmd            = c("total_bmd"),
  vat_mass             = c("vat_mass"),
  vo2_peak_L           = c("vo2max_l","vo2max_l_","vo2max_l__"),
  vo2_peak_mlkg        = c("vo2max","vo2_ml_kg_min"),
  o2_pulse             = c("o2_pulse","o2_pulse_avg"),
  leg_press_1rm_kg     = c("lp_1rm_kg","lp_avg_load_kg"), # latter ≠ 1RM; proxy only if you accept it
  leg_ext_1rm_kg       = c("le_1rm_kg","le_avg_load_kg"),
  chest_press_1rm_kg   = c("cp_1rm_kg","cp_avg_load_kg"),
  handgrip_max         = c("peak_handgrip")
)

# ---------- mapping helpers ----------
build_column_map <- function(expected_names, observed_names, max_edit_dist = 2) {
  norm_expected <- normalize_name(expected_names)
  norm_observed <- normalize_name(observed_names)
  
  map <- data.frame(
    expected      = expected_names,
    expected_norm = norm_expected,
    match_type    = NA_character_,
    observed      = NA_character_,
    observed_norm = NA_character_,
    edit_dist     = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # exact matches
  for (i in seq_along(norm_expected)) {
    idx <- match(norm_expected[i], norm_observed)
    if (!is.na(idx)) {
      map$match_type[i]  <- "exact"
      map$observed[i]    <- observed_names[idx]
      map$observed_norm[i] <- norm_observed[idx]
      map$edit_dist[i]   <- 0L
    }
  }
  
  # synonym matches where still unmatched
  for (i in which(is.na(map$match_type))) {
    ekey <- map$expected_norm[i]
    syns <- synonyms[[ekey]]
    if (!is.null(syns) && length(syns) > 0) {
      # pick first present synonym
      hits <- match(syns, norm_observed)
      if (any(!is.na(hits))) {
        idx <- hits[which(!is.na(hits))[1]]
        map$match_type[i]  <- "synonym"
        map$observed[i]    <- observed_names[idx]
        map$observed_norm[i] <- norm_observed[idx]
        map$edit_dist[i]   <- 0L
      }
    }
  }
  
  # fuzzy matches (Levenshtein) for remaining
  for (i in which(is.na(map$match_type))) {
    d <- adist(map$expected_norm[i], norm_observed)
    j <- which.min(d)
    if (length(j) == 1 && is.finite(d[1, j]) && d[1, j] <= max_edit_dist) {
      map$match_type[i]  <- "fuzzy"
      map$observed[i]    <- observed_names[j]
      map$observed_norm[i] <- norm_observed[j]
      map$edit_dist[i]   <- as.integer(d[1, j])
    }
  }
  
  # mark missing
  map$match_type[is.na(map$match_type)] <- "missing"
  
  map
}

# apply map to a data.frame to produce harmonized columns in expected order
apply_column_map <- function(df, col_map) {
  out <- as.list(rep(NA, nrow(col_map)))
  names(out) <- col_map$expected
  for (i in seq_len(nrow(col_map))) {
    if (col_map$match_type[i] %in% c("exact","synonym","fuzzy")) {
      out[[i]] <- df[[ col_map$observed[i] ]]
    } else {
      out[[i]] <- NA
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# your three files (paths as in your output)
files <- c(
  "human-main-ha-adu_clinical_acute_bout_ds_crf-curated_v1.0.txt",
  "human-main-ha-adu_clinical_biospecimen_ds_crf-curated_v1.1.txt",
  "human-main-ha-adu_clinical_screening_ds_crf-curated_v1.0.txt"
)

# root dir you showed
root <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_sets/"

library(data.table)  # for fread only; rest is base

maps <- setNames(vector("list", length(files)), files)
harmonized <- setNames(vector("list", length(files)), files)

for (k in seq_along(files)) {
  p <- file.path(root, files[k])
  message("[DEBUG] reading: ", p)
  df <- as.data.frame(fread(p))
  
  # build map
  m <- build_column_map(expected, colnames(df), max_edit_dist = 2)
  maps[[k]] <- m
  
  # print a compact summary
  message("[DEBUG] summary for: ", files[k])
  message("  exact:   ", sum(m$match_type == "exact"))
  message("  synonym: ", sum(m$match_type == "synonym"))
  message("  fuzzy:   ", sum(m$match_type == "fuzzy"))
  message("  missing: ", sum(m$match_type == "missing"))
  
  # construct harmonized dataset in expected order (NAs where missing)
  hz <- apply_column_map(df, m)
  
  # add specialization_group (best-effort default; override as needed)
  # heuristic: screening file -> endurance; acute_bout (with cp/le/op/lp load/volume cols) -> resistance
  spec <- if (grepl("screening", files[k])) "endurance" else if (grepl("acute_bout", files[k])) "resistance" else NA
  hz$specialization_group <- spec
  
  harmonized[[k]] <- hz
}

# combine only rows that actually belong to HA cohort (all three are HA); rbind with fill (base approach)
# ensure same column order first
for (k in seq_along(harmonized)) {
  harmonized[[k]] <- harmonized[[k]][, c(expected, "specialization_group")]
}
ha_combined <- do.call(rbind, harmonized)
rownames(ha_combined) <- NULL

# -------------------- config --------------------
# set this TRUE if you want to group on (pid, Timepoint, visit_code)
use_visit_code_in_key <- T

# numeric equality tolerance (e.g., 1e-6 means 0.1000001 == 0.1)
num_tol <- 1e-6

# -------------------- helpers --------------------
is_numish <- function(x) {
  # treat columns as numeric if most non-NA entries can be parsed as numbers
  nx <- suppressWarnings(as.numeric(as.character(x)))
  ok <- sum(!is.na(nx))
  total <- sum(!is.na(x))
  if (total == 0) return(FALSE)
  ok / total > 0.8
}

to_num_or_char <- function(x) {
  if (is_numish(x)) {
    return(suppressWarnings(as.numeric(as.character(x))))
  } else {
    return(as.character(x))
  }
}

all_equalish <- function(x, tol = num_tol) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(TRUE)
  # if numeric-like, compare numerically
  if (is_numish(x)) {
    xn <- suppressWarnings(as.numeric(as.character(x)))
    xn <- xn[!is.na(xn)]
    if (length(xn) <= 1) return(TRUE)
    return(max(xn) - min(xn) <= tol)
  } else {
    # character exact match
    ux <- unique(as.character(x))
    return(length(ux) == 1)
  }
}

consensus_or_na <- function(x, tol = num_tol) {
  if (all_equalish(x, tol = tol)) {
    # return the (only) non-NA value (or NA if all NA)
    y <- x[!is.na(x)]
    if (length(y) == 0) return(NA)
    # prefer a numeric if numish
    if (is_numish(y)) {
      return(suppressWarnings(as.numeric(as.character(y[1]))))
    } else {
      return(as.character(y[1]))
    }
  } else {
    return(NA)
  }
}

# collapse one data.frame group into a single row + conflict listing
collapse_group <- function(df_group, tol = num_tol) {
  cols <- colnames(df_group)
  out <- vector("list", length(cols))
  names(out) <- cols
  
  conflicts <- character(0)
  
  for (cl in cols) {
    val <- consensus_or_na(df_group[[cl]], tol = tol)
    if (is.na(val) && sum(!is.na(df_group[[cl]])) > 1) {
      # more than one non-NA in group but could not agree -> conflict
      conflicts <- c(conflicts, cl)
    }
    out[[cl]] <- val
  }
  
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out$n_rows_merged <- nrow(df_group)
  out$conflict_cols <- if (length(conflicts) == 0) "" else paste(conflicts, collapse = ";")
  out
}

# -------------------- main collapsing routine --------------------
collapse_ha_by_key <- function(ha, use_visit_code = FALSE, tol = num_tol) {
  stopifnot(all(c("pid","Timepoint") %in% colnames(ha)))
  key_cols <- c("pid", "Timepoint")
  if (use_visit_code && "visit_code" %in% colnames(ha)) {
    key_cols <- c(key_cols, "visit_code")
  }
  
  # ensure key columns exist as character
  for (k in key_cols) {
    ha[[k]] <- as.character(ha[[k]])
  }
  
  # split indices by key
  key_str <- do.call(paste, c(ha[key_cols], sep = "||"))
  idx_list <- split(seq_len(nrow(ha)), key_str)
  
  # iterate groups
  message("[DEBUG] total groups: ", length(idx_list))
  res <- vector("list", length(idx_list))
  ki <- 0L
  for (nm in names(idx_list)) {
    ki <- ki + 1L
    if (ki %% 1000 == 0) message("[DEBUG] processed groups: ", ki)
    gidx <- idx_list[[nm]]
    df_g <- ha[gidx, , drop = FALSE]
    
    # make comparisons fair: coerce each column to numeric where appropriate
    for (cl in colnames(df_g)) {
      df_g[[cl]] <- to_num_or_char(df_g[[cl]])
    }
    
    res[[ki]] <- collapse_group(df_g, tol = tol)
  }
  
  collapsed <- do.call(rbind, res)
  # ensure key columns are preserved from their collapsed values
  # (they should be identical/consensus by construction)
  collapsed <- collapsed[order(collapsed$pid, collapsed$Timepoint), ]
  
  # diagnostics
  has_conf <- nchar(collapsed$conflict_cols) > 0
  message("[DEBUG] collapsed rows: ", nrow(collapsed))
  message("[DEBUG] groups with conflicts: ", sum(has_conf), " (", round(100*mean(has_conf),1), "%)")
  
  # per-column conflict rate
  all_conf_cols <- unlist(strsplit(paste(collapsed$conflict_cols[has_conf], collapse = ";"), ";", fixed = TRUE))
  all_conf_cols <- all_conf_cols[nzchar(all_conf_cols)]
  col_conf_table <- sort(table(all_conf_cols), decreasing = TRUE)
  
  list(
    data = collapsed,
    per_column_conflicts = col_conf_table
  )
}

# -------------------- run on your ha_combined --------------------
# assumes ha_combined already in your env
message("[DEBUG] starting collapse on ha_combined (group = pid, Timepoint)")
res_basic <- collapse_ha_by_key(ha_combined, use_visit_code = use_visit_code_in_key, tol = num_tol)

ha_collapsed <- res_basic$data
conf_counts <- res_basic$per_column_conflicts

# quick sanity summaries
message("[DEBUG] original rows: ", nrow(ha_combined))
message("[DEBUG] collapsed rows: ", nrow(ha_collapsed))
message("[DEBUG] top 15 conflicting columns (if any):")
print(utils::head(conf_counts, 15))

colMeans(is.na(ha_collapsed))
