library(data.table)

#download data if it is not already downloaded
# data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_sets/"
data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_dictionary/"
# data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/raw/data_sets/"
data_dir <- "~/data/human-main/phenotype/human-main-ha-adu/raw/data_dictionary/"
dir.create(data_dir, showWarnings = FALSE, recursive = T)
# path_to_file <- "gs://motrpac-data-hub/human-main/phenotype/human-main-ha-adu/curated/data_sets/*"
# path_to_file <- "gs://motrpac-data-hub/human-main/phenotype/human-main-ha-adu/curated/data_dictionary/*"
path_to_file <- "gs://motrpac-data-hub/human-main/phenotype/human-main-ha-adu/raw/data_dictionary/*"
# path_to_file <- "gs://motrpac-data-hub/human-main/phenotype/human-main-ha-adu/raw/data_sets/*"
# system(paste0("gsutil -m cp -r ", path_to_file, " ", data_dir))
  
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

#language model output

# ----------- paths -----------
cur_data_dir  <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_sets/"
cur_dict_dir  <- "~/data/human-main/phenotype/human-main-ha-adu/curated/data_dictionary/"
raw_data_dir  <- "~/data/human-main/phenotype/human-main-ha-adu/raw/data_sets/"
raw_dict_dir  <- "~/data/human-main/phenotype/human-main-ha-adu/raw/data_dictionary/"
out_dir       <- "~/data/human-main/phenotype/human-main-ha-adu/derived/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------- expected schema from your sedentary data -----------
# assumes sd_ph is already in your workspace
expected <- colnames(sd_ph)

# ----------- helpers: normalization + small utilities -----------
normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("%", "pct", x)
  x <- gsub("kg|_kg|lbs|_lbs|cm|_cm|mmhg", "", x)  # strip common unit tokens
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}
dataset_short <- function(s) gsub("\\..*$","", s)

# ----------- synonyms for names (left = normalized expected) -----------
synonyms <- list(
  pid                  = c("pid","rec_id"),
  timepoint            = c("timepoint","timepointorder","timepointversion","phase","tp"),
  study_label          = c("study","study_label"),
  sex                  = c("sex","gender"),
  age_years            = c("calculatedage","age","age_years"),
  visit_code           = c("visit_code","visitcode","visit"),
  body_weight          = c("wtkg","weight","body_weight"),
  body_height_cm       = c("htcm","height","body_height_cm"),
  bmi                  = c("bmi","body_mass_index"),
  waist_circum_cm      = c("wccm","waist_circum","waist","waist_circum_cm"),
  bp_sys               = c("sys_bp","sbp","systolic","avg_ecgsys_cpet"),
  bp_dia               = c("dias_bp","dbp","diastolic","avg_ecgdia_cpet"),
  pulse_rest           = c("hr_min","rest_hr","bpm_rest"),
  vo2cart1             = c("avg_vo2cart_cpet","vo2_avg","vo2_cart1"),
  vo2cart2             = c("avg_vco2cart_cpet","vco2_avg","vco2_cart"),
  watts_max            = c("wattend_cpet","watts_peak","max_watts"),
  wtkg_cpet            = c("wtkg_cpet"),
  hr_end_cpet          = c("hrend_cpet","hr_end_cpet"),
  grip1                = c("grip1"),
  grip2                = c("grip2"),
  grip3                = c("grip3"),
  knee_iso_torque      = c("peaktorq_iske","knee_isometric_torque","peaktorqavg_iske"),
  leg_press_1rm_lbs    = c("lp_1rm_lbs","leg_press_1rm_lbs"),
  leg_ext_1rm_lbs      = c("le_1rm_lbs","leg_extension_1rm_lbs"),
  chest_press_1rm_lbs  = c("cp_1rm_lbs","chest_press_1rm_lbs"),
  lean_mass            = c("lean_mass","total_lean_mass","lean_body_mass"),
  fat_mass             = c("fat_mass","total_fat_mass"),
  total_bmd            = c("total_bmd","whole_body_bmd","bmd_total"),
  vat_mass             = c("vat_mass","visceral_adipose_tissue_mass","vat_vol","vat_volume"),
  vo2_peak_L           = c("vo2max_l","vo2_peak_l"),
  vo2_peak_mlkg        = c("vo2max","vo2_ml_kg_min","vo2_peak_mlkg"),
  o2_pulse             = c("o2_pulse","o2_pulse_avg"),
  leg_press_1rm_kg     = c("lp_1rm_kg"),
  leg_ext_1rm_kg       = c("le_1rm_kg"),
  chest_press_1rm_kg   = c("cp_1rm_kg"),
  handgrip_max         = c("peak_handgrip","handgrip_max"),
  specialization_group = c("specialization_group","training_specialization","endurance_resistance_group")
)

# description keywords (helps find BMD/DEXA/etc.)
desc_keywords <- list(
  total_bmd       = c("bmd","bone mineral density","dexa","whole body"),
  lean_mass       = c("lean mass","dexa","fat free","whole body"),
  fat_mass        = c("fat mass","adipose","dexa","whole body"),
  vat_mass        = c("visceral","vat","adipose","abdominal"),
  handgrip_max    = c("handgrip","grip","dynamometer","max"),
  knee_iso_torque = c("knee","isometric","torque"),
  leg_press_1rm_kg   = c("leg press","1rm","one repetition","maximum","kg"),
  leg_ext_1rm_kg     = c("leg extension","1rm","one repetition","maximum","kg"),
  chest_press_1rm_kg = c("chest press","1rm","one repetition","maximum","kg")
)

# ----------- curated subset → harmonize → collapse -----------

curated_files <- c(
  "human-main-ha-adu_clinical_acute_bout_ds_crf-curated_v1.0.txt",
  "human-main-ha-adu_clinical_biospecimen_ds_crf-curated_v1.1.txt",
  "human-main-ha-adu_clinical_screening_ds_crf-curated_v1.0.txt"
)

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
  # exact
  for (i in seq_along(norm_expected)) {
    idx <- match(norm_expected[i], norm_observed)
    if (!is.na(idx)) {
      map$match_type[i]  <- "exact"
      map$observed[i]    <- observed_names[idx]
      map$observed_norm[i] <- norm_observed[idx]
      map$edit_dist[i]   <- 0L
    }
  }
  # synonym
  for (i in which(is.na(map$match_type))) {
    ekey <- map$expected_norm[i]
    syns <- synonyms[[ekey]]
    if (!is.null(syns) && length(syns) > 0) {
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
  # fuzzy
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
  map$match_type[is.na(map$match_type)] <- "missing"
  map
}

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

harmonized <- setNames(vector("list", length(curated_files)), curated_files)
for (k in seq_along(curated_files)) {
  p <- file.path(cur_data_dir, curated_files[k])
  message("[DEBUG] reading curated: ", p)
  df <- as.data.frame(fread(p))
  m  <- build_column_map(expected, colnames(df), max_edit_dist = 2)
  
  message("[DEBUG] map summary (", curated_files[k], "): ",
          " exact=",   sum(m$match_type=="exact"),
          " synonym=", sum(m$match_type=="synonym"),
          " fuzzy=",   sum(m$match_type=="fuzzy"),
          " missing=", sum(m$match_type=="missing"))
  
  hz <- apply_column_map(df, m)
  spec <- if (grepl("screening", curated_files[k])) "endurance"
  else if (grepl("acute_bout", curated_files[k])) "resistance"
  else NA
  hz$specialization_group <- spec
  harmonized[[k]] <- hz[, c(expected, "specialization_group")]
}

ha_combined <- do.call(rbind, harmonized)
rownames(ha_combined) <- NULL

# ----------- collapse by key -----------
use_visit_code_in_key <- TRUE
num_tol <- 1e-6

is_numish <- function(x) {
  nx <- suppressWarnings(as.numeric(as.character(x)))
  ok <- sum(!is.na(nx)); total <- sum(!is.na(x))
  if (total == 0) return(FALSE)
  ok / total > 0.8
}
to_num_or_char <- function(x) if (is_numish(x)) suppressWarnings(as.numeric(as.character(x))) else as.character(x)
all_equalish <- function(x, tol=num_tol) {
  x <- x[!is.na(x)]; if (length(x) <= 1) return(TRUE)
  if (is_numish(x)) { xn <- suppressWarnings(as.numeric(as.character(x))); xn <- xn[!is.na(xn)]
  if (length(xn) <= 1) return(TRUE); return(max(xn) - min(xn) <= tol)
  } else { length(unique(as.character(x))) == 1 }
}
consensus_or_na <- function(x, tol=num_tol) if (all_equalish(x,tol)) { y <- x[!is.na(x)]; if (!length(y)) NA else if (is_numish(y)) suppressWarnings(as.numeric(as.character(y[1]))) else as.character(y[1]) } else NA

collapse_group <- function(df_group, tol=num_tol) {
  cols <- colnames(df_group); out <- vector("list", length(cols)); names(out) <- cols
  conflicts <- character(0)
  for (cl in cols) {
    val <- consensus_or_na(df_group[[cl]], tol)
    if (is.na(val) && sum(!is.na(df_group[[cl]])) > 1) conflicts <- c(conflicts, cl)
    out[[cl]] <- val
  }
  out <- as.data.frame(out, stringsAsFactors=FALSE)
  out$n_rows_merged <- nrow(df_group)
  out$conflict_cols  <- if (length(conflicts)) paste(conflicts, collapse=";") else ""
  out
}

collapse_ha_by_key <- function(ha, use_visit_code=FALSE, tol=num_tol) {
  stopifnot(all(c("pid","Timepoint") %in% colnames(ha)))
  key_cols <- c("pid","Timepoint")
  if (use_visit_code && "visit_code" %in% colnames(ha)) key_cols <- c(key_cols,"visit_code")
  for (k in key_cols) ha[[k]] <- as.character(ha[[k]])
  key_str <- do.call(paste, c(ha[key_cols], sep="||"))
  idx_list <- split(seq_len(nrow(ha)), key_str)
  message("[DEBUG] collapsing groups: ", length(idx_list))
  res <- vector("list", length(idx_list)); i <- 0L
  for (nm in names(idx_list)) {
    i <- i + 1L; if (i %% 1000 == 0) message("[DEBUG] processed groups: ", i)
    g <- ha[idx_list[[nm]], , drop=FALSE]
    for (cl in colnames(g)) g[[cl]] <- to_num_or_char(g[[cl]])
    res[[i]] <- collapse_group(g, tol)
  }
  collapsed <- do.call(rbind, res)
  collapsed <- collapsed[order(collapsed$pid, collapsed$Timepoint), ]
  has_conf <- nchar(collapsed$conflict_cols) > 0
  message("[DEBUG] collapsed rows: ", nrow(collapsed), " | groups with conflicts: ",
          sum(has_conf), " (", round(100*mean(has_conf),1), "%)")
  list(data=collapsed,
       per_column_conflicts=sort(table(unlist(strsplit(paste(collapsed$conflict_cols[has_conf], collapse=";"), ";", fixed=TRUE))), decreasing=TRUE))
}

res_collapse  <- collapse_ha_by_key(ha_combined, use_visit_code_in_key, num_tol)
ha_collapsed  <- res_collapse$data
message("[DEBUG] original rows: ", nrow(ha_combined), " → collapsed: ", nrow(ha_collapsed))

# ----------- build unified catalogs: columns + BOTH dictionaries -----------
sniff_read_delim <- function(fp) {
  con <- file(fp, open="r"); on.exit(close(con))
  first <- readLines(con, n=1)
  sep <- if (grepl("\t", first)) "\t" else ","
  read.csv(fp, sep=sep, stringsAsFactors=FALSE, quote="", check.names=FALSE)
}

read_dict_dir <- function(dirpath) {
  files <- list.files(dirpath, full.names = TRUE)
  canonical <- c(
    "field_name","data_set_name","version","external_data_release","data_type",
    "field_name_description","calculated_created_variable","categorical_values",
    "category_definitions"
  )
  
  alias_map <- list(
    field_name                 = c("field_name","field","variable","variable_name","name"),
    data_set_name              = c("data_set_name","dataset_name","dataset","data_set","form_name","form","instrument"),
    version                    = c("version","ver"),
    external_data_release      = c("external_data_release","release_level","external_release","release"),
    data_type                  = c("data_type","datatype","type"),
    field_name_description     = c("field_name_description","field_description","variable_description",
                                   "description","item_label","question_text","label"),
    calculated_created_variable= c("calculated_created_variable","calc_created_var","calculated_variable",
                                   "derived","computed","created_variable"),
    categorical_values         = c("categorical_values","categories","values","levels","codes"),
    category_definitions       = c("category_definitions","category_defs","value_labels","labels","definitions")
  )
  
  normalize_header <- function(nm) {
    nm <- tolower(trimws(nm))
    nm <- gsub("[^a-z0-9]+","_", nm)
    gsub("^_|_$","", nm)
  }
  
  apply_aliases <- function(df) {
    names(df) <- normalize_header(names(df))
    # rename aliases to canonical
    for (canon in names(alias_map)) {
      if (!(canon %in% names(df))) {
        hit <- alias_map[[canon]][alias_map[[canon]] %in% names(df)]
        if (length(hit)) {
          names(df)[match(hit[1], names(df))] <- canon
        }
      }
    }
    # add any missing canonical columns as NA
    for (canon in canonical) if (!(canon %in% names(df))) df[[canon]] <- NA
    df <- df[, canonical, drop = FALSE]
    df
  }
  
  out_list <- vector("list", length(files))
  names(out_list) <- basename(files)
  
  for (i in seq_along(files)) {
    f <- files[i]
    # fread handles csv/tsv and oddities; fill=TRUE handles ragged rows
    df <- tryCatch(
      as.data.frame(fread(f, sep = "auto", quote = "", fill = TRUE, showProgress = FALSE)),
      error = function(e) {
        warning("[WARN] could not read dict: ", basename(f), " | ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(df)) next
    df <- apply_aliases(df)
    df$dict_file <- basename(f)
    out_list[[i]] <- df
  }
  
  # bind whatever we managed to load
  do.call(rbind, out_list[!vapply(out_list, is.null, logical(1))])
}

dict_cur <- read_dict_dir(cur_dict_dir)
dict_raw <- read_dict_dir(raw_dict_dir)
dict_cur$source <- "curated"
dict_raw$source <- "raw"
dict_all <- rbind(dict_cur, dict_raw)
if (!nrow(dict_all)) stop("No dictionary rows found; check dict directories.")

dict_all$field_name_norm  <- normalize_name(dict_all$field_name)
dict_all$data_set_norm    <- normalize_name(dict_all$data_set_name)
dict_all$desc_norm        <- tolower(dict_all$field_name_description)

# column catalogs from BOTH data_dirs (don’t read full files, just headers)
# robust header reader: fread first; fallback to first non-empty line
get_colnames_safe <- function(fp) {
  # try fast path
  cn <- tryCatch(colnames(fread(fp, nrows = 0L, showProgress = FALSE)),
                 error = function(e) NULL)
  if (!is.null(cn) && length(cn) > 0) return(cn)
  
  # fallback: sniff first non-empty, non-comment line
  con <- file(fp, open = "r"); on.exit(close(con))
  line <- NULL
  for (i in 1:50) {
    l <- suppressWarnings(readLines(con, n = 1L, warn = FALSE))
    if (length(l) == 0) break
    s <- trimws(l)
    if (nchar(s) == 0) next
    if (grepl("^(#|//)", s)) next
    line <- s
    break
  }
  if (is.null(line)) return(character(0))
  
  # strip UTF-8 BOM if present
  line <- sub("^\ufeff", "", line)
  
  delims <- c("\t", ",", "|", ";")
  counts <- sapply(delims, function(d) length(strsplit(line, d, fixed = TRUE)[[1]]))
  delim <- delims[which.max(counts)]
  
  parts <- strsplit(line, delim, fixed = TRUE)[[1]]
  parts <- gsub('^"|"$', '', parts)
  parts <- trimws(parts)
  parts
}

# build long-form catalog; gracefully skip unreadable files
catalog_from_dir <- function(dirpath, source_label) {
  files <- list.files(dirpath, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)
  # keep obvious delimited text files
  files <- files[grepl("\\.(txt|csv|tsv)$", basename(files), ignore.case = TRUE)]
  if (!length(files)) {
    return(data.frame(
      source = character(0),
      data_set_file = character(0),
      data_set_base = character(0),
      data_set_name = character(0),
      colname = character(0),
      col_norm = character(0),
      data_set_norm = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  rows <- vector("list", length(files))
  n <- 0L
  for (f in files) {
    cols <- tryCatch(get_colnames_safe(f), error = function(e) character(0))
    if (length(cols) == 0) {
      message("[DEBUG] skipping (no header detected): ", basename(f))
      next
    }
    n <- n + 1L
    base <- basename(f)
    dsn  <- gsub("\\..*$", "", base)
    df <- data.frame(
      source        = source_label,
      data_set_file = f,
      data_set_base = base,
      data_set_name = dsn,
      colname       = cols,
      stringsAsFactors = FALSE
    )
    df$col_norm      <- normalize_name(df$colname)
    df$data_set_norm <- normalize_name(df$data_set_name)
    rows[[n]] <- df
  }
  
  if (n == 0L) {
    return(data.frame(
      source = character(0),
      data_set_file = character(0),
      data_set_base = character(0),
      data_set_name = character(0),
      colname = character(0),
      col_norm = character(0),
      data_set_norm = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  do.call(rbind, rows[seq_len(n)])
}

cat_cur <- catalog_from_dir(cur_data_dir, "curated")
cat_raw <- catalog_from_dir(raw_data_dir, "raw")
cat_long <- rbind(cat_cur, cat_raw)

# attach dictionary description where field+dataset align; else fallback by field only
match1 <- match(paste(cat_long$col_norm, cat_long$data_set_norm),
                paste(dict_all$field_name_norm, dict_all$data_set_norm))
desc1 <- ifelse(is.na(match1), "", dict_all$desc_norm[match1])
dt1   <- ifelse(is.na(match1), "", dict_all$data_type[match1])

ifelse_empty <- function(a, b) ifelse(nzchar(a), a, b)

match2 <- match(cat_long$col_norm, dict_all$field_name_norm)
desc2  <- ifelse(is.na(match2), "", dict_all$desc_norm[match2])
dt2    <- ifelse(is.na(match2), "", dict_all$data_type[match2])

cat_long$desc <- ifelse_empty(desc1, desc2)
cat_long$data_type <- ifelse_empty(dt1, dt2)

# ----------- rank candidates for filling missing vars -----------
lev_dist <- function(a,b) { if (length(a)==0 || length(b)==0) return(Inf); as.integer(adist(a,b)) }

score_candidate <- function(expected_key_norm, col_norm, desc, syns, keywords, source_label) {
  nm_score   <- as.numeric(identical(expected_key_norm, col_norm))
  syn_score  <- as.numeric(col_norm %in% syns)
  d          <- lev_dist(expected_key_norm, col_norm)
  fuzzy_score<- if (is.finite(d)) max(0, (2 - min(d,2)))/2 else 0  # 1.0 for d=0, 0.5 for d=1, 0 for >=2
  kw_score   <- 0
  if (!is.na(desc) && nchar(desc) > 0 && length(keywords) > 0) {
    dn <- tolower(desc)
    hits <- sapply(keywords, function(k) grepl(k, dn, fixed=FALSE))
    kw_score <- sum(hits)/max(1, length(keywords))
  }
  # tiny bias toward curated if tie
  curated_bonus <- if (identical(source_label, "curated")) 0.05 else 0
  2*nm_score + 1.5*syn_score + 1.0*fuzzy_score + 1.0*kw_score + curated_bonus
}

rank_candidates_for <- function(expected_var) {
  e_norm <- normalize_name(expected_var)
  syns <- synonyms[[e_norm]]; if (is.null(syns)) syns <- character(0)
  kws  <- desc_keywords[[e_norm]]; if (is.null(kws)) kws <- character(0)
  scores <- mapply(function(cn, ds, src) score_candidate(e_norm, cn, ds, syns, kws, src),
                   cn=cat_long$col_norm, ds=cat_long$desc, src=cat_long$source)
  out <- data.frame(
    expected=expected_var,
    source=cat_long$source,
    data_set_base=cat_long$data_set_base,
    data_set_name=cat_long$data_set_name,
    colname=cat_long$colname,
    col_norm=cat_long$col_norm,
    score=as.numeric(scores),
    desc=ifelse(is.na(cat_long$desc),"",cat_long$desc),
    stringsAsFactors=FALSE
  )
  out <- out[out$score > 0, ]
  out[order(-out$score), ]
}

# choose targets = columns with >= 50% NA in collapsed table
na_rates <- round(colMeans(is.na(ha_collapsed[, intersect(colnames(ha_collapsed), expected), drop=FALSE])), 2)
targets  <- names(na_rates)[na_rates >= 0.50]
message("[DEBUG] targeting (>=50% NA): ", paste(targets, collapse=", "))

shortlists <- lapply(targets, rank_candidates_for); names(shortlists) <- targets
all_short  <- do.call(rbind, shortlists)
write.csv(all_short, file=file.path(out_dir, "candidate_columns_ranked_curated_and_raw.csv"), row.names=FALSE)
for (nm in names(shortlists)) {
  write.csv(head(shortlists[[nm]], 20),
            file=file.path(out_dir, paste0("candidates_", normalize_name(nm), "_top20.csv")),
            row.names=FALSE)
}
message("[DEBUG] wrote candidate CSVs in: ", out_dir)

# ----------- auto-merge best picks (non-destructive fill of NAs) -----------
best_picks <- do.call(rbind, lapply(names(shortlists), function(nm) head(shortlists[[nm]], 1)))
best_picks <- best_picks[!is.na(best_picks$score) & best_picks$score > 0, ]

# index files → dir so we can read from either curated or raw
file_dir_lookup <- rbind(
  data.frame(source="curated", data_set_base=list.files(cur_data_dir), dir=cur_data_dir, stringsAsFactors=FALSE),
  data.frame(source="raw",     data_set_base=list.files(raw_data_dir), dir=raw_data_dir,  stringsAsFactors=FALSE)
)

read_selected_columns <- function(source_label, basefile, needed_cols) {
  ddir <- file_dir_lookup$dir[file_dir_lookup$source==source_label & file_dir_lookup$data_set_base==basefile]
  if (!length(ddir)) return(NULL)
  fp <- file.path(ddir[1], basefile)
  have <- tryCatch(colnames(data.table::fread(fp, nrows=0L)), error=function(e) character(0))
  sel <- intersect(needed_cols, have)
  if (!length(sel)) return(NULL)
  as.data.frame(data.table::fread(fp, select=sel, showProgress=FALSE))
}

# merging util: fill NAs only; prefer (pid, visit_code); fallback to pid
merge_in_df <- function(dst, src, cols_map, by_keys) {
  # cols_map: named list/vector (dst_name = src_name)
  map_src <- as.character(unlist(cols_map, use.names = FALSE))
  map_dst <- names(cols_map)
  
  # drop NA/empty
  keep_m <- !(is.na(map_src) | map_src == "" | is.na(map_dst) | map_dst == "")
  map_src <- map_src[keep_m]; map_dst <- map_dst[keep_m]
  
  # never rename keys (either as source or destination)
  key_cols <- by_keys
  not_key_idx <- !(map_src %in% key_cols | map_dst %in% key_cols)
  map_src <- map_src[not_key_idx]
  map_dst <- map_dst[not_key_idx]
  
  if (length(map_src) == 0) {
    # nothing to bring in besides keys; still need to ensure keys present
    map_src <- character(0); map_dst <- character(0)
  }
  
  # Subset source to needed columns
  keep <- unique(c(by_keys, map_src))
  has  <- intersect(keep, names(src))
  if (length(has) == 0) {
    message("[DEBUG] merge_in_df: no requested columns in source; skipping")
    return(dst)
  }
  src2 <- src[, has, drop = FALSE]
  
  # if pid missing but rec_id present, make pid
  if (!("pid" %in% names(src2)) && ("rec_id" %in% names(src))) {
    src2$pid <- src$rec_id
  }
  
  # ensure keys exist; if visit_code missing, fallback to pid only
  if (!all(by_keys %in% names(src2))) {
    if ("pid" %in% names(src2)) {
      by_keys <- "pid"
      message("[DEBUG] merge_in_df: falling back to by='pid' (", paste(setdiff(c("pid","visit_code"), names(src2)), collapse=", "), " missing)")
    } else {
      message("[DEBUG] merge_in_df: no merge keys present; skipping")
      return(dst)
    }
  }
  
  # coerce keys to character in both frames
  for (k in by_keys) {
    if (k %in% names(src2)) src2[[k]] <- as.character(src2[[k]])
    if (k %in% names(dst))  dst[[k]]  <- as.character(dst[[k]])
  }
  
  # rename non-key mappings (only if present)
  if (length(map_src)) {
    for (i in seq_along(map_src)) {
      sc <- map_src[i]; dc <- map_dst[i]
      if (sc %in% names(src2) && !(dc %in% names(src2))) {
        names(src2)[names(src2) == sc] <- dc
      }
    }
  }
  
  # de-dup colnames defensively
  names(src2) <- make.unique(names(src2), sep="__dup")
  
  # Revalidate keys after rename (should still be there because we never rename keys)
  if (!all(by_keys %in% names(src2))) {
    message("[DEBUG] merge_in_df: keys disappeared after rename; skipping")
    return(dst)
  }
  
  # merge
  m <- merge(dst, src2, by = by_keys, all.x = TRUE, suffixes = c("", ".src"))
  
  # fill only where dst is NA
  if (length(map_dst)) {
    for (dc in map_dst) {
      sc <- paste0(dc, ".src")
      if (sc %in% names(m)) {
        idx <- is.na(m[[dc]]) & !is.na(m[[sc]])
        if (any(idx)) m[[dc]][idx] <- m[[sc]][idx]
        m[[sc]] <- NULL
      }
    }
  }
  m
}

ha_enhanced <- ha_collapsed
sources_split <- split(best_picks, paste(best_picks$source, best_picks$data_set_base, sep="|"))

for (key in names(sources_split)) {
  picks <- sources_split[[key]]
  
  # drop NA/empty colnames
  picks <- subset(picks, !is.na(colname) & nzchar(colname))
  
  # hard guard: never use keys as source for a non-key expected
  key_cols <- c("pid","visit_code")
  picks <- subset(picks, !(colname %in% key_cols & !(expected %in% key_cols)))
  
  parts <- strsplit(key, "\\|")[[1]]
  src_label <- parts[1]; basefile <- parts[2]
  needed <- unique(c(picks$colname, "pid","rec_id","visit_code"))
  
  src_df <- read_selected_columns(src_label, basefile, needed)
  if (is.null(src_df)) next
  
  by_keys <- c("pid","visit_code")
  if (!all(by_keys %in% names(src_df)) && !("pid" %in% names(src_df) || "rec_id" %in% names(src_df))) next
  if (!all(by_keys %in% names(src_df))) by_keys <- "pid"
  
  cols_map <- structure(as.list(picks$colname), names = picks$expected)
  message("[DEBUG] merging from [", src_label, "] ", basefile, " by: ", paste(by_keys, collapse=", "))
  ha_enhanced <- merge_in_df(ha_enhanced, src_df, cols_map, by_keys)
}


# ----------- quick missingness deltas -----------
miss_before <- colMeans(is.na(ha_collapsed[, intersect(colnames(ha_collapsed), expected), drop=FALSE]))
miss_after  <- colMeans(is.na(ha_enhanced[, intersect(colnames(ha_enhanced), expected), drop=FALSE]))
delta <- round(miss_before - miss_after, 3)
improve <- sort(delta[delta != 0], decreasing=TRUE)
message("[DEBUG] missingness reduction (top 15):")
print(utils::head(improve, 15))

# Optionally write outputs
# write.csv(ha_enhanced, file=file.path(out_dir, "ha_collapsed_enhanced.csv"), row.names=FALSE)

# ----------- inspection helpers -----------
inspect_group <- function(pid_val, timepoint_val, visit_code_val=NULL, src=ha_combined) {
  rows <- src$pid == pid_val & src$Timepoint == timepoint_val
  if (!is.null(visit_code_val) && "visit_code" %in% names(src)) rows <- rows & src$visit_code == visit_code_val
  message("[DEBUG] rows in group: ", sum(rows)); src[rows, , drop=FALSE]
}
collapsed_conflicts_for <- function(pid_val, timepoint_val, visit_code_val=NULL, src=ha_enhanced) {
  rows <- src$pid == pid_val & src$Timepoint == timepoint_val
  if (!is.null(visit_code_val) && "visit_code" %in% names(src)) rows <- rows & src$visit_code == visit_code_val
  if (!any(rows)) return(character(0))
  cc <- src$conflict_cols[which(rows)[1]]; if (!nzchar(cc)) character(0) else unique(strsplit(cc, ";", fixed=TRUE)[[1]])
}

colMeans(is.na(ha_enhanced))

# ----------- minimal unit tests -----------
# dictionaries have expected columns
if (nrow(dict_all)) stopifnot(all(c("field_name","data_set_name") %in% colnames(dict_all)))
# ranking returns something for a classic DEXA field like total_bmd (likely present)
rbmd <- rank_candidates_for("total_bmd"); stopifnot(nrow(rbmd) >= 1)
# merge fills NA but not overwrite
x <- data.frame(pid=c("1","2"), visit_code=c("V1","V1"), total_bmd=c(NA, 0.99), stringsAsFactors=FALSE)
y <- data.frame(pid=c("1","2"), visit_code=c("V1","V1"), bmd_total=c(1.05, 0.88), stringsAsFactors=FALSE)
cm <- c(total_bmd="bmd_total")
z <- merge_in_df(x, y, cm, c("pid","visit_code")); stopifnot(isTRUE(all.equal(z$total_bmd, c(1.05,0.99))))
message("[DEBUG] unit tests passed.")
