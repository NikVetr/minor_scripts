############################################################
##  Dictionary explorer  –  robustness patch
############################################################
library(data.table)
library(MotrpacHumanPreSuspensionData)

#precompiled clinical table
MotrpacHumanPreSuspensionData::cln_curated_ex_performance_testing

## 1.  target patterns  (unchanged) -----------------------
patterns <- list(
  vo2_abs_L          = c("vo2 L","absolute vo2","vo2cart"),
  vo2_max_mlkg       = c("vo2 mL/kg","ml/kg/min","vo2max"),
  vt1_vo2            = c("VT1","ventilatory threshold 1"),
  vt2_vo2            = c("VT2","ventilatory threshold 2"),
  ve_vco2_slope      = c("VE/VCO2","ventilatory efficiency"),
  o2_pulse           = c("O2 pulse"),
  watts_max          = c("max watts","wattend","wattexer"),
  pulse_rest         = c("resting heart rate","hrend","pulse"),
  bp_sys             = c("systolic","sysend","sysbp"),
  bp_dia             = c("diastolic","diaend","diabp"),
  handgrip           = c("grip","hand grip"),
  knee_iso_torque    = c("peak torque","torque","nm","isometric knee"),
  leg_press_1rm      = c("leg press","1rm leg press"),
  leg_ext_1rm        = c("leg extension","1rm leg extension"),
  chest_press_1rm    = c("chest press","1rm chest"),
  lean_mass          = c("lean mass"),
  fat_mass           = c("fat mass"),
  dxa_bmd            = c("BMD","bone mineral density"),
  body_weight        = c("weight","wtkg"),
  body_height_cm     = c("height","htcm"),
  bmi                = c("bmi"),
  waist_circum_cm    = c("waist circumference"),
  steps_avg          = c("average steps","steps"),
  vecmag_avg         = c("vector magnitude")
)

#slightly more inclusive patterns
patterns$vo2_abs_L        <- c(patterns$vo2_abs_L,  "vo2cart", "peakvo2", "vo2_l")
patterns$vo2_max_mlkg     <- c(patterns$vo2_max_mlkg, "mlkg", "vo2max", "vo2_mlkg")
patterns$watts_max        <- c(patterns$watts_max, "peak watts")
patterns$handgrip         <- c(patterns$handgrip, "griplbs", "hgrip", "^gs")
patterns$steps_avg        <- c(patterns$steps_avg, "avgstep", "stepavg")
patterns$vecmag_avg       <- c(patterns$vecmag_avg, "vmag", "avgvm")
patterns$pulse_rest       <- c(patterns$pulse_rest, "hrrest", "resthr", "^hr0")
patterns$o2_pulse       <- c("o2 pulse",  "o2pulse", "o2_pulse")
patterns$vo2_max_mlkg   <- c("mlkg", "vo2max", "peakvo2", "vo2_mlkg")
patterns$vt1_vo2        <- c("vt1", "threshold 1")
patterns$vt2_vo2        <- c("vt2", "threshold 2")
patterns$ve_vco2_slope  <- c("ve/vco2", "vco2 slope", "vevco2")
patterns$vo2_max_mlkg <- c("mlkg", "vo2max", "peakvo2", "vo2_mlkg")
patterns$o2_pulse     <- c("o2pulse", "o2_pulse", "o2 pulse")
patterns$steps_avg    <- c("steps_per_day", "stepsavg", "avgsteps")
patterns$vecmag_avg   <- c("vector magnitude", "avgvm", "vmag")
patterns$vo2_max_mlkg <- c("mlkg", "vo2max", "vo2_mlkg", "peakvo2", "vo2maxmlkg")
patterns$steps_avg   <- c("steps_per_day", "avg_steps", "stepsavg")
patterns$vecmag_avg  <- c(patterns$vecmag_avg, "avg_vector_magnitude", "vmag")

## 2.  safe search function -------------------------------
search_dict <- function(dict, pat_vec) {
  if (is.null(dict) || !"field_name" %in% names(dict)) return(NULL)
  pat <- paste(pat_vec, collapse = "|")
  dict[
    grepl(pat, field_name_description, ignore.case = TRUE) |
      grepl(pat, field_name,             ignore.case = TRUE),
    .(field_name, field_name_description)
  ]
}

## 3.  iterate through all cln_raw tables -----------------
all_raw <- c(
  ls("package:MotrpacHumanPreSuspensionData", pattern = "^cln_raw_"),
  "cln_curated_accel_derived_variables_baseline"   # <- add this
)
results <- vector("list", length(patterns)); names(results) <- names(patterns)

for (var in names(patterns)) {
  cat(paste0("(", var, ") "))
  hits_all <- NULL
  for (tbl in all_raw) {
    data(list = tbl, envir = environment())
    dict <- get(tbl)$dict
    hit  <- search_dict(as.data.table(dict), patterns[[var]])
    if (!is.null(hit) && nrow(hit)) {
      hit[, table := tbl]
      hits_all <- rbind(hits_all, hit)
    }
  }
  if (!is.null(hits_all)) {
    setorder(hits_all, table, field_name)
    results[[var]] <- hits_all
  } else {
    results[[var]] <- NULL
  }
}

## 4.  print report --------------------------------------
names(results)
for (var in names(results)) {
  cat("\n\n", var, "\n\n", sep = "")
  if (is.null(results[[var]]) || nrow(results[[var]]) == 0L) {
    cat("  (no matches)\n")
  } else {
    print(results[[var]])
  }
}

#and recover the relevant variables
source_specs <- list(
  
  ## height / weight / waist, BP–HR
  cln_raw_height_weight_waist_circumference = c(
    pid             = "pid",
    visit_code      = "visit_code",
    body_weight     = "wtkgavg_hwwt",
    body_height_cm  = "htcmavg_hwwt",
    bmi             = "calcbmi_hwwt",
    waist_circum_cm = "wccmavg_hwwt"
  ),
  
  cln_raw_bld_pressure_heart_rate = c(
    pid        = "pid",
    bp_sys     = "sysavg_bphr",
    bp_dia     = "diasavg_bphr",
    pulse_rest = "pulseavg_bphr"
  ),
  
  ## CPET
  cln_raw_cardiopulmonary_ex_test = c(
    pid         = "pid",
    vo2cart1    = "vo2cart1_cpet",
    vo2cart2    = "vo2cart2_cpet",
    watts_max   = "wattend_cpet",
    wtkg_cpet   = "wtkg_cpet",
    hr_end_cpet = "hrend_cpet"          # keep if you’ll compute O2‑pulse
  ),
  
  ## strength
  cln_raw_grip_strength = c(
    pid    = "pid",
    grip1  = "trial1kg_grip",
    grip2  = "trial2kg_grip",
    grip3  = "trial3kg_grip"
  ),
  
  cln_raw_isometric_knee_extension = c(
    pid              = "pid",
    knee_iso_torque  = "peaktorq_iske"
  ),
  
  cln_raw_acute_resistance_ex_test = c(
    pid                = "pid",
    leg_press_1rm_lbs  = "lplbs1_acre",
    leg_ext_1rm_lbs    = "lelbs1_acre",
    chest_press_1rm_lbs= "cplbs1_acre"
  ),
  
  ## DXA (GE scanners)
  cln_raw_dxa_scan_results_ge = c(
    pid          = "pid",
    lean_mass    = "total_lean_mass",
    fat_mass     = "total_fat_mass",
    total_bmd    = "total_bmd",
    vat_mass     = "vat_mass"     # visceral adipose tissue
  ),
  
  ## accelerometer
  cln_curated_accel_derived_variables_baseline = c(
    pid         = "pid",
    steps_avg   = "daily_total_steps",
    vecmag_avg  = "daily_total_vm"
  )
)

pieces <- vector("list", length(source_specs))

i <- 1
for (tbl in names(source_specs)) {
  data(list = tbl, envir = environment())            # load object
  dt <- as.data.table(get(tbl)$data)                 # each motrdat has $data
  cols <- source_specs[[tbl]]
  
  sel <- dt[ pid %in% sed_pid &
               ( !"visit_code" %in% names(dt) | visit_code %in% baseline_visit ),
             ..cols ]                                # select & rename
  setnames(sel, old = cols, new = names(cols))
  pieces[[i]] <- sel
  i <- i + 1
}

clin_wide <- Reduce(function(x, y) merge(x, y, by = "pid", all = TRUE),
                    pieces)

## derive VO2, strength kg etc.
clin_wide[, vo2_peak_L     := pmax(vo2cart1, vo2cart2, na.rm = TRUE)]
clin_wide[, vo2_peak_mlkg  := vo2_peak_L / wtkg_cpet * 1000]
clin_wide[, o2_pulse       := vo2_peak_L / hr_end_cpet]
clin_wide[, `:=`(
  leg_press_1rm_kg    = leg_press_1rm_lbs  * 0.453592,
  leg_ext_1rm_kg      = leg_ext_1rm_lbs    * 0.453592,
  chest_press_1rm_kg  = chest_press_1rm_lbs* 0.453592,
  handgrip_max        = pmax(grip1, grip2, grip3, na.rm = TRUE)
)]

clin_wide <- as.data.frame(clin_wide)
clin_wide <- clin_wide[,(colMeans(apply(clin_wide, 2, is.na)) < 0.99)]

fwrite(clin_wide, "~/data/smontgom/motrpac_precovid_sedentary_clinical_obs.txt")
