library(jsonlite)
library(data.table)
library(MotrpacHumanPreSuspensionData)

# 1. Load Config
setwd("~/repos/precovid-analyses/")
config <- fromJSON('~/motrpac_config.json')
local_repo <- config$local_path

# 1. Manually set the path to the BIC QC repo
# (Adjust this string if your folder name is slightly different)
bic_qc_path <- "~/repos/motrpac-bic-norm-qc"

# Verify it exists
if(!dir.exists(bic_qc_path)) stop(paste("Directory not found at:", bic_qc_path))

# 2. Source the required function libraries
source(file.path(bic_qc_path, "tools/metabolomics_data_parsing_functions.R"))
source(file.path(bic_qc_path, "tools/unsupervised_normalization_functions.R"))
source(file.path(bic_qc_path, "tools/qc_report_aux_functions.R"))
source(file.path(bic_qc_path, "tools/MetabolomicsLibrary.R"))

# --- 1. CLEAN SLATE ---
# Define your base working directory
base_work_dir <- file.path(local_repo, "tmp", "HA_hilic_repro")

# Delete it if it exists to ensure no phantom files
if(dir.exists(base_work_dir)) unlink(base_work_dir, recursive = TRUE)
dir.create(base_work_dir, recursive = TRUE)

# --- 2. DOWNLOAD TO STAGING ---
# We download to a "staging" folder first so we don't confuse the parser
staging_dir <- file.path(base_work_dir, "staging")
dir.create(staging_dir)

target_bucket <- "gs://motrpac-data-hub/human-precovid/results/metabolomics-untargeted/t02-plasma/metab-u-hilicpos/"
print("Downloading data...")
system(paste("gsutil -m cp -r", target_bucket, staging_dir))

# Check download success
downloaded <- list.files(staging_dir, recursive = T, full.names = T)
if(length(downloaded) == 0) stop("Download failed! No files in staging.")
print(paste("Downloaded", length(downloaded), "files."))

# --- 3. STRUCTURE THE DATA ---
# Create the specific hierarchy the parser demands
final_data_dir <- file.path(base_work_dir, "data")
struct_dir <- file.path(final_data_dir, "t02-plasma", "metab-u-hilicpos")
dir.create(struct_dir, recursive = TRUE)

# Move files from staging to the structured folder
# We filter for only the .txt files we need
txt_files <- downloaded[grepl("\\.txt$", downloaded)]
if(length(txt_files) == 0) stop("No .txt files found in downloaded data.")

file.copy(from = txt_files, 
          to = file.path(struct_dir, basename(txt_files)))

# --- 4. PREPARE THE OBJECT ---
# Normalize paths to handle Mac's /private/var quirks or symlinks
real_base_path <- normalizePath(final_data_dir, winslash = "/", mustWork = TRUE)
real_files <- normalizePath(list.files(struct_dir, full.names = TRUE), winslash = "/", mustWork = TRUE)

# Create the object
metab_obj = list(
  downloaded_files = real_files, 
  local_path = real_base_path
)

# --- 5. VERIFY MATH ---
# Test the "Folder Math" before running the function
rel_path <- gsub(metab_obj$local_path, "", metab_obj$downloaded_files[1], fixed = TRUE)
# Remove leading slash if present
rel_path <- sub("^/", "", rel_path)

print(paste("Base Path:", metab_obj$local_path))
print(paste("File Path:", metab_obj$downloaded_files[1]))
print(paste("Resulting Relative Path:", rel_path))

# Check if the math worked
if(!grepl("^t02-plasma", rel_path)) {
  stop("Folder Math failed! The relative path does not start with t02-plasma.")
}








# 1. Add the trailing slash to local_path
# This ensures that when we subtract the path, we subtract the slash too.
if (!grepl("/$", metab_obj$local_path)) {
  metab_obj$local_path <- paste0(metab_obj$local_path, "/")
}

# 2. Verify the Math again
# The result MUST NOT start with a slash
rel_path <- gsub(metab_obj$local_path, "", metab_obj$downloaded_files[1], fixed = TRUE)
print(paste("New Relative Path:", rel_path))

# 3. Run Parser
if(!grepl("^/", rel_path)){
  print("Leading slash removed. Parsing...")
  metab_datasets = read_metabolomics_datasets_from_download_obj(metab_obj, TRUE)
  metabolomics_parsed_datasets = metab_datasets$metabolomics_parsed_datasets
  
  print("--- SUCCESS ---")
  print(names(metabolomics_parsed_datasets))
} else {
  stop("Still has a leading slash! Check the paths manually.")
}
