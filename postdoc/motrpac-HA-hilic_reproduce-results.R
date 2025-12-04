# ha_hilic_analysis.R -- demonstration script for processing MoTrPAC metabolomics data
#
# This script illustrates how to download, quality‑check and analyze hydrophilic interaction
# liquid chromatography (HILIC) metabolomics data from the MoTrPAC human pre‑COVID dataset.
# The example uses the T02 plasma HILIC positive ion (metab‑u‑hilicpos) dataset for the
# sedentary (SED) cohort, but the same steps can be adapted for the high‑activity (HA)
# cohort.  It assumes you have cloned the `precovid‑analyses` repository and configured
# a `motrpac_config.json` file in your home directory as described in the repository
# README.  See the README for details about the config file and directory layout.

## 1. Load configuration and required packages

# Read the config file.  The config needs at least the following elements:
#   - local_path: absolute path to the local clone of `precovid‑analyses` (no trailing slash)
#   - gsutil_path: path to the gsutil binary on your machine (e.g. "/usr/bin/gsutil")
#   - gsutil_user: optional user argument if needed by your installation
library(rjson)
config_path <- '~/motrpac_config.json'
if (!file.exists(config_path)) stop('Config file not found at ', config_path)
config <- rjson::fromJSON(file = config_path)

source(file.path(config$local_path, "library", "library-general.R"))
source(file.path(config$local_path, "QC", "generate_metab_qc_norm.R"))
source(file.path(config$motrpac_bic_norm_qc_repo_path, "tools", "unsupervised_normalization_functions.R"))
source(file.path(config$motrpac_bic_norm_qc_repo_path, "tools", "metabolomics_data_parsing_functions.R"))


# download SED metabolomics results expected by the QC script
dest_sed <- file.path(config$local_path, "data", "tmp", "metabolomics_qc_norm")
dir.create(dest_sed, recursive = TRUE, showWarnings = FALSE)
system2(config$gsutil_path, c(
  "-m", "cp", "-r",
  "gs://motrpac-data-hub/human-precovid/results/metabolomics-untargeted/t02-plasma/", dest_sed
))

library(MotrpacHumanPreSuspensionData) # make sure this is installed
source(file.path(config$local_path, "QC", "generate_metab_qc_norm.R"))
generate_metabolomics_qc_norm(
  repo_local_dir = config$local_path,
  bic_norm_qc_path = config$motrpac_bic_norm_qc_repo_path
)
