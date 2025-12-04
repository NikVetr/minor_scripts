library(arrow)
library(data.table)

# copy everything locally (or narrow the wildcard if you want)
data_dir <- "~/data/ha_biospecimen/"
dir.create(data_dir, showWarnings = FALSE)
system(paste0("gsutil -m cp -r gs://motrpac-data-hub/human/phenotype/curated/biospecimen/* ", data_dir))

# read all parquet/csv into one data.table
files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
read_one <- function(p) {
  if (grepl("\\.parquet$", p, ignore.case = TRUE)) as.data.table(read_parquet(p)) else fread(p)
}
ha_raw <- rbindlist(lapply(files, read_one), use.names = TRUE, fill = TRUE)

# start from ha_raw
setDT(ha_raw)

# simple "first non-NA" helper that won't error if all NA
first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 0) x[1] else NA
}

# One row per pid
ha_demo <- ha_raw[!is.na(pid), .(
  activity_level = "Adult Highly Active",
  sex            = first_non_na(sex),
  age_years      = as.numeric(first_non_na(calculatedAge)),
  bmi            = as.numeric(first_non_na(bmi)),
  modality       = first_non_na(randomGroupCode),  # e.g. ADUEndur / ADUResist
  visitcode_chr  = first_non_na(visitcode)
), by = pid]

fwrite(ha_demo, "~/data/motrpac_ha_clinical_obs.txt")

