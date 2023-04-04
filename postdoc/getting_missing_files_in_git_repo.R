files <- c("/Volumes/SSD500GB/gtex-pipeline/expression_genotypes/GTEx_v8_noChr_noRSID.bim",
"/Volumes/SSD500GB/gtex-pipeline/GTEx_Analysis_v8_eQTL_covariates/",
"/Volumes/SSD500GB/gtex-pipeline/GTEx_Analysis_v8_eQTL_expression_matrices/",
"/Volumes/SSD500GB/gtex-pipeline/log2-normalized-expression/log2-normalized-expression_*.expression.bed.gz",
"~/data/smontgom/41467_2021_23579_MOESM6_ESM.csv",
"~/data/smontgom/est_gcor_mat.RData",
"~/data/smontgom/GENES_HUMAN.txt",
"~/data/smontgom/GENES_RAT.txt",
"~/data/smontgom/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt",
"~/data/smontgom/GTEx_v8_ExpressionScores/tissues/",
"~/data/smontgom/gwas_metadata.csv",
"~/data/smontgom/imputed_gwas_hg38_1.1/",
"~/data/smontgom/meta_analysis_results.RData",
"~/data/smontgom/old_dea_deseq_20201121/*_training-dea_20201121.RData",
"~/data/smontgom/open-targets_tissue-x-disease_*",
"~/data/smontgom/opentargets/associationByOverallDirect.csv",
"~/data/smontgom/opentargets/associationByOverallDirect.csv",
"~/data/smontgom/opentargets/associationByOverallIndirect.csv",
"~/data/smontgom/PANTHER17_human_rat_ref_genome_orthologs.tsv",
"~/data/smontgom/relative_expression_motrpac_gtex",
"~/data/smontgom/RGD_ORTHOLOGS_20201001.txt",
"~/data/smontgom/RSID_POS_MAP_*.txt",
"~/data/smontgom/zcor_transcriptome_pass1b.tsv",
"~/repos/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.*.bim",
"~/repos/ldsc/baseline/baseline.*.annot",
"~/repos/ldsc/baseline/baseline.*.annot",
"~/repos/ldsc/custom_genesets/cluster_*.chr_*.annot",
"~/repos/ldsc/ENSG_coord.txt")

wildcards <- grepl("\\*", files)
single_files <- files[!wildcards]
single_sizes <- setNames(file.info(single_files)$size, single_files)

multi_files <- files[wildcards]
headless <- function(x, n) head(x, length(x) - n)
folders <- paste0(sapply(lapply(strsplit(multi_files, "/"), headless, 1), paste0, collapse = "/"), "/")
wild_filenames <- sapply(strsplit(multi_files, "/"), tail, 1)

multi_sizes <- lapply(setNames(seq_along(folders), folders), function(i){
  files_in_folder <- list.files(folders[i])
  matches <- grepl(glob2rx(wild_filenames[i]), files_in_folder)
  matching_filenames <- files_in_folder[matches]
  matching_paths <- paste0(folders[i], matching_filenames)
  setNames(file.info(matching_paths)$size, matching_paths)
})

#folder contents
folder_sizes <- sapply(multi_sizes, sum)
hist(log10(folder_sizes) / 3, breaks = 100:500/100); abline(v = 8/3, lwd = 2, col = 2)

#indiv files in folders
hist(log10(unlist(multi_sizes)) / 3, breaks = 100:300/100); abline(v = 8/3, lwd = 2, col = 2)

#indiv files
hist(log10(unlist(single_sizes)) / 3, breaks = 100:300/100); abline(v = 8/3, lwd = 2, col = 2)

#total size
sum(c(single_sizes, unlist(multi_sizes))) / 1E9
sum(single_sizes) / 1E9 #github can accommodate, but goal is <1GB total and <50MB per file
sort(single_sizes / 1E6)
sum(sort(folder_sizes / 1E9, T)[-1])
sort(single_sizes / 1E9, T)
head(sort(c(single_sizes, unlist(multi_sizes)) / 1E9, T))

#game plan -- we have 63 GB of data total, put all single files < 0.1 GB (or GiB? 1E9 bytes) onto github, everything else onto Zenodo incl stuff in folders

#single files 

github_files_local <- names(single_sizes)[single_sizes / 1E9 < 0.1]
github_data_dir <- "~/repos/MoTrPAC_Complex_Traits/data/"

github_file_map <- cbind(from = single_sizes, to = paste0(github_data_dir, gsub("*../", "", github_files_local)))