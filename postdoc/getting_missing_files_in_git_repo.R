#### functions ####
buffer_list <- function(x){
  max_entries <- max(sapply(x, length))
  lapply(x, function(i) c(i, rep("", max_entries - length(i))))
}

headless <- function(x, n) head(x, length(x) - n)

getNestedDirs <- function(path, make = T) {
  
  dirs <- path
  while (dirname(dirs[length(dirs)]) != dirs[length(dirs)]) {
    dirs <- c(dirs, dirname(dirs[length(dirs)]))
  }
  dirs <- rev(dirs)
  
  if(make){
    for(i in 1:length(dirs)){
      if(!dir.exists(dirs[i])){
        dir.create(dirs[i])
      }
    }
  } else {
    return(gsub("//", "/", paste0(dirs, "/")))
  }
  
}

copy_directory <- function(from, to) { #copies contents, but not parent directory
  
  if(!dir.exists(from)){
    warning(paste0("Directory <<", from, ">> does not exist."))
    return(F)
  }
  
  if (!dir.exists(to)) {
    getNestedDirs(to)
  }
  
  files <- list.files(from, full.names = TRUE, recursive = FALSE)
  
  for (file in files) {
    
    file_name <- basename(file)
    dest_path <- file.path(to, file_name)
    
    if (file.info(file)$isdir) {
      copy_directory(file, dest_path)
    } else {
      file.copy(file, dest_path)
    }
  }
  
  return(T
         )
}

dir.copy <- function(from, to){
  for(i in 1:length(from)){
    if(basename(to[i]) != basename(from[i])){
      dest_dir <- file.path(to[i], basename(from[i]))
      copy_directory(from[i], dest_dir)
    } else {
      copy_directory(from[i], to[i])  
    } 
  }
}

move_file <- function(from, to) {
  
  print(paste0("Now moving <<", from, ">>"))
  
  #check that source file exists
  if(!file.exists(from)){
    warning(paste0("File <<", from, ">> does not exist."))
    return(F)
  }
  
  #check that target directory exists
  if (!dir.exists(dirname(to))) {
    getNestedDirs(dirname(to))
  }
  
  # Copy the file to the new location
  if(!file.exists(to) || (file.info(to)$size != file.info(from)$size)){
    result <- file.copy(from, to, overwrite = T)
  } else {
    result <- T
  }
  success <- result & file.exists(to) & (file.info(from)$size == file.info(to)$size)
  
  # If the copy was successful, delete the original file
  if(file.access(from, 2) != -1){
    if (success) {
      file.remove(from)
    } else {
      warning("File was not successfully copied.")
    }  
  } else {
    warning(paste0("File <<", from, ">> ", " does not have write permissions."))
    return(F)
  }
  
  
  return(success)
  
}

file.move <- function(from, to){
  for(i in seq_along(from)){
    if(basename(to[i]) != basename(from[i])){
      dest_path <- file.path(to[i], basename(from[i]))
      move_file(from[i], dest_path)
    } else {
      move_file(from[i], to[i])
    } 
  }
}


equal.directories <- function(dir1, dir2, hash = F){
  
  #check both dirs exist
  if(!dir.exists(dir1)){
    warning(paste0("Dir <<", dir1, ">> does not exist."))
    return(F)
  }
  
  if(!dir.exists(dir2)){
    warning(paste0("Dir <<", dir2, ">> does not exist."))
    return(F)
  }
  
  #fix names of directories (number of slashes)
  dir1 <- gsub("/+", "/", paste0(dir1, "/"))
  dir2 <- gsub("/+", "/", paste0(dir2, "/"))
  
  # List all files and directories in both directories
  files1 <- list.files(dir1, full.names = TRUE, recursive = TRUE)
  files2 <- list.files(dir2, full.names = TRUE, recursive = TRUE)
  files1 <- gsub("/+", "/", files1)
  files2 <- gsub("/+", "/", files2)
  
  # Create relative paths for comparison
  relative1 <- gsub(paste0("^", normalizePath(dir1)), "", files1)
  relative2 <- gsub(paste0("^", normalizePath(dir2)), "", files2)
  
  # Check if both directories have the same files and subdirectories
  if (!identical(sort(relative1), sort(relative2))) {
    return(FALSE)
  }
  
  # Check if all files are of the same size
  for (i in seq_along(files1)) {
    if (file.info(files1[i])$isdir) {
      next
    }
    
    if(hash){ #more reliable, but slower (less change of collision)
      
      if (!identical(digest::digest(files1[i], file = TRUE), digest::digest(files2[match(relative1[i], relative2)], file = TRUE))) {
        return(FALSE)
      }
      
    } else {
      
      if (!identical(file.info(files1[i])$size, file.info(files2[match(relative1[i], relative2)])$size)) {
        return(FALSE)
      }
      
    }
    
  }
  
  return(TRUE)
}

move_directory <- function(from, to, per_file = T) {
  
  cat(paste0("\n\nNow moving <<", from, ">>\n\n"))
  
  if(!dir.exists(from)){
    warning(paste0("Directory <<", from, ">> does not exist."))
    return(F)
  }
  
  if (!dir.exists(to)) {
    getNestedDirs(to)
  }
  
  files <- list.files(from, full.names = TRUE, recursive = FALSE)
  
  for (file in files) {
    
    file_name <- basename(file)
    dest_path <- file.path(to, file_name)
    
    if (file.info(file)$isdir) {
      
      copy_directory(file, dest_path)
      
    } else {
      
      if(per_file){
        #saves on storage, deleting each individual file as you go
        #but if there's an error in the middle, your directories may be bonked
        file.move(file, dest_path)
      } else {
        file.copy(file, dest_path)  
      }
      
    }
  }
  
  if(!per_file){
    success <- dir.exists(to) & equal.directories(from, to)
    if(file.access(from, 2) != -1){
      if (success) {
        unlink(normalizePath(from), recursive = TRUE)
      } else {
        if(!dir.exists(to)){"Destination dir does not exist."}
        if(!equal.directories(from, to)){"Source and dest dirs not equal."}
        warning("Directory was not successfully copied.")    }  
    } else {
      warning(paste0("Directory <<", from, ">> ", " does not have write permissions."))
    }
  
    return(success)
  }
  
}

dir.move <- function(from, to, per_file = F){
  for(i in 1:length(from)){
    if(basename(to[i]) != basename(from[i])){
      dest_dir <- file.path(to[i], basename(from[i]))
      move_directory(from[i], dest_dir, per_file)
    } else {
      move_directory(from[i], to[i], per_file)  
    } 
  }
}


#### start processing ####
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

#read output from ~/scripts/montgomery_lab/editing_companion_filepaths.R
files <- unlist(read.csv(file = "~/data/smontgom/pass1b_companion_paper_external_filenames.csv")[,2])

wildcards <- grepl("\\*", files)
ends_in_slash <- substr(files, nchar(files), nchar(files)) == "/"

#process single files w/ sizes
single_files <- files[!wildcards & !ends_in_slash]
single_sizes <- setNames(file.info(single_files)$size, single_files)

#get wildcard file info
multi_files <- files[wildcards]
folders <- paste0(sapply(lapply(strsplit(multi_files, "/"), headless, 1), 
                         paste0, collapse = "/"), "/")
wild_filenames <- sapply(strsplit(multi_files, "/"), tail, 1)
folder_names <- folders
folder_names_excl <- c("smontgom", "SSD500GB", "repos") #to create new folder names in target dir
for(i in 1:length(folder_names_excl)){
  folder_names <- gsub(paste0(".*", folder_names_excl[i]), "", folder_names)
}

multi_sizes <- lapply(setNames(seq_along(folders), multi_files), function(i){
  files_in_folder <- list.files(folders[i])
  matches <- grepl(glob2rx(wild_filenames[i]), files_in_folder)
  if(sum(matches) == 0){print(paste0("error in #", i, "; ", folders[i], wild_filenames[i]))}
  matching_filenames <- files_in_folder[matches]
  matching_paths <- paste0(folders[i], matching_filenames)
  setNames(file.info(matching_paths)$size, matching_paths)
})
wild_folder_sizes <- sapply(multi_sizes, sum) 

#get whole folder info
multi_files_folders <- files[ends_in_slash]
multi_file_folder_names <- multi_files_folders
for(i in 1:length(folder_names_excl)){
  multi_file_folder_names <- gsub(paste0(".*", folder_names_excl[i]), "", multi_file_folder_names)
}
folder_sizes <- setNames(file.info(multi_files_folders)$size, multi_files_folders)


#folder contents
hist(log10(wild_folder_sizes) / 3); abline(v = 8/3, lwd = 2, col = 2)

#indiv files in folders
hist(log10(unlist(multi_sizes)) / 3); abline(v = 8/3, lwd = 2, col = 2)

#indiv files
hist(log10(unlist(single_sizes)) / 3); abline(v = 8/3, lwd = 2, col = 2)

#total size
sum(c(single_sizes, unlist(multi_sizes))) / 1E9
sum(single_sizes) / 1E9 #github can accommodate, but goal is <1GB total and <50MB per file
sort(single_sizes / 1E6)
sum(sort(wild_folder_sizes / 1E9, T)[-1])
sort(single_sizes / 1E9, T)
head(sort(c(single_sizes, unlist(multi_sizes)) / 1E9, T))

#check that everything is accounted for
length(files)
length(single_sizes) + length(multi_sizes) + length(wild_folder_sizes)
names(single_sizes) %in% files
names(multi_sizes) %in% files
names(folder_sizes) %in% files
files[!(files %in% 
    c(names(single_sizes), names(multi_sizes), names(folder_sizes)))]

#create reference material for later script editing
files_map <- data.frame(old_path = files, wild = wildcards, folders = ends_in_slash, single = !ends_in_slash & !wildcards)
files_map$new_dir <- ""

#game plan -- we have 63 GB of data total, no wait ~200GB, put all single files < 0.1 GB (or GiB? 1E9 bytes) onto github, 
#everything else onto Zenodo incl all stuff in folders

#single files internal
github_data_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/"
github_files_local <- names(single_sizes)[single_sizes / 1E9 < 0.1]
files_map$new_dir[files_map$old_path %in% github_files_local] <- github_data_dir
github_files_local_names <- github_files_local
for(i in 1:length(folder_names_excl)){
  github_files_local_names <- gsub(paste0(".*", paste0(folder_names_excl[i], "/")), "", github_files_local_names)
}
github_file_map <- data.frame(from = github_files_local, to = paste0(github_data_dir, github_files_local_names))

#parsing the hierarchy of the folders (to iteratively construct everything needed in target dir)
folder_names_single <- github_files_local_names
folder_hier_single <- do.call(rbind, buffer_list(strsplit(github_files_local_names, "/"))) #hier stands for hierarchy
ndepth_single <- apply(folder_hier_single, 1, function(x) sum(x != ""))
folder_hier_single <- matrix(folder_hier_single[ndepth_single > 1,], ncol = ncol(folder_hier_single))
ndepth_single <- apply(folder_hier_single, 1, function(x) sum(x != ""))
folder_hier_single <- matrix(folder_hier_single[!apply(folder_hier_single, 1, 
                            function(x) all(x == "")),], ncol = ncol(folder_hier_single))
folder_hier_single <- matrix(folder_hier_single[!duplicated(folder_hier_single),], ncol = ncol(folder_hier_single))

folders_to_make_single <- lapply(1:ncol(folder_hier_single), function(i){
  if(i==1){
    unique(paste0(github_data_dir, folder_hier_single[,1:i]))
  } else {
    unique(paste0(github_data_dir, apply(matrix(folder_hier_single[,1:i], ncol = ncol(folder_hier_single)), 1, paste0, collapse = "/"))[ndepth_single > i])
  }
})

#make the necessary folders
for(i in 1:length(folders_to_make_single)){
  if(length(folders_to_make_single[[i]]) > 0){
    for(j in 1:length(folders_to_make_single[[i]])){
      if(!dir.exists(folders_to_make_single[[i]][j])){
      dir.create(folders_to_make_single[[i]][j])  
      }
    }  
  }
}

#single files external
github_data_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"
github_files_external <- names(single_sizes)[single_sizes / 1E9 >= 0.1]
files_map$new_dir[files_map$old_path %in% github_files_external] <- github_data_dir
github_files_external_names <- github_files_external
for(i in 1:length(folder_names_excl)){
  github_files_external_names <- gsub(paste0(".*", paste0(folder_names_excl[i], "/")), "", github_files_external_names)
}
github_file_map_external <- data.frame(from = github_files_external, to = paste0(github_data_dir, github_files_external_names))

folder_names_single_external <- github_files_external_names
folder_hier_single_external <- do.call(rbind, buffer_list(strsplit(github_files_external_names, "/")))
ndepth_single_external <- apply(folder_hier_single_external, 1, function(x) sum(x != ""))
folder_hier_single_external <- folder_hier_single_external[ndepth_single_external > 1,]
ndepth_single_external <- apply(folder_hier_single_external, 1, function(x) sum(x != ""))
folder_hier_single_external <- folder_hier_single_external[!apply(folder_hier_single_external, 1, function(x) all(x == "")),]
folder_hier_single_external <- folder_hier_single_external[!duplicated(folder_hier_single_external),]

folders_to_make_single_external <- lapply(1:ncol(folder_hier_single_external), function(i){
  if(i==1){
    unique(paste0(github_data_dir, folder_hier_single_external[,1:i]))
  } else {
    unique(paste0(github_data_dir, apply(folder_hier_single_external[,1:i], 1, paste0, collapse = "/"))[ndepth_single_external > i])
  }
})

for(i in 1:length(folders_to_make_single_external)){
  if(length(folders_to_make_single_external[[i]]) > 0){
    for(j in 1:length(folders_to_make_single_external[[i]])){
      if(!dir.exists(folders_to_make_single_external[[i]][j])){
        dir.create(folders_to_make_single_external[[i]][j])  
      }
    }  
  }
}


#multifiles (wildcards)
github_data_dir <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/"
folder_names #this has the wildcard parent folders
wild_filenames #this has the wildcard file names corresponding to the above

folder_hier <- do.call(rbind, buffer_list(strsplit(folder_names, "/")))[,-1]
folder_hier <- folder_hier[!apply(folder_hier, 1, function(x) all(x == "")),]
folder_hier <- folder_hier[!duplicated(folder_hier),]
ndepth <- apply(folder_hier, 1, function(x) sum(x != ""))
folders_to_make <- lapply(1:ncol(folder_hier), function(i){
  if(i==1){
    unique(paste0(github_data_dir, folder_hier[,1:i]))
  } else {
    unique(paste0(github_data_dir, apply(folder_hier[,1:i], 1, paste0, collapse = "/"))[ndepth >= i])
  }
})

for(i in 1:length(folders_to_make)){
  for(j in 1:length(folders_to_make[[i]])){
    if(!dir.exists(folders_to_make[[i]][j])){
      dir.create(folders_to_make[[i]][j])
    }
  }
}

github_file_map_wildfiles <- lapply(1:length(folders), function(i){
  folder_files <- list.files(folders[i])
  matches <- folder_files[grep(pattern = glob2rx(wild_filenames[i]), x = folder_files)]
  data.frame(from = paste0(folders[i], matches), 
             to = paste0(substr(github_data_dir, 1, nchar(github_data_dir)-1), folder_names[i], matches)
  )
})

github_file_map_wildfiles <- do.call(rbind, github_file_map_wildfiles)

#multifiles (folders)
multi_files_folders #these are the full local paths of the below folders
multi_file_folder_names #these are the target folders

folder_hier_for_folders <- do.call(rbind, buffer_list(strsplit(multi_file_folder_names, "/")))[,-1]
folder_hier_for_folders <- folder_hier_for_folders[!apply(folder_hier_for_folders, 1, function(x) all(x == "")),]
folder_hier_for_folders <- folder_hier_for_folders[!duplicated(folder_hier_for_folders),]
ndepth_for_folders <- apply(folder_hier_for_folders, 1, function(x) sum(x != ""))
folders_to_make_for_folders <- lapply(1:ncol(folder_hier_for_folders), function(i){
  if(i==1){
    (paste0(github_data_dir, folder_hier_for_folders[,1:i])[ndepth_for_folders > i])
  } else {
    (paste0(github_data_dir, apply(folder_hier_for_folders[,1:i], 1, paste0, collapse = "/"))[ndepth_for_folders > i])
  }
})
for(i in 1:(length(folders_to_make_for_folders)-1)){
  for(j in 1:length(folders_to_make_for_folders[[i]])){
    if(!dir.exists(folders_to_make_for_folders[[i]][j])){
      dir.create(folders_to_make_for_folders[[i]][j])    
    }
  }
}

github_file_map_folders <- data.frame(from = multi_files_folders, 
             to = paste0(substr(github_data_dir, 1, nchar(github_data_dir)-1), multi_file_folder_names)
)

#check what already exists in the target directories
file.exists(github_file_map$to)
file.exists(github_file_map_external$to)
file.exists(github_file_map_wildfiles$to)
file.exists(github_file_map_folders$to)

file.exists(github_file_map$from)
file.exists(github_file_map_external$from)
file.exists(github_file_map_wildfiles$from)
file.exists(github_file_map_folders$from)

sum(file.info(github_file_map$from)$size) / 1E9
sum(file.info(github_file_map_external$from)$size) / 1E9
sum(file.info(github_file_map_wildfiles$from)$size) / 1E9
sum(file.info(github_file_map_folders$from)$size) / 1E9

sum(file.info(github_file_map$from)$size[!file.exists(github_file_map$to)]) / 1E9
sum(file.info(github_file_map_external$from)$size[!file.exists(github_file_map_external$to)]) / 1E9
sum(file.info(github_file_map_wildfiles$from)$size[!file.exists(github_file_map_wildfiles$to)]) / 1E9
sum(file.info(github_file_map_folders$from)$size[!file.exists(github_file_map_folders$to)]) / 1E9

#finish file map and write to disk
files_map$new_dir[files_map$new_dir == ""] <- github_data_dir
# write.csv(files_map, file = "~/data/smontgom/pass1b_companion_paper_external_filenames_map.csv")

#### now actually move the files over ####
# file.copy(from = github_file_map$from, to = github_file_map$to)
# file.copy(from = github_file_map_external$from, to = github_file_map_external$to)
# file.copy(from = github_file_map_wildfiles$from, to = github_file_map_wildfiles$to)
# mcopy_directory(from = github_file_map_folders$from, github_file_map_folders$to)

#switch to external drive for this repo
old <- "~/repos/MoTrPAC_Complex_Traits/data/"
new <- "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/"

github_file_map$to <- gsub(old, new, github_file_map$to)
github_file_map_external$to <- gsub(old, new, github_file_map_external$to)
github_file_map_wildfiles$to <- gsub(old, new, github_file_map_wildfiles$to)
github_file_map_folders$to <- gsub(old, new, github_file_map_folders$to)

file.move(from = github_file_map$from, to = github_file_map$to)
file.move(from = github_file_map_external$from, to = github_file_map_external$to)
file.move(from = github_file_map_wildfiles$from, to = github_file_map_wildfiles$to)
dir.move(from = github_file_map_folders$from, github_file_map_folders$to, per_file = F)

#restart earlier move
dir.move(from = github_file_map_folders$from[!file.exists(github_file_map_folders$to)], 
         github_file_map_folders$to[!file.exists(github_file_map_folders$to)], per_file = F)


equal.directories(github_file_map_folders$from[7], github_file_map_folders$to[7])
equal.directories(github_file_map_folders$from[7], github_file_map_folders$to[7], T)

#confirm everything has been transferred
all(file.exists(github_file_map$to),
file.exists(github_file_map_external$to),
file.exists(github_file_map_folders$to),
file.exists(github_file_map_wildfiles$to))

#do some extra cleanup
extra_files <- unlist(read.csv(file = "~/data/smontgom/pass1b_companion_paper_external_extra_filenames.csv")[,2])
not_found <- extra_files[!(file.exists(extra_files) | 
                         sapply(extra_files, function(fi) any(grepl(pattern = glob2rx(fi), x = list.files(dirname(fi), full.names = T)))))]
extra_files <- extra_files[file.exists(extra_files)]

#check if the not_found files have already been transferred
root_dirs <- c("~/data/smontgom/",
               "~/repos/",
               "/Volumes/SSD500GB/")
not_found_roots <- root_dirs[apply(sapply(root_dirs, function(rd) startsWith(not_found, rd)), 1, which)]
transfer_dir <- cbind(external = sapply(seq_along(not_found_roots), function(i){
  poss_filename <- gsub(pattern = not_found_roots[i], 
                        replacement = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/external/", 
                        x = not_found[i])
  file.exists(poss_filename) | any(grepl(pattern = glob2rx(poss_filename), x = list.files(dirname(poss_filename), full.names = T)))
}), 
  internal = sapply(seq_along(not_found_roots), function(i){
    poss_filename <- gsub(pattern = not_found_roots[i], 
                          replacement = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/internal/", 
                          x = not_found[i])
    file.exists(poss_filename) | any(grepl(pattern = glob2rx(poss_filename), x = list.files(dirname(poss_filename), full.names = T)))
  }))
rownames(transfer_dir) <- not_found

not_transferred <- not_found[!apply(transfer_dir, 1, any)]
already_transferred <- setdiff(not_found, not_transferred)
not_transferred

#transfer the extra files remaining
root_dirs <- c("~/data/smontgom/",
               "~/repos/",
               "/Volumes/SSD500GB/")
to_transfer <- setdiff(extra_files, root_dirs)
to_transfer_info <- file.info(to_transfer)
any(to_transfer_info$isdir) #no directories
any(grepl("\\*", to_transfer)) #no wildcards
to_transfer_info$destination <- paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/", c("external/", "internal/"))[
  (to_transfer_info$size / 1E9 < 0.1) + 1]
to_transfer_info$to_transfer_roots <- root_dirs[apply(sapply(root_dirs, function(rd) startsWith(to_transfer, rd)), 1, which)]
to_transfer_info$from <- rownames(to_transfer_info)
rownames(to_transfer_info) <- NULL
to_transfer_info$to <- paste0(to_transfer_info$destination, 
                               sapply(1:nrow(to_transfer_info), function(i)
                                 gsub(to_transfer_info$to_transfer_roots[i], "", to_transfer_info$from[i]))
                               )
file.move(from = to_transfer_info$from, to = to_transfer_info$to)

#write map for other script
extra_map <- data.frame(old_path = rownames(transfer_dir),
                        wild = F, folders = F, single = T,
                        new_dir = paste0("/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/", 
                                         c("external/", "internal/"))[apply(transfer_dir, 1, which)])

extra_map <- rbind(extra_map, 
                   data.frame(old_path = to_transfer_info$from,
                                         wild = F, folders = F, single = T,
                                         new_dir = to_transfer_info$destination),
                   data.frame(old_path = root_dirs,
                              wild = F, folders = F, single = T,
                              new_dir = "/Volumes/2TB_External/MoTrPAC_Complex_Traits/data/"))
write.csv(extra_map, file = "~/data/smontgom/pass1b_companion_paper_external_extra_filenames_map.csv")
