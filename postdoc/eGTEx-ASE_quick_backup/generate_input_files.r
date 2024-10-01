#specifications
setwd("~/repos/egtex-ase")
sample_prop_tail_to_remove <- 0

#### Load Packages ####
library(cmdstanr)
library(parallel)
library(future)
library(posterior)
library(tidyverse)
library(fs)
library(data.table)
library(broom)
library(VGAM)
library(matrixStats)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(S4Vectors)
library(GenomicFeatures)
library(GenomicRanges)

#### Set up future multisession ####
if(!inherits(plan(), "multisession")){
  plan(multisession)  
}

#### specify functions ####
source("scripts/functions.R")
source("~/repos/Stan2R/R/functions.R")

#### create useful variables ####
stan_model_dir <- "Stan/models/"
stan_data_dir <- "Stan/data/"
stan_output_dir <- "Stan/output/"
stan_progress_dir <- "Stan/progress/"
tissue_codes_df <- fread("tissue_codes.txt")
tissue_codes <- setNames(tissue_codes_df$code, tissue_codes_df$name)

#### Compile Input for Bayes (mmPCRSeq) ####

#load in GTEx counts (mmPCRSeq)
tissueColors = fread("gtex_tissue_colors.txt", select = c(1:5))
targetsHG38 = fread("targeted-hg38.bed", col.names = c("chr","pos","pos2","gene","pool","id","ensg"))
targetsHG38$id = paste0("chr",targetsHG38$id)

counts = fread("eGTEx-mmPCR-ASE-aggr.txt.gz")
names(counts)[1:5] = c("tissue","sample","run","chr","pos")

counts$variantID = paste0(counts$chr, "_",
                          counts$pos, "_",
                          counts$refAllele, "_", 
                          counts$altAllele, "_b38")

countsGtex = counts %>% filter(variantID %in% targetsHG38$id & tissue != "OVC")
countsGtex = left_join(countsGtex, tissueColors[,1:2], by=c("tissue" = "tissue_site_detail"))
countsGtex = left_join(countsGtex, targetsHG38[,c(4,6)], by=c("variantID" = "id"))
countsGtex = separate(countsGtex, gene, c("gene","num"),"_")

#load in OVC counts
ovcids = fread("OVCids.txt", col.names = c("stage","grade","individual"))
idCorrection = fread("OVC-sampleIDs-correction.txt", col.names = c("individual","sample"))
ovcids = left_join(ovcids, idCorrection, by="individual")
ocvTestCount = counts %>% filter(variantID %in% targetsHG38$id & 
                                   tissue == "OVC" & 
                                   altCount/totalCount < (1-sample_prop_tail_to_remove) & altCount/totalCount > sample_prop_tail_to_remove)
ocvTestCount = left_join(ocvTestCount, targetsHG38[,c(4,6)], by=c("variantID" = "id"))
ocvTestCount = separate(ocvTestCount, gene, c("gene","num"),"_")
ocvTestCount = left_join(ocvTestCount, ovcids, by = "sample")
ocvTestCount$grade[which(is.na(ocvTestCount$grade) & ocvTestCount$stage == "Benign")] = "Benign,\nUnstaged"
ocvTestCount$grade[which(is.na(ocvTestCount$grade) & is.na(ocvTestCount$stage))] = "Benign,\nUnstaged"
ocvTestCount$grade[which(ocvTestCount$grade == "GB")] = "B"

#split data according to tissue (for model fitting)
tissue_codes <- setNames(tolower(gsub("\\s+", "-", gsub("-", " ", unique(countsGtex$tissue)))), unique(countsGtex$tissue))
tissue_codes_df <- data.frame(code = tissue_codes, name = names(tissue_codes))
fwrite(tissue_codes_df, file = "tissue_codes.txt")
remove_dupe_hetsites <- F
write_tissdat_to_disk <- T
if(write_tissdat_to_disk){
  
  if(!exists("deQTL")){
    deQTL <- fread("countsGtex_lead-eqtl-genotype.csv")
    deQTL <- deQTL[deQTL$genotype == "0/1", c("tissue", "individual", "gene")]
    deQTL <- split(deQTL, deQTL$tissue)
  }
  
  for(tissue_j in unique(countsGtex$tissue)){
    tissue_code <- tissue_codes[tissue_j]
    cat(paste0("(", tissue_code, ") "))
    
    data_GTEx_tiss <- subset(countsGtex, tissue == tissue_j & 
                          altCount/totalCount < (1-sample_prop_tail_to_remove) & altCount/totalCount > sample_prop_tail_to_remove)
    dt <- data_GTEx_tiss[,c("gene", "tissue", "variantID", "refCount", "altCount", "totalCount", "individual")]
    genes_tiss <- unique(dt$gene)
    indivs_tiss <- unique(dt$individual)
    
    #subset to remove pseudoreplicates (same gene, difft het SNP)
    dt$gene_i <- match(dt$gene, genes_tiss)
    dt$indiv_j <- match(dt$individual, indivs_tiss)
    dt$gene_indiv <- paste0(dt$gene, "-", dt$individual)
    
    #combine both pseudoreplicated and dereplicated datasets
    if(remove_dupe_hetsites){
      dt <- split(dt, dt$gene_indiv)
      dt <- do.call(rbind, lapply(dt, function(x){x[which.max(x$totalCount),]}))
    }
    
    #create data object for inference
    dat_tiss <- list(
      n = nrow(dt),
      n_gene = length(genes_tiss),
      n_indiv = length(indivs_tiss),
      gene_i = dt$gene_i,
      indiv_j = dt$indiv_j,
      count = dt$altCount,
      total = dt$totalCount
    )
    
    #write data to disk for later inference
    fwrite(data.frame(index = 1:length(genes_tiss), gene_symbol = genes_tiss), 
           file = paste0(stan_data_dir, tissue_code, "_", sample_prop_tail_to_remove, "_gene-map.csv"))
    fwrite(data.frame(index = 1:length(indivs_tiss), individual = indivs_tiss), 
           file = paste0(stan_data_dir, tissue_code, "_", sample_prop_tail_to_remove, "_indiv-map.csv"))
    fwrite(dt, file = paste0(stan_data_dir, tissue_code, "_", sample_prop_tail_to_remove, ".csv"))
    write_stan_json(dat_tiss, file = paste0(stan_data_dir, tissue_code, "_", sample_prop_tail_to_remove, "_all.json"), 
                    always_decimal = FALSE)
    
    #now split the data by gene and write that to disk
    stan_data_gene_dir <- paste0(stan_data_dir, "genes/")
    base_file_name <- paste0(stan_data_gene_dir, tissue_code, "_", sample_prop_tail_to_remove)
    tiss_deQTL <- deQTL[[tissue_j]]
    
    #write data to disk too for later reference
    for(gene_i in genes_tiss){
      
      file_name <- paste0(base_file_name,  "_", gene_i)  
      
      #subset data to individual gene
      dt_sub <- dt[dt$gene == gene_i,]
      
      #get het indiv\
      het_indivs <- unique(tiss_deQTL$individual[tiss_deQTL$gene == gene_i])
      dt_sub$eQTL_het <- dt_sub$individual %in% het_indivs
      indivs <- unique(dt_sub$individual)
        
      #compile data for Stan
      dat_tiss_gene <- list(
        n = nrow(dt_sub),
        n_indiv = length(unique(indivs)),
        indiv_j = match(dt_sub$individual, indivs),
        count = dt_sub$altCount,
        total = dt_sub$totalCount,
        eQTL_het = indivs %in% het_indivs + 1,
        bounds = c(0.55, 0.95)
      )
      
      fwrite(dt_sub, file = paste0(file_name, ".csv"))
      write_stan_json(dat_tiss_gene, file = paste0(file_name, ".json"), 
                      always_decimal = FALSE)
      
    }
    
  }
}

#subset data for comparison
countsGtex <- countsGtex[countsGtex$tissue == "Ovary",]
data_GTEx <- subset(countsGtex, tissue_site_detail_abbr == "OVARY" & 
                      altCount/totalCount < (1-sample_prop_tail_to_remove) & altCount/totalCount > sample_prop_tail_to_remove)
data_GTEx <- data_GTEx[,c("gene", "tissue", "refCount", "altCount", "totalCount", "individual", "variantID")]
data_GTEx$grade <- "GTEx"
data_OVC <- ocvTestCount[,c("gene", "tissue", "refCount", "altCount", "totalCount", "individual.x", "grade", "variantID")]
colnames(data_OVC)[colnames(data_OVC) == "individual.x"] <- "individual"
genes <- intersect(data_GTEx$gene, data_OVC$gene)

data_OVC <- data_OVC[data_OVC$gene %in% genes,]
data_GTEx <- data_GTEx[data_GTEx$gene %in% genes,]
d <- rbind(data_OVC, data_GTEx)

#sum up gene-indiv obs
d$gene_indiv <- paste0(d$gene, "-", d$individual)

#not doing this because data are not appropriately phased
# dsplit <- split(d, d$gene_indiv)
# dsplit <- lapply(dsplit, function(x){
#   out <- x[1,]
#   if(nrow(x) > 1){
#     out[,c("refCount","altCount","totalCount")] <- apply(x[,c("refCount","altCount","totalCount")], 2, sum)
#   }
#   return(out)
# })
# d <- do.call(rbind, dsplit)

d$gene_i <- match(d$gene, genes)
indivs <- unique(d$individual)
fwrite(x = data.frame(index = 1:length(indivs), gene_symbol = indivs), 
       file = paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove, "_indiv-map.csv"))
d$indiv_j <- match(d$individual, indivs)
fwrite(x = data.frame(index = 1:length(genes), gene_symbol = genes), 
       file = paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove, "_gene-map.csv"))
d$gtex <- 2 - (d$grade == "GTEx") #so = 1 if in GTEx, = 2 if cancer
grades <- setNames(object = c(1, 2, 2, 2, 3, 4, 5), 
                   nm = c("GTEx", "Benign,\nUnstaged", "B", "X", "I", "II", "III"))
# grades <- setNames(object = c(1:7), 
#                    nm = c("GTEx", "Benign,\nUnstaged", "B", "X", "I", "II", "III"))
d$grade_k <- grades[d$grade]

#subset d to remove pseudoreplicates (same gene, diff't het SNP)
remove_dupe_hetsites <- F
if(remove_dupe_hetsites){
  d <- split(d, d$gene_indiv)
  hist(1 - sapply(d, function(x){max(x$totalCount) / sum(x$totalCount)}), xlab = "prop obs lost from discarding hets")
  d <- do.call(rbind, lapply(d, function(x){x[which.max(x$totalCount),]}))
}

well_sampled_genes <- F
if(well_sampled_genes){
  min_n_per_grade <- 3
  d_split_genes <- split(d[,c("indiv_j", "grade_k")], d$gene_i)
  n_grade_per_gene <- do.call(rbind, lapply(d_split_genes, function(x) 
    table(factor(x$grade_k[!duplicated(x)], levels = 1:5))
    ))
  n_min_per_grade <- apply(n_grade_per_gene, 1, min)
  table(n_min_per_grade)
  passing_gene_i <- as.numeric(which(n_min_per_grade >= min_n_per_grade))
  mean(d$gene_i %in% passing_gene_i)
  d <- d[d$gene_i %in% passing_gene_i,]
  genes <- unique(d$gene)
  indivs <- unique(d$individual)
  d$gene_i <- match(d$gene, genes)
  d$indiv_j <- match(d$individual, indivs)
}


#write to disk
fwrite(x = data.frame(index = 1:length(indivs), gene_symbol = indivs), 
       file = paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove,
                     ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                     "_indiv-map.csv"))
fwrite(x = data.frame(index = 1:length(genes), gene_symbol = genes), 
       file = paste0(stan_data_dir, "ovc-vs-gtex_", sample_prop_tail_to_remove,
                     ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                     "_gene-map.csv"))

dat <- list(
  n = nrow(d),
  n_gene = length(genes),
  n_indiv = length(indivs),
  n_grade = length(unique(d$grade_k)),
  gene_i = d$gene_i,
  indiv_j = d$indiv_j,
  grade_k = d$grade_k,
  indiv_x_gene_l = match(d$gene_indiv, unique(d$gene_indiv)),
  gtex = d$gtex,
  count = d$altCount,
  total = d$totalCount
)
rm(counts)

#save the data to disk
write_stan_json(dat, file = paste0(stan_data_dir, "OVC-vs-GTEx", "_", 
                                   sample_prop_tail_to_remove,
                                   ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                                   ifelse(remove_dupe_hetsites, "", "_all"), ".json"), 
                always_decimal = FALSE)
fwrite(d, file = paste0(stan_data_dir, "OVC-vs-GTEx", "_", sample_prop_tail_to_remove,
                        ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                        ".csv"))

#write individual gene data to disk
stan_data_genes_dir <- paste0(stan_data_dir, "genes/")
for(gene_i in genes){
  cat("(", gene_i, ") ")
  gd <- d[d$gene == gene_i,]
  gindivs <- unique(gd$individual)
  gd$indiv_j <- match(gd$individual, gindivs)
  indiv_grade <- data.frame(gd$indiv_j, gd$grade_k)
  indiv_grade <- indiv_grade[!duplicated(indiv_grade),]
  gdat <- list(
    n = nrow(gd),
    n_indiv = length(gindivs),
    n_grade = length(unique(d$grade_k)),
    indiv_j = gd$indiv_j,
    grade_k = indiv_grade$gd.grade_k,
    ovc = gd$gtex,
    count = gd$altCount,
    total = gd$totalCount,
    bounds = c(0.55, 0.95)
  )
  gene_path <- paste0(stan_data_genes_dir, "OVC-vs-GTEx", "_", 
                      sample_prop_tail_to_remove,
                      ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""),
                      ifelse(remove_dupe_hetsites, "", "_all"), 
                      "_", gene_i)
  write_stan_json(gdat, file = paste0(gene_path, ".json"), 
                  always_decimal = FALSE)
  fwrite(gd, file = paste0(gene_path, ".csv"))
}

#now permute count and write that to disk
dat_perm <- dat
dat_perm_inds <- sample(1:dat_perm$n)
dat_perm$count <- dat_perm$count[dat_perm_inds]
dat_perm$total <- dat_perm$total[dat_perm_inds]
write_stan_json(dat_perm, file = paste0(stan_data_dir, "OVC-vs-GTEx_perm", 
                                        ifelse(remove_dupe_hetsites, "", "_all"), ".json"), 
                always_decimal = FALSE)


#### Compile Input for Bayes (RNASeq) ####
rnaTargets <- dir_ls("data/rnaseq-targets/",glob = "*")
names(rnaTargets) <- gsub("data/rnaseq-targets/|.v8.wasp_corrected.ase_table.tsv","",names(rnaTargets))
rnaCounts <- rnaTargets %>% map_dfr(fread, .id = "individual")
rnaCounts <- rnaCounts[rnaCounts$GENOTYPE_WARNING == 0 &
                         rnaCounts$LOW_MAPABILITY == 0 &
                         rnaCounts$MAPPING_BIAS_SIM == 0,]

tissueColors = fread("gtex_tissue_colors.txt", select = c(1:5))
tissue_codes_df <- fread(file = "tissue_codes.txt")
tissue_codes <- setNames(tissue_codes_df$code, tissue_codes_df$name)
rnaCounts$TISSUE_NAME <- setNames(tissueColors$tissue_site_detail, 
                                  tissueColors$tissue_site_detail_abbr)[
                                    rnaCounts$TISSUE_ID]
rnaCounts$tissue_code <- tissue_codes[rnaCounts$TISSUE_NAME]
rnaCounts <- rnaCounts[,c("individual", "TISSUE_ID", "REF_COUNT", 
                          "ALT_COUNT", "TOTAL_COUNT", "TISSUE_NAME", 
                          "tissue_code", "VARIANT_ID", "GENE_ID")]

#map from ensembl IDs to gene symbols
targetsHG38 <- fread("targeted-hg38.bed", col.names = c("chr","pos","pos2","gene","pool","id","ensg"))
targetsHG38$id <- paste0("chr",targetsHG38$id)
targetsHG38$ensembl <- do.call(rbind, strsplit(targetsHG38$ensg, "\\."))[,1]
targetsHG38$gene <- do.call(rbind, strsplit(targetsHG38$gene, "_"))[,1]
rnaCounts <- rnaCounts[rnaCounts$GENE_ID %in% targetsHG38$ensembl,]
rnaCounts$gene <- setNames(targetsHG38$gene, targetsHG38$ensembl)[rnaCounts$GENE_ID]

#medulla omitted previously for some reason (no color code I guess)
rnaCounts$TISSUE_NAME[rnaCounts$TISSUE_ID == "KDNMDL"] <- "Kidney - Medulla"
rnaCounts$tissue_code[rnaCounts$TISSUE_ID == "KDNMDL"] <- "kidney-medulla"

colnames(rnaCounts) <- c("individual" = "individual", 
                         "TISSUE_ID" = "tissue_id", 
                         "REF_COUNT" = "refCount", 
                         "ALT_COUNT" = "altCount", 
                         "TOTAL_COUNT" = "totalCount", 
                         "TISSUE_NAME" = "tissue",
                         "VARIANT_ID" = "variantID",
                         "gene" = "gene",
                         "GENE_ID" = "ensembl",
                         "tissue_code" = "tissue_codes")[colnames(rnaCounts)]

stan_data_dir_rnaseq <- paste0(stan_data_dir, "rnaseq/")
remove_dupe_hetsites <- F
write_tissdat_to_disk <- F
if(write_tissdat_to_disk){
  
  if(!exists("deQTL")){
    deQTL <- fread("countsGtex_lead-eqtl-genotype.csv")
    deQTL <- deQTL[deQTL$genotype == "0/1", c("tissue", "individual", "gene")]
    deQTL <- split(deQTL, deQTL$tissue)
  }
  
  for(tissue_j in unique(rnaCounts$tissue)){
    tissue_code <- tissue_codes[tissue_j]
    cat(paste0("(", tissue_code, ") "))
    
    data_GTEx_tiss <- subset(rnaCounts, tissue == tissue_j & 
                               altCount/totalCount < (1-sample_prop_tail_to_remove) & altCount/totalCount > sample_prop_tail_to_remove)
    dt <- data_GTEx_tiss[,c("gene", "tissue", "variantID", "refCount", "altCount", "totalCount", "individual")]
    genes_tiss <- unique(dt$gene)
    indivs_tiss <- unique(dt$individual)
    
    #subset to remove pseudoreplicates (same gene, difft het SNP)
    dt$gene_i <- match(dt$gene, genes_tiss)
    dt$indiv_j <- match(dt$individual, indivs_tiss)
    dt$gene_indiv <- paste0(dt$gene, "-", dt$individual)
    
    #dereplicate dataset if necessary
    if(remove_dupe_hetsites){
      dt <- split(dt, dt$gene_indiv)
      dt <- do.call(rbind, lapply(dt, function(x){x[which.max(x$totalCount),]}))
    }
    
    #create data object for inference
    dat_tiss <- list(
      n = nrow(dt),
      n_gene = length(genes_tiss),
      n_indiv = length(indivs_tiss),
      gene_i = dt$gene_i,
      indiv_j = dt$indiv_j,
      count = dt$altCount,
      total = dt$totalCount
    )
    
    #write data to disk for later inference
    fwrite(data.frame(index = 1:length(genes_tiss), gene_symbol = genes_tiss), 
           file = paste0(stan_data_dir_rnaseq, tissue_code, "_", sample_prop_tail_to_remove, "_gene-map.csv"))
    fwrite(data.frame(index = 1:length(indivs_tiss), individual = indivs_tiss), 
           file = paste0(stan_data_dir_rnaseq, tissue_code, "_", sample_prop_tail_to_remove, "_indiv-map.csv"))
    fwrite(dt, file = paste0(stan_data_dir_rnaseq, tissue_code, "_", sample_prop_tail_to_remove, ".csv"))
    write_stan_json(dat_tiss, file = paste0(stan_data_dir_rnaseq, tissue_code, "_", sample_prop_tail_to_remove, "_all.json"), 
                    always_decimal = FALSE)
    
    #now split the data by gene and write that to disk
    stan_data_gene_dir <- paste0(stan_data_dir_rnaseq, "genes/")
    base_file_name <- paste0(stan_data_gene_dir, tissue_code, "_", sample_prop_tail_to_remove)
    tiss_deQTL <- deQTL[[tissue_j]]
    
    #write data to disk too for later reference
    for(gene_i in genes_tiss){
      
      file_name <- paste0(base_file_name,  "_", gene_i)  
      
      #subset data to individual gene
      dt_sub <- dt[dt$gene == gene_i,]
      
      #get het indiv\
      het_indivs <- unique(tiss_deQTL$individual[tiss_deQTL$gene == gene_i])
      dt_sub$eQTL_het <- dt_sub$individual %in% het_indivs
      indivs <- unique(dt_sub$individual)
      
      #compile data for Stan
      dat_tiss_gene <- list(
        n = nrow(dt_sub),
        n_indiv = length(unique(indivs)),
        indiv_j = match(dt_sub$individual, indivs),
        count = dt_sub$altCount,
        total = dt_sub$totalCount,
        eQTL_het = indivs %in% het_indivs + 1,
        bounds = c(0.55, 0.95)
      )
      
      fwrite(dt_sub, file = paste0(file_name, ".csv"))
      write_stan_json(dat_tiss_gene, file = paste0(file_name, ".json"), 
                      always_decimal = FALSE)
      
    }
    
  }
}

#now x references to cancer

#first do GTEx
rnaCountsGTEx <- rnaCounts[rnaCounts$tissue == "Ovary",]
data_GTEx <- subset(rnaCountsGTEx, tissue_id == "OVARY" & 
                      altCount/totalCount < (1-sample_prop_tail_to_remove) & 
                      altCount/totalCount > sample_prop_tail_to_remove)
data_GTEx <- data_GTEx[,c("variantID", "gene", "tissue", "refCount", "altCount", "totalCount", "individual")]
data_GTEx$grade <- "GTEx"

#now OVC
ovc_countsRNA <- fread("RNAseq_OVBM_ase_counts.tsv")
ovc_countsRNA$individual <- ovc_countsRNA$OTB_Ids
ovc_countsRNA$Grade_Differentiation[is.na(ovc_countsRNA$Grade_Differentiation)] <- "Benign,\nUnstaged"
ovc_countsRNA$grade <- ovc_countsRNA$Grade_Differentiation
ovc_countsRNA$tissue <- "OVC"
ovc_countsRNA$refCount <- ovc_countsRNA$AD1
ovc_countsRNA$altCount <- ovc_countsRNA$AD2
ovc_countsRNA$totalCount <- ovc_countsRNA$sum
ovc_countsRNA$variantID <- paste0(ovc_countsRNA$CHROM, "_",
                                  ovc_countsRNA$POS, "_",
                                  ovc_countsRNA$REF, "_", 
                                  ovc_countsRNA$ALT, "_b38")
mean(ovc_countsRNA$variantID %in% data_GTEx$variantID)
mean(data_GTEx$variantID %in% ovc_countsRNA$variantID)
common_variants <- intersect(ovc_countsRNA$variantID, data_GTEx$variantID)

#subset data for comparison
ovc_countsRNA <- ovc_countsRNA[ovc_countsRNA$variantID %in% common_variants,]
data_GTEx <- data_GTEx[data_GTEx$variantID %in% common_variants,]
ovc_countsRNA$gene <- setNames(data_GTEx$gene, data_GTEx$variantID)[ovc_countsRNA$variantID]
data_OVC <- ovc_countsRNA[,c("gene", "tissue", "refCount", "altCount", "totalCount", "individual", "grade", "variantID")]
data_GTEx <- data_GTEx[,c("gene", "tissue", "refCount", "altCount", "totalCount", "individual", "grade", "variantID")]

#get gene map
genes <- intersect(data_GTEx$gene, data_OVC$gene)
data_OVC <- data_OVC[data_OVC$gene %in% genes,]
data_GTEx <- data_GTEx[data_GTEx$gene %in% genes,]
d <- rbind(data_OVC, data_GTEx)

#sum up gene-indiv obs
d$gene_indiv <- paste0(d$gene, "-", d$individual)
d$gene_i <- match(d$gene, genes)
indivs <- unique(d$individual)
fwrite(x = data.frame(index = 1:length(indivs), gene_symbol = indivs), 
       file = paste0(stan_data_dir_rnaseq, "ovc-vs-gtex_", sample_prop_tail_to_remove, "_indiv-map_rnaseq.csv"))
d$indiv_j <- match(d$individual, indivs)
fwrite(x = data.frame(index = 1:length(genes), gene_symbol = genes), 
       file = paste0(stan_data_dir_rnaseq, "ovc-vs-gtex_", sample_prop_tail_to_remove, "_gene-map_rnaseq.csv"))
d$gtex <- 2 - (d$grade == "GTEx") #so = 1 if in GTEx, = 2 if cancer
d$grade[which(d$grade == "GB")] <- "B"

grades <- setNames(object = c(1, 2, 2, 2, 2, 3, 4, 5), 
                   nm = c("GTEx", "Benign,\nUnstaged", "GB", "B", "X", "I", "II", "III"))
# grades <- setNames(object = c(1:7), 
#                    nm = c("GTEx", "Benign,\nUnstaged", "B", "X", "I", "II", "III"))
d$grade_k <- grades[d$grade]

#subset d to remove pseudoreplicates (same gene, diff't het SNP)
remove_dupe_hetsites <- F
if(remove_dupe_hetsites){
  d <- split(d, d$gene_indiv)
  hist(1 - sapply(d, function(x){max(x$totalCount) / sum(x$totalCount)}), xlab = "prop obs lost from discarding hets")
  d <- do.call(rbind, lapply(d, function(x){x[which.max(x$totalCount),]}))
}

dat <- list(
  n = nrow(d),
  n_gene = length(genes),
  n_indiv = length(indivs),
  n_grade = length(unique(d$grade_k)),
  gene_i = d$gene_i,
  indiv_j = d$indiv_j,
  grade_k = d$grade_k,
  indiv_x_gene_l = match(d$gene_indiv, unique(d$gene_indiv)),
  gtex = d$gtex,
  count = d$altCount,
  total = d$totalCount
)

#save the data to disk
write_stan_json(dat, file = paste0(stan_data_dir_rnaseq, "OVC-vs-GTEx", "_", 
                                   sample_prop_tail_to_remove, 
                                   ifelse(remove_dupe_hetsites, "", "_all"), "_rnaseq.json"), 
                always_decimal = FALSE)


#write individual gene data to disk
stan_data_dir_rnaseq <- paste0(stan_data_dir, "rnaseq/")
stan_data_gene_dir_rnaseq <- paste0(stan_data_dir_rnaseq, "genes/")
for(gene_i in genes){
  cat("(", gene_i, ") ")
  gd <- d[d$gene == gene_i,]
  gindivs <- unique(gd$individual)
  gd$indiv_j <- match(gd$individual, gindivs)
  indiv_grade <- data.frame(gd$indiv_j, gd$grade_k)
  indiv_grade <- indiv_grade[!duplicated(indiv_grade),]
  gdat <- list(
    n = nrow(gd),
    n_indiv = length(gindivs),
    n_grade = length(unique(d$grade_k)),
    indiv_j = gd$indiv_j,
    grade_k = indiv_grade$gd.grade_k,
    ovc = gd$gtex,
    count = gd$altCount,
    total = gd$totalCount,
    bounds = c(0.55, 0.95)
  )
  gene_path <- paste0(stan_data_gene_dir_rnaseq, "OVC-vs-GTEx", "_", 
                      sample_prop_tail_to_remove,
                      ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""),
                      ifelse(remove_dupe_hetsites, "", "_all"), 
                      "_", gene_i)
  write_stan_json(gdat, file = paste0(gene_path, "_rnaseq.json"), 
                  always_decimal = FALSE)
  fwrite(gd, file = paste0(gene_path, "_rnaseq.csv"))
}



#### Hyperprior for beta prob plot ####
nsamps <- 2E4
conc_prior <- c(mean = 3, sd = 1)
loc_prior <- c(mean = 0, sd = 1)
mix_prob <- runif(nsamps)
point_mass_width <- 0.02
concs <- exp(rnorm(nsamps, mean = conc_prior["mean"], 
                   sd = conc_prior["sd"]))
locs <- 0.5 + invlogit(rnorm(nsamps, mean = loc_prior["mean"], 
                             sd = loc_prior["sd"])) / 2
alphas <- concs * locs
betas <- concs * (1-locs)
prob_range <- 0:500/500
prob_range <- prob_range - (prob_range-median(prob_range))/1E2
denss <- lapply(1:nsamps, function(i)
  dbeta(prob_range, alphas[i], betas[i]) * (1-mix_prob[i]))

beta_mix_dens <- apply(do.call(rbind, denss), 2, mean)
beta_mix_dens <- beta_mix_dens / sum(beta_mix_dens) / diff(prob_range[1:2]) * (1-median(mix_prob))
beta_mix_dens <- (beta_mix_dens + rev(beta_mix_dens)) / 2
point_mass_dens <- mix_prob / point_mass_width

# trim <- function(x, n) x[(n+1):(length(x)-n)]
plot(prob_range, beta_mix_dens, type = "l",
     ylim = c(0, max(beta_mix_dens) * 2.25),
     xlab = "ASE Binomial Probability",
     ylab = "Density",
     main = latex2exp::TeX("E(Density) of the Intercept Prior ($\\alpha$)"))
polygon(x = c(prob_range, rev(prob_range)),
        y = c(beta_mix_dens, rep(0, length(beta_mix_dens))),
        col = adjustcolor(1, 0.2), border = adjustcolor(1, 0.1))
segments(0.5, 0, 0.5, max(beta_mix_dens) * 1.75, lwd = 2)

plot(NULL, NULL, xlim = c(0,1), ylim = c(0,10) , xlab = "ASE Binomial Probability",
     ylab = "Density")
for(i in 1:min(c(nsamps, 500))){
  polygon(x = c(prob_range, rev(prob_range)),
          y = (c(denss[[i]], rep(0, length(prob_range))) +
                 c(rev(denss[[i]]), rep(0, length(prob_range))))/2,
          col = adjustcolor(1, 0.005), border = adjustcolor(1, 0.005))
}
for(i in 1:100){
  segments(0.5, 0, 0.5, runif(1, 0, par("usr")[4]), lwd = 4,
           col = adjustcolor(1, 0.005))
}

#### Annotated PDF plot ####
grDevices::cairo_pdf(filename = paste0("~/repos/egtex-ase/figures/ase-model-priors.pdf"), 
                     width = 1200 / 72, height = 400 / 72, family="Arial Unicode MS", pointsize = 19)
par(mfrow = c(1,2), mar = c(5,5,3,2))

#samples from the prior
plot(NULL, NULL, xlim = c(0,1), ylim = c(0, max(beta_mix_dens) * 10),
     xlab = "ASE Binomial Probability",
     ylab = "Density",
     main = latex2exp::TeX("Samples from the Intercept Prior ($\\alpha$)"))
for(i in 1:min(c(nsamps, 500))){
  polygon(x = c(prob_range, rev(prob_range)),
          y = (c(denss[[i]], rep(0, length(prob_range))) +
                 c(rev(denss[[i]]), rep(0, length(prob_range))))/2,
          col = adjustcolor(1, 0.005), border = adjustcolor(1, 0.05))
}
for(i in 1:500){
  segments(0.5, 0, 0.5, runif(1, 0, par("usr")[4]), lwd = 4,
           col = adjustcolor(1, 0.005))
}

#average density
plot(prob_range, beta_mix_dens, type = "l",
     ylim = c(0, max(beta_mix_dens) * 2.5),
     xlab = "ASE Binomial Probability",
     ylab = "Density",
     main = latex2exp::TeX("E(Density) of the Intercept Prior ($\\alpha$)"))
polygon(x = c(prob_range, rev(prob_range)),
        y = c(beta_mix_dens, rep(0, length(beta_mix_dens))),
        col = adjustcolor(1, 0.2), border = adjustcolor(1, 0.1))
arrows(0.5, 0, 0.5, max(beta_mix_dens) * 1.8, lwd = 3, length = 0.1)


cols_annot <- c(pi = rgb(140, 21, 21, maxColorValue = 255), 
                mu = rgb(21, 140, 21, maxColorValue = 255), 
                nu = rgb(21, 21, 140, maxColorValue = 255))

#annotations for mu
arrowlwd <- 3
arrows(x0 = prob_range[which.max(beta_mix_dens)], 
       x1 = prob_range[which.max(beta_mix_dens)] , 
       y0 = max(beta_mix_dens) * 2, 
       y1 = max(beta_mix_dens) * 1.1, lwd = arrowlwd, 
       col = cols_annot["mu"], length = 0.1)
arrows(x0 = 1-prob_range[which.max(beta_mix_dens)], 
       x1 = 1-prob_range[which.max(beta_mix_dens)] , 
       y0 = max(beta_mix_dens) * 2, 
       y1 = max(beta_mix_dens) * 1.1, lwd = arrowlwd, 
       col = cols_annot["mu"], length = 0.1)
segments(x0 = prob_range[which.max(beta_mix_dens)], 
         x1 = 1-prob_range[which.max(beta_mix_dens)] , 
         y0 = max(beta_mix_dens) * 2, 
         y1 = max(beta_mix_dens) * 2, lwd = arrowlwd, 
         col = "white", length = 0.1)
segments(x0 = prob_range[which.max(beta_mix_dens)], 
         x1 = 1-prob_range[which.max(beta_mix_dens)] , 
         y0 = max(beta_mix_dens) * 2, 
         y1 = max(beta_mix_dens) * 2, lwd = arrowlwd, 
         col = cols_annot["mu"], length = 0.1)
segments(x0 = 0.5, 
         x1 = 0.5, 
         y0 = max(beta_mix_dens) * 2, 
         y1 = max(beta_mix_dens) * 2.1, lwd = arrowlwd, 
         col = cols_annot["mu"])
text(x = 0.5, y = max(beta_mix_dens) * 2.1, col = adjustcolor(cols_annot["mu"], 1), pos = 3,
     labels = latex2exp::TeX("$\\mu$"), xpd = NA, cex = 1.25)

#annotations for vu
arrows(x0 = prob_range[which.max(beta_mix_dens)] - 0.1, 
       x1 = prob_range[which.max(beta_mix_dens)] + 0.1, 
       y0 = max(beta_mix_dens) * 0.5, 
       y1 = max(beta_mix_dens) * 0.5, lwd = arrowlwd, 
       col = cols_annot["nu"], code = 3, length = 0.1)
arrows(x0 = 1 - prob_range[which.max(beta_mix_dens)] - 0.1, 
       x1 = 1 - prob_range[which.max(beta_mix_dens)] + 0.1, 
       y0 = max(beta_mix_dens) * 0.5, 
       y1 = max(beta_mix_dens) * 0.5, lwd = arrowlwd, 
       col = cols_annot["nu"], code = 3, length = 0.1)
segments(x0 = 1 - prob_range[which.max(beta_mix_dens)], 
         x1 = 1 - prob_range[which.max(beta_mix_dens)], 
         y0 = max(beta_mix_dens) * 0.5, 
         y1 = max(beta_mix_dens) * 0.4, lwd = arrowlwd, 
         col = cols_annot["nu"])
text(x = 1 - prob_range[which.max(beta_mix_dens)], 
     y = max(beta_mix_dens) * 0.4, 
     col = adjustcolor(cols_annot["nu"], 1), pos = 1,
     labels = latex2exp::TeX("$\\nu$"), 
     xpd = NA, cex = 1.25)
segments(x0 = prob_range[which.max(beta_mix_dens)], 
         x1 = prob_range[which.max(beta_mix_dens)], 
         y0 = max(beta_mix_dens) * 0.5, 
         y1 = max(beta_mix_dens) * 0.4, lwd = arrowlwd, 
         col = cols_annot["nu"])
text(x = prob_range[which.max(beta_mix_dens)], 
     y = max(beta_mix_dens) * 0.4, 
     col = adjustcolor(cols_annot["nu"], 1), pos = 1,
     labels = latex2exp::TeX("$\\nu$"), 
     xpd = NA, cex = 1.25)

#annotations for pi
arrows(x0 = 0.5325, 
       x1 = 0.5325, 
       y0 = max(beta_mix_dens) * 1.6, 
       y1 = max(beta_mix_dens) * 1, lwd = arrowlwd, 
       col = cols_annot["pi"], code = 3, length = 0.1)
segments(x0 = 0.5325, 
         x1 = 0.555, 
         y0 = max(beta_mix_dens) * 1.3, 
         y1 = max(beta_mix_dens) * 1.3, lwd = arrowlwd, 
         col = cols_annot["pi"])
text(x = 0.555, 
     y = max(beta_mix_dens) * 1.3, 
     col = adjustcolor(cols_annot["pi"], 1), pos = 4,
     labels = latex2exp::TeX("$\\pi$"), 
     xpd = NA, cex = 1.25)


dev.off()



#### Beta ASE plot ####
nsamps <- 5E3
conc_prior <- c(mean = 3, sd = 1)
loc_prior <- c(mean = 1, sd = 0.5)
mix_prob <- invlogit(rnorm(nsamps)-1)
point_mass_width <- 0.02
concs <- exp(rnorm(nsamps, mean = conc_prior["mean"], 
                   sd = conc_prior["sd"]))
locs <- 0.5 + invlogit(rnorm(nsamps, mean = loc_prior["mean"], 
                             sd = loc_prior["sd"])) / 2
alphas <- concs * locs
betas <- concs * (1-locs)
prob_range <- 0:500/500
prob_range <- prob_range - (prob_range-median(prob_range))/1E2
denss <- lapply(1:nsamps, function(i)
  dbeta(prob_range, alphas[i], betas[i]) * (1-mix_prob[i]))

beta_mix_dens <- apply(do.call(rbind, denss), 2, mean)
beta_mix_dens <- beta_mix_dens / sum(beta_mix_dens) / diff(prob_range[1:2]) * (1-median(mix_prob))
beta_mix_dens <- (beta_mix_dens + rev(beta_mix_dens)) / 2
point_mass_dens <- mix_prob / point_mass_width


ase_conc_prior <- c(mean = 5, sd = 1)
ase_concs <- exp(rnorm(nsamps, mean = ase_conc_prior["mean"], 
                   sd = conc_prior["sd"]))
ase_locs <- rep(0.5, nsamps)
ase_alphas <- ase_concs * ase_locs
ase_betas <- ase_concs * (1-ase_locs)
ase_denss <- lapply(1:nsamps, function(i)
  dbeta(prob_range, ase_alphas[i], ase_betas[i]) * (mix_prob[i]))

ase_beta_mix_dens <- apply(do.call(rbind, ase_denss), 2, mean)
ase_beta_mix_dens <- ase_beta_mix_dens / sum(ase_beta_mix_dens) / diff(prob_range[1:2]) * (median(mix_prob))
ase_beta_mix_dens <- (ase_beta_mix_dens + rev(ase_beta_mix_dens)) / 2



# trim <- function(x, n) x[(n+1):(length(x)-n)]
plot(prob_range, beta_mix_dens, type = "l",
     ylim = c(0, max(c(beta_mix_dens, ase_beta_mix_dens)) * 2.25),
     xlab = "ASE Binomial Probability",
     ylab = "Density",
     main = latex2exp::TeX("E(Density) of the Intercept Prior ($\\alpha$)"))
polygon(x = c(prob_range, rev(prob_range)),
        y = c(beta_mix_dens, rep(0, length(beta_mix_dens))),
        col = adjustcolor(1, 0.2), border = adjustcolor(1, 0.1))

lines(prob_range, ase_beta_mix_dens)
polygon(x = c(prob_range, rev(prob_range)),
        y = c(ase_beta_mix_dens, rep(0, length(ase_beta_mix_dens))),
        col = adjustcolor(1, 0.2), border = adjustcolor(1, 0.1))


plot(NULL, NULL, xlim = c(0,1), ylim = c(0,10) , xlab = "ASE Binomial Probability",
     ylab = "Density")
for(i in 1:min(c(nsamps, 500))){
  polygon(x = c(prob_range, rev(prob_range)),
          y = (c(denss[[i]], rep(0, length(prob_range))) +
                 c(rev(denss[[i]]), rep(0, length(prob_range))))/2,
          col = adjustcolor(1, 0.005), border = adjustcolor(1, 0.005))
}
for(i in 1:min(c(nsamps, 500))){
  polygon(x = c(prob_range, rev(prob_range)),
          y = (c(ase_denss[[i]], rep(0, length(prob_range))) +
                 c(rev(ase_denss[[i]]), rep(0, length(prob_range))))/2,
          col = adjustcolor(1, 0.005), border = adjustcolor(1, 0.005))
}

#### Annotated PDF plot ####
grDevices::cairo_pdf(filename = paste0("~/repos/egtex-ase/figures/ase-model-priors.pdf"), 
                     width = 1200 / 72, height = 400 / 72, family="Arial Unicode MS", pointsize = 19)
par(mfrow = c(1,2), mar = c(5,5,3,2))

#samples from the prior
plot(NULL, NULL, xlim = c(0,1), ylim = c(0,max(c(beta_mix_dens, ase_beta_mix_dens)) * 5) , xlab = "ASE Binomial Probability",
     ylab = "Density", main = latex2exp::TeX("Samples from the Intercept Prior ($\\alpha$)"))
for(i in 1:min(c(nsamps, 500))){
  polygon(x = c(prob_range, rev(prob_range)),
          y = (c(denss[[i]], rep(0, length(prob_range))) +
                 c(rev(denss[[i]]), rep(0, length(prob_range))))/2,
          col = adjustcolor(1, 0.005), border = adjustcolor(1, 0.05))
}
for(i in 1:min(c(nsamps, 500))){
  polygon(x = c(prob_range, rev(prob_range)),
          y = (c(ase_denss[[i]], rep(0, length(prob_range))) +
                 c(rev(ase_denss[[i]]), rep(0, length(prob_range))))/2,
          col = adjustcolor(1, 0.005), border = adjustcolor(1, 0.05))
}

#average density
plot(prob_range, beta_mix_dens, type = "l",
     ylim = c(0, max(c(beta_mix_dens, ase_beta_mix_dens)) * 1.5),
     xlab = "ASE Binomial Probability",
     ylab = "Density",
     main = latex2exp::TeX("E(Density) of the Intercept Prior ($\\alpha$)"))
polygon(x = c(prob_range, rev(prob_range)),
        y = c(beta_mix_dens, rep(0, length(beta_mix_dens))),
        col = adjustcolor(1, 0.2), border = adjustcolor(1, 0.1))

lines(prob_range, ase_beta_mix_dens)
polygon(x = c(prob_range, rev(prob_range)),
        y = c(ase_beta_mix_dens, rep(0, length(ase_beta_mix_dens))),
        col = adjustcolor(1, 0.2), border = adjustcolor(1, 0.1))

cols_annot <- c(pi = rgb(140, 21, 21, maxColorValue = 255), 
                mu = rgb(21, 140, 21, maxColorValue = 255), 
                nu = rgb(21, 21, 140, maxColorValue = 255))

#annotations for mu
arrowlwd <- 3
arrows(x0 = prob_range[which.max(beta_mix_dens)], 
       x1 = prob_range[which.max(beta_mix_dens)] , 
       y0 = max(beta_mix_dens) * 2, 
       y1 = max(beta_mix_dens) * 1.1, lwd = arrowlwd, 
       col = cols_annot["mu"], length = 0.1)
arrows(x0 = 1-prob_range[which.max(beta_mix_dens)], 
       x1 = 1-prob_range[which.max(beta_mix_dens)] , 
       y0 = max(beta_mix_dens) * 2, 
       y1 = max(beta_mix_dens) * 1.1, lwd = arrowlwd, 
       col = cols_annot["mu"], length = 0.1)
segments(x0 = prob_range[which.max(beta_mix_dens)], 
         x1 = 1-prob_range[which.max(beta_mix_dens)] , 
         y0 = max(beta_mix_dens) * 2, 
         y1 = max(beta_mix_dens) * 2, lwd = arrowlwd, 
         col = "white")
segments(x0 = prob_range[which.max(beta_mix_dens)], 
         x1 = 1-prob_range[which.max(beta_mix_dens)] , 
         y0 = max(beta_mix_dens) * 2, 
         y1 = max(beta_mix_dens) * 2, lwd = arrowlwd, 
         col = cols_annot["mu"])
segments(x0 = 0.5, 
         x1 = 0.5, 
         y0 = max(beta_mix_dens) * 2, 
         y1 = max(beta_mix_dens) * 2.1, lwd = arrowlwd, 
         col = cols_annot["mu"])
text(x = 0.5, y = max(beta_mix_dens) * 2.1, col = adjustcolor(cols_annot["mu"], 1), pos = 3,
     labels = latex2exp::TeX("$\\mu$"), xpd = NA, cex = 1.25)

#annotations for vu
vu_arrow_w <- 0.1
arrows(x0 = prob_range[which.max(beta_mix_dens)] - vu_arrow_w/2, 
       x1 = prob_range[which.max(beta_mix_dens)] + vu_arrow_w/2, 
       y0 = max(beta_mix_dens) * 0.5, 
       y1 = max(beta_mix_dens) * 0.5, lwd = arrowlwd, 
       col = cols_annot["nu"], code = 3, length = 0.1)
arrows(x0 = 1 - prob_range[which.max(beta_mix_dens)] - vu_arrow_w/2, 
       x1 = 1 - prob_range[which.max(beta_mix_dens)] + vu_arrow_w/2, 
       y0 = max(beta_mix_dens) * 0.5, 
       y1 = max(beta_mix_dens) * 0.5, lwd = arrowlwd, 
       col = cols_annot["nu"], code = 3, length = 0.1)
segments(x0 = 1 - prob_range[which.max(beta_mix_dens)], 
         x1 = 1 - prob_range[which.max(beta_mix_dens)], 
         y0 = max(beta_mix_dens) * 0.5, 
         y1 = max(beta_mix_dens) * 0.4, lwd = arrowlwd, 
         col = cols_annot["nu"])
text(x = 1 - prob_range[which.max(beta_mix_dens)], 
     y = max(beta_mix_dens) * 0.4, 
     col = adjustcolor(cols_annot["nu"], 1), pos = 1,
     labels = latex2exp::TeX("$\\nu$"), 
     xpd = NA, cex = 1.25)
segments(x0 = prob_range[which.max(beta_mix_dens)], 
         x1 = prob_range[which.max(beta_mix_dens)], 
         y0 = max(beta_mix_dens) * 0.5, 
         y1 = max(beta_mix_dens) * 0.4, lwd = arrowlwd, 
         col = cols_annot["nu"])
text(x = prob_range[which.max(beta_mix_dens)], 
     y = max(beta_mix_dens) * 0.4, 
     col = adjustcolor(cols_annot["nu"], 1), pos = 1,
     labels = latex2exp::TeX("$\\nu$"), 
     xpd = NA, cex = 1.25)

#annotations for pi
arrows(x0 = 0.5, 
       x1 = 0.5, 
       y0 = max(beta_mix_dens) * 1.3, 
       y1 = max(beta_mix_dens) * 0.8, lwd = arrowlwd, 
       col = cols_annot["pi"], code = 3, length = 0.1)
segments(x0 = 0.5, 
         x1 = 0.555, 
         y0 = max(beta_mix_dens) * 1.05, 
         y1 = max(beta_mix_dens) * 1.05, lwd = arrowlwd, 
         col = cols_annot["pi"])
text(x = 0.555, 
     y = max(beta_mix_dens) * 1.05, 
     col = adjustcolor(cols_annot["pi"], 1), pos = 4,
     labels = latex2exp::TeX("$\\pi$"), 
     xpd = NA, cex = 1.25)

dev.off()

#### plot prior predictive ####

model_index <- 8
model_names <- list(
  "ase-gtex_inf-priors_no-indiv", #1
  "ase-gtex_inf-priors", #2
  "ase-gtex", #3
  "ase-gtex_centered", #4
  "ase-gtex_offset-conc", #5
  "ase-gtex_nearly-balanced", #6
  "ase-gtex_nearly-balanced_LoH", #7
  "ase-gtex_nearly-balanced_LoH_nested", #8
  "ase-gtex_nearly-balanced_LoH_noIndiv", #9
  "ase-gtex_nearly-balanced_LoH_nested_noIndiv", #10
  "ase-gtex_nearly-balanced_LoH_indivXgene", #11
  "ase-gtex_nearly-balanced_LoH_nested_indivXgene", #12
  "ase-gtex_nearly-balanced_LoH_nested_centered", #13
  "ase-gtex_nearly-balanced_LoH_nested_conc-scale-bounds", #14
  
  "ase-gtex_nearly-balanced_LoH_nested_sample-probs", #15
  "ase-gtex_nearly-balanced_LoH_nested_sample-probs-truncnorm" #16
)
model_name <- model_names[[model_index]]
model_name <- "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-less_inform_conc"
sample_prop_tail_to_remove <- 0

model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- cmdstan_model(model_path)
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
tiss_j <- "artery-aorta"
dat <- jsonlite::fromJSON(paste0(stan_data_dir, tiss_j, "_", sample_prop_tail_to_remove, "_all.json"))


n_prior_samps <- 200
nobs <- dat$n
nsubsamp <- 1E3

prior_samps <- mclapply(1:n_prior_samps, function(i){
  
  mcprint(paste0(i, " "))
  
  r_code <- parse_Stan(model_lines, 
                       dat = dat, 
                       samps = NA, 
                       output_file = NA, 
                       sample_index = NA, 
                       post_pred_sim = F, 
                       sim = TRUE)
  
  #evaluate code
  stan_env <- new.env()
  eval(parse(text = r_code), envir = stan_env)
  prior_predictive_sim <- stan_env$out
  
  #sample from appropriate beta distributions for model 8
  
  probs_slab <- rbeta(nobs, prior_predictive_sim$shape_slab_alpha, prior_predictive_sim$shape_slab_beta)
  # probs_slab <- c(probs_slab, 1-probs_slab)
  probs_loh <- rbeta(nobs, prior_predictive_sim$shape_loh_alpha, prior_predictive_sim$shape_loh_beta)
  # probs_loh <- c(probs_loh, 1-probs_loh)
  probs_ab <- rbeta(nobs, prior_predictive_sim$shapes_ab, prior_predictive_sim$shapes_ab)
  # probs_ab <- c(probs_ab, 1-probs_ab) #not necessary (prior already symmetric) but it's a cheap way to increase n
  three_probs <- cbind(probs_ab, probs_loh, probs_slab)
  
  dist_picked <- apply(prior_predictive_sim$mix_prob, 1, function(x) sample(1:3, 1, prob = c(x[1], x[2] + x[3], x[4] + x[5])))
  obs_subsamp_inds <- sort(sample(1:nobs, size = nsubsamp, replace = F))
  probs_picked <- abs(sample(c(0,1), nsubsamp, replace = T) - three_probs[cbind(obs_subsamp_inds, dist_picked[obs_subsamp_inds])])
  
  return(cbind(dist_picked[obs_subsamp_inds], probs_picked))
  
}, mc.cores = 8)

prior_samps <- do.call(rbind, prior_samps)
prior_samps_split <- split(prior_samps[,2], f = prior_samps[,1])
breaks <- c(0:50)/50
hist(prior_samps_split[[1]], breaks = breaks)
hist(prior_samps_split[[2]], breaks = breaks)
hist(prior_samps_split[[3]], breaks = breaks)

hist(unlist(prior_samps_split), breaks = breaks, freq = F)
plot(density(unlist(prior_samps_split)))

density_bounded_logit <- function(x, bounds = c(0,1)){
  #useful functions
  logit <- function(p) log(p/(1-p))
  logit_prime <- function(p) (1/p + 1 / (1-p))
  invlogit <- function(x) {exp(x) / (1+exp(x))}
  
  #rescale to (0,1)
  xs <- x - bounds[1]
  xs <- xs / bounds[2]
  
  #transform to unbounded space
  xub <- logit(xs)
  
  #compute kde
  logit_dens <- density(xub)
  
  #transform back to (0,1) space
  invlogit_dens <- logit_dens
  invlogit_dens$x <- invlogit(logit_dens$x)
  invlogit_dens$y <- logit_dens$y * logit_prime(invlogit_dens$x)
  
  #transform back to (a,b) space
  invlogit_dens$x <- invlogit_dens$x * bounds[2] + bounds[1]
  invlogit_dens$y <- invlogit_dens$y / bounds[2]
  
  return(invlogit_dens)
  
}


density_bounded_reflect <- function(x, bounds = c(0,1)){
  
  #transform to reflected space
  xr <- c(x, bounds[1] - x, bounds[2] + (bounds[2] - x))
  
  #compute kde
  refl_dens <- density(xr, from = bounds[1], to = bounds[2])
  
  return(refl_dens)
  
}


plot(density_bounded_reflect(unlist(prior_samps_split)))

bounded_kdes <- lapply(prior_samps_split, density_bounded_reflect)
xlocs <- bounded_kdes[[i]]$x
ydens <- NULL
for(i in seq_along(bounded_kdes)){
  bounded_kdes[[i]]$y <- bounded_kdes[[i]]$y * length(prior_samps_split[[i]]) / sum(sapply(prior_samps_split, length))
  ydens <- cbind(ydens, bounded_kdes[[i]]$y)
}

border_cols <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
fill_cols <- adjustcolor(border_cols, 0.5)
ifelse2 <- function(test, yes, no) if(test){return(yes)}else{return(no)}
plot_multidens <- function(x, y, fill_cols = "lightgrey", border_cols = "black", stacked = T, ...){
  ny <- ncol(y)
  if(length(border_cols) < ny){
    border_cols <- rep(border_cols, length.out = ny)
  }
  if(length(fill_cols) < ny){
    fill_cols <- rep(fill_cols, length.out = ny)
  }
  
  if(stacked){
    csum_y <- t(apply(y[,c(3,2,1)], 1, cumsum))
    plot(x, csum_y[,ny], type = "l", col = border_cols[ny], ylim = c(0, max(csum_y[,ny])), ...)
    for(i in 1:ny){
      polygon(x = c(x, rev(x)), y = c(csum_y[,i], ifelse2(i==1, rep(0,length(x)), csum_y[,i-1])), 
              col = fill_cols[i], border = border_cols[i], lwd = 2)  
      
    }  
  } else {
    plot(x, y[,ny], type = "l", col = border_cols[ny], ylim = c(0, max(y)), ...)
    for(i in 1:ny){
      polygon(x = c(x, rev(x)), y = c(y[,i], rep(0,length(x))), 
              col = fill_cols[i], border = border_cols[i], lwd = 2)
    }
  }
}


plot_multidens(x = xlocs, y = ydens[,c(3,2,1)], fill_cols = fill_cols, border_cols = border_cols,
               xlab = "draw probability", ylab = "density", stacked = T)
legend(x = "topright", col=  border_cols[1:3], pt.cex = 2, box.col = "black", bg = "grey97", pt.bg = fill_cols[1:3], 
       legend = c("Allelic Imbalance", "Loss of Heterozygosity", "Allelic Balance"), pch = 22)
lines(xlocs, csum_ydens[,2], type = "l")
lines(xlocs, csum_ydens[,1], type = "l")
polygon()

# apply(prior_predictive_sim$mix_prob, 2, mean)

# 1 = probability in allelic balance
# 2 = probability not in AB but in LoH (tail 1)
# 3 = probability not in AB but in LoH (tail 2)
# 4 = probability not in AB but in slab (tail 1)
# 5 = probability not in AB but in slab (tail 2)

#### prior predictive for OVC vs GTEx ####

use_rnaseq <- F
if(use_rnaseq){
  stan_data_dir <- "Stan/data/rnaseq/"
} else {
  stan_data_dir <- "Stan/data/"
}
model_index <- 25
model_names <- list(
  "beta-binomial-joint-conc", #1
  "beta-binomial-joint-conc-loc", #2
  "beta-binomial-joint-conc-loc-noncentered", #3
  "beta-binomial-joint-conc-loc-noncentered-simple", #4, does not fit
  "beta-binomial-joint-loc", #5
  "beta-binomial-joint-conc-loc-indiv", #6
  "beta-binomial-simple-mixture", #7
  "beta-binomial-simplest-mixture", #8
  "beta-binomial-simpler-mixture", #9
  "beta-binomial-joint-conc-loc-indiv-pointmix", #10
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered", #11
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_grade", #12
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene", #13
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid", #14
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid-identifiability", #15
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_iid_nearly-balanced", #16
  "beta-binomial-joint-conc-loc-indiv-pointmix-noncentered-tumor_gradeXgene_cumulative_nearly-balanced", #17
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade", #18
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc", #19
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc-correct_grade_REfs", #20
  "ovc-ase-gtex_IID_nearly-balanced_LoH_nested-tumor_grade", #21
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc_wider-init", #22
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc_IID-Diff-Split", #23
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-inform_conc_wider-init_cumul-add", #24
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-less_inform_conc", #25
  "ovc-ase-gtex_nearly-balanced_LoH_nested-tumor_grade-less_inform_conc-no_indiv" #26
)

model_name <- model_names[[model_index]]
model_path <- paste0(stan_model_dir, model_name, ".stan")
mod <- cmdstan_model(model_path, cpp_options = list(STAN_THREADS = T, stan_threads = TRUE))
model_lines <- readLines(model_path)
model_string <- paste0(model_lines, collapse = "\n")
cat(model_path)
# cat(model_string)

#load in the data
data_index <- 1
well_sampled_genes <- F
remove_dupe_hetsites <- F
min_n_per_grade <- 3
data_names <- list(
  "OVC-vs-GTEx_0", #1
  "OVC-vs-GTEx_perm_0", #2
  "OVC-vs-GTEx_sim_0", #3
  "OVC-vs-GTEx_0_all" #4
)
data_name <- paste0(data_names[[data_index]],
                    ifelse(well_sampled_genes, paste0("_min-", min_n_per_grade, "-per-grade"), ""), 
                    ifelse(use_rnaseq, "_rnaseq", ""),
                    ifelse(remove_dupe_hetsites, "", "_all"))
dat <- jsonlite::fromJSON(paste0(stan_data_dir, data_name, ".json"))

n_prior_samps <- 100
nobs <- dat$n
nsubsamp <- 1E3

prior_samps <- mclapply(1:n_prior_samps, function(i){
  
  mcprint(paste0(i, " "))
  
  r_code <- parse_Stan(stan_code = model_lines, 
                       dat = dat, 
                       samps = NA, 
                       output_file = NA, 
                       sample_index = NA, 
                       post_pred_sim = F, 
                       sim = TRUE)
  # open_script(r_code)
  
  #evaluate code
  stan_env <- new.env()
  eval(parse(text = r_code), envir = stan_env)
  prior_predictive_sim <- stan_env$out
  
  #sample from appropriate beta distributions for model 25
  five_probs <- prior_predictive_sim$mix_prob
  three_probs <- data.frame(AB = five_probs[,1], 
                   loh = five_probs[,2] * 2, 
                   slab = five_probs[,4] * 2)
  beta_probs <- data.frame(AB = rbeta(n = nobs, shape1 = prior_predictive_sim$shapes_ab, 
                                      shape2 = prior_predictive_sim$shapes_ab),
                           loh = rbeta(n = nobs, shape1 = prior_predictive_sim$shape_loh_alpha, 
                                      shape2 = prior_predictive_sim$shape_loh_beta),
                           slab = rbeta(n = nobs, shape1 = prior_predictive_sim$shape_slab_alpha, 
                                      shape2 = prior_predictive_sim$shape_slab_beta))
  slab_flip <- sample(c(T,F), nobs, T)
  loh_flip <- sample(c(T,F), nobs, T)
  beta_probs$slab[slab_flip] <- 1 - beta_probs$slab[slab_flip]
  beta_probs$loh[loh_flip] <- 1 - beta_probs$loh[loh_flip]
  
  sampled_component <- apply(three_probs, 1, function(p) sample(1:3, 1, prob = p))
  sampled_probability <- beta_probs[cbind(x = 1:nobs, y = sampled_component)]
  subsamp_inds <- sort(sample(1:nobs, nsubsamp))
  
  return(data.frame(datind = subsamp_inds, 
                    mixcomp = sampled_component[subsamp_inds], 
                    betaprob = sampled_probability[subsamp_inds]))
  
}, mc.cores = 2)


prior_samps <- do.call(rbind, prior_samps)
prior_samps_split <- split(prior_samps$betaprob, f = prior_samps$mixcomp)
breaks <- c(0:50)/50
hist(prior_samps_split[[1]], breaks = breaks)
hist(prior_samps_split[[2]], breaks = breaks)
hist(prior_samps_split[[3]], breaks = breaks)

hist(unlist(prior_samps_split), breaks = breaks, freq = F)
plot(density(unlist(prior_samps_split)))

density_bounded_logit <- function(x, bounds = c(0,1)){
  #useful functions
  logit <- function(p) log(p/(1-p))
  logit_prime <- function(p) (1/p + 1 / (1-p))
  invlogit <- function(x) {exp(x) / (1+exp(x))}
  
  #rescale to (0,1)
  xs <- x - bounds[1]
  xs <- xs / bounds[2]
  
  #transform to unbounded space
  xub <- logit(xs)
  
  #compute kde
  logit_dens <- density(xub)
  
  #transform back to (0,1) space
  invlogit_dens <- logit_dens
  invlogit_dens$x <- invlogit(logit_dens$x)
  invlogit_dens$y <- logit_dens$y * logit_prime(invlogit_dens$x)
  
  #transform back to (a,b) space
  invlogit_dens$x <- invlogit_dens$x * bounds[2] + bounds[1]
  invlogit_dens$y <- invlogit_dens$y / bounds[2]
  
  return(invlogit_dens)
  
}


density_bounded_reflect <- function(x, bounds = c(0,1)){
  
  #transform to reflected space
  xr <- c(x, bounds[1] - x, bounds[2] + (bounds[2] - x))
  
  #compute kde
  refl_dens <- density(xr, from = bounds[1], to = bounds[2])
  
  return(refl_dens)
  
}


plot(density_bounded_reflect(unlist(prior_samps_split)))

bounded_kdes <- lapply(prior_samps_split, density_bounded_reflect)
xlocs <- bounded_kdes[[i]]$x
ydens <- NULL
for(i in seq_along(bounded_kdes)){
  bounded_kdes[[i]]$y <- bounded_kdes[[i]]$y * length(prior_samps_split[[i]]) / sum(sapply(prior_samps_split, length))
  ydens <- cbind(ydens, bounded_kdes[[i]]$y)
}

border_cols <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
fill_cols <- adjustcolor(border_cols, 0.5)
ifelse2 <- function(test, yes, no) if(test){return(yes)}else{return(no)}
plot_multidens <- function(x, y, fill_cols = "lightgrey", border_cols = "black", stacked = T, ...){
  ny <- ncol(y)
  if(length(border_cols) < ny){
    border_cols <- rep(border_cols, length.out = ny)
  }
  if(length(fill_cols) < ny){
    fill_cols <- rep(fill_cols, length.out = ny)
  }
  
  if(stacked){
    csum_y <- t(apply(y[,c(3,2,1)], 1, cumsum))
    plot(x, csum_y[,ny], type = "l", col = border_cols[ny], ylim = c(0, max(csum_y[,ny])), ...)
    for(i in 1:ny){
      polygon(x = c(x, rev(x)), y = c(csum_y[,i], ifelse2(i==1, rep(0,length(x)), csum_y[,i-1])), 
              col = fill_cols[i], border = border_cols[i], lwd = 2)  
      
    }  
  } else {
    plot(x, y[,ny], type = "l", col = border_cols[ny], ylim = c(0, max(y)), ...)
    for(i in 1:ny){
      polygon(x = c(x, rev(x)), y = c(y[,i], rep(0,length(x))), 
              col = fill_cols[i], border = border_cols[i], lwd = 2)
    }
  }
}


####plot_multihist####

plot_multihist <- function(input, fill_cols = "lightgrey", 
                           border_cols = "black", stacked = T, nbins = 100,
                           xlab = "Value", ylab = "Mass", lab.cex = 1, label_blobs = F,
                           ...){
  
  #get colors settled
  nhist <- length(input)
  if(length(border_cols) < nhist){
    border_cols <- rep(border_cols, length.out = nhist)
  }
  if(length(fill_cols) < nhist){
    fill_cols <- rep(fill_cols, length.out = nhist)
  }
  
  #process histograms input
  nbreaks <- nbins + 1
  input_range <- range(unlist(input))
  breaks <- seq(input_range[1], input_range[2], length.out = nbreaks)
  hists_raw <- do.call(cbind, lapply(input, function(x) 
    hist(x, breaks = breaks, plot = F)$density))
  
  #add up across rows if hists are stacked
  if(stacked){
    hists_mass <- t(apply(hists_raw, 1, cumsum))
  } else {
    hists_mass <- hists_raw
  }
  
  #get polygon coordinates
  hists <- data.frame(lb = breaks[-(nbreaks)], 
                      ub = breaks[-1],
                      dens = hists_mass)
  eps <- 1E-6
  hist_polys <- lapply(1:nhist, function(i){
    
    hmas <- hists_raw[,i]
    hcmas <- hists_mass[,i]
    rlehmas <- rle(hmas > 1E-6)
    csum_rlehmas <- cumsum(rlehmas$lengths)
    hchunks <- lapply(1:length(rlehmas$lengths), function(chi){
      if(rlehmas$values[chi]){
        if(chi == 1){
          return(1:rlehmas$lengths[1])
        } else {
          return((csum_rlehmas[chi-1]+1):(csum_rlehmas[chi]))
        }
      } else {
        return(NULL)
      }
    })
    hchunks <- hchunks[!sapply(hchunks, is.null)]
    nchunks <- length(hchunks)
    
    poly_coords <- lapply(1:nchunks, function(chi){
      bini <- hchunks[[chi]]
      nbins_chunk <- length(bini)
      brini <- c(bini, tail(bini,1) + 1)
      if(i == 1 || !stacked){
        hx = c(breaks[brini], 
               rep(rev(breaks[brini]), each = 2)[-1])
        hx = hx[-length(hx)]
        hy = c(rep(0, nbins_chunk+1), 
               rep(rev(hists_mass[hchunks[[chi]],i]), each = 2))
      } else {
        hx = c(rep(breaks[brini], each = 2)[-1], 
               rep(rev(breaks[brini])[-1], each = 2))
        hy = c(rep(hists_mass[bini,i-1], each = 2), 
               rep(rev(hists_mass[bini,i]), each = 2),
               hists_mass[bini[1],i-1])
        length(hx)
        length(hy)
      }
      return(data.frame(x = hx, y = hy ))  
    })
    
    return(poly_coords)
    
  })
  
  #start plotting
  plot.new()
  plot.window(xlim = range(c(0, breaks)), ylim = range(c(0, hists_raw)))
  
  #plot histograms
  for(i in 1:nhist){
    hp <- hist_polys[[i]]
    nchunks <- length(hp)
    for(j in 1:nchunks){
      polygon(hp[[j]]$x, 
              hp[[j]]$y, 
              col = fill_cols[i],
              border = border_cols[i])  
    }
    
  }
  
  #label blobs
  if(!is.null(names(input)) & label_blobs){
    for(i in 1:nhist){
      cathist <- do.call(rbind, hist_polys[[i]])
      labx <- mean(cathist$x[which.max(cathist$y) + 0:1])
      dispy <- diff(par("usr")[3:4]) / 20
      laby <- max(cathist$y)
      text(x = labx, y = laby + dispy,
           labels = names(input)[i],
           col = border_cols[i], xpd = NA)  
      segments(x0 = labx, x1 = labx, y0 = laby, 
               y1 = laby + dispy - strheight(names(input)[i], units = "user") / 2,
               col = border_cols[i], xpd = NA)
    }
  }

  #histogram axes
  axis(1, at = pretty(range(c(0, breaks))), 
       labels = pretty(range(c(0, breaks))), line = 0.1)
  axis(2, line = 0.1)
  
  #axis labels
  mtext(side = 1, 
        text = xlab, 
        line = 2.5, cex = lab.cex)
  mtext(side = 2, 
        text = ylab, 
        line = 2.5, cex = lab.cex)
 
}

border_cols <- c("darkred", "darkorchid4", "darkblue")
fill_cols <- adjustcolor(border_cols, 0.5)
plot_multihist(prior_samps_split[3:1], 
               fill_cols = fill_cols, border_cols = border_cols, 
               nbins = 50)


#### Het-site Consistency ####

#start with just one tissue
tissue_code <- tissue_codes[1]
sample_prop_tail_to_remove <- 0
dt <- fread(file = paste0(stan_data_dir, tissue_code, "_", sample_prop_tail_to_remove, ".csv"))
dt_removed <- fread(file = paste0(stan_data_dir, tissue_code, "_", sample_prop_tail_to_remove, "_removed.csv"))
dt_dupes <- rbind(dt[dt$gene_indiv %in% dt_removed$gene_indiv,], dt_removed)

dt_dupes$count_hdev <- abs(dt_dupes$refCount - dt_dupes$totalCount / 2)
dt_dupes$prop_hdev <- dt_dupes$count_hdev / dt_dupes$totalCount * 2
dt_dupes$pmean_prop_hdev <- abs((dt_dupes$refCount+1) / (dt_dupes$totalCount+2) - 0.5) / 0.5
dt_dupes$q95_prop_hdev <- qbeta_dev(p = 0.95, shape1 = dt_dupes$refCount+1, shape2 = dt_dupes$altCount+1) / 0.5

dt_dupesplit <- split(dt_dupes, dt_dupes$gene_indiv)
dupemat <- do.call(rbind, lapply(dt_dupesplit, function(x){
  y <- as.data.frame(x[,c("gene", "individual", "totalCount", "prop_hdev", "q95_prop_hdev", "pmean_prop_hdev", "variantID")])
  combins <- combn(1:nrow(y), m = 2)
  out <- data.frame(individual = y$individual[1], 
                    gene = y$gene[1],
                    s1 = y[combins[1,],"variantID"],
                    s2 = y[combins[2,],"variantID"],
                    p1 = y[combins[1,],"prop_hdev"],
                    p2 = y[combins[2,],"prop_hdev"],
                    pmp1 = y[combins[1,],"pmean_prop_hdev"],
                    pmp2 = y[combins[2,],"pmean_prop_hdev"],
                    pq95p1 = y[combins[1,],"q95_prop_hdev"],
                    pq95p2 = y[combins[2,],"q95_prop_hdev"],
                    n1 = y[combins[1,],"totalCount"],
                    n2 = y[combins[2,],"totalCount"])
  return(out)
        
}))
dupemat$min_n <- apply(dupemat[,c("n1", "n2")], 1, min)
dupemat$max_p <- apply(dupemat[,c("pq95p1", "pq95p2")], 1, max)

pt_col_inds <- ceiling(log10(dupemat$min_n) / max(log10(dupemat$min_n)) * 100)
pt_col_inds <- pt_col_inds - min(pt_col_inds) + 1
col_scale <- (viridisLite::viridis(max(pt_col_inds)))

par(mar = c(5,5,4,7))
plot(dupemat[,c("pmp1", "pmp2")], pch = 19, col = adjustcolor(col_scale[pt_col_inds], 0.2), main = dt_dupes$tissue[1],
     xlab = "allelic balance at het site #1", ylab = "allelic balance at het site #2")
add_continuous_legend(col_scale, x = 1.075, y = 1,
                      labels = dupemat$min_n, 
                      positions = (pt_col_inds - min(pt_col_inds)) / max(pt_col_inds - min(pt_col_inds)),
                      log_scale = T, left_below = F, n_labels = 5, main = "n")

#correlations within pairs go down as you increase sample size?
prob_threshes <- 0:40/100
corr_prop_thresh <- sapply(prob_threshes, function(prop_thresh) 
  ccc(logit(dupemat[(1-dupemat$max_p) > prop_thresh, c("pq95p1", "pq95p2")])))
plot(logit(dupemat[,c("pq95p1", "pq95p2")]), pch = 19, col = adjustcolor(col_scale[pt_col_inds], 0.2),
     xlab = "logit(allelic balance at het site #1)", ylab = "logit(allelic balance at het site #2)")
rect(par("usr")[1]-diff(par("usr")[1:2]), par("usr")[3]-diff(par("usr")[3:4]), logit(0.9), logit(0.9), border = 2, lwd = 2, lty = 2)
text(x = par("usr")[1] + diff(par("usr")[1:2]) / 3, 
     y = par("usr")[4] - strheight("= exclusion threshold"), labels = "= exclusion threshold", col = 2, pos = 4, font = 2)
segments(x0 = par("usr")[1] + diff(par("usr")[1:2]) / 40, 
         x1 = par("usr")[1] + diff(par("usr")[1:2]) / 3, 
         y0 = par("usr")[4] - strheight("= exclusion threshold"), 
         y1 = par("usr")[4] - strheight("= exclusion threshold"), col = 2, lwd = 2, lty = 2)
plot(prob_threshes, corr_prop_thresh, type = "l", lwd = 2,
     xlab = "tail exclusion proportion", ylab = "concordance correlation coefficient", col = 2)

#split by gene or individual?
dupemat_g <- split(dupemat, dupemat$gene)
by_gene_ccc <- sapply(dupemat_g, function(x) ccc(logit(x[,c("pq95p1", "pq95p2")])))
by_gene_pc <- sapply(dupemat_g, function(x) cor(logit(x[,c("pq95p1", "pq95p2")]))[1,2])

par(mfrow = c(2,2), mar = c(0,4,8,1))
bgccc_hist <- hist(by_gene_ccc, breaks = seq(-1,1,length.out = 50), plot = F)
barplot(bgccc_hist$counts, axes = TRUE, space = 0, 
        ylab = "frequency")
par(mar = c(0,0,8,8))
plot.new()
text(x = 0.5, y = 0.5, paste0("Within-gene distribution of paired correlations,\n(within-individual, multiple het sites),\n", dt_dupes$tissue[1]), xpd = NA)
par(mar = c(4,4,1,1))
plot(by_gene_ccc, by_gene_pc, xlab = "concordance correlation coefficient", ylab = "pearson's correlation coefficient",
     col = adjustcolor(1, 0.8), pch = 19, xlim = c(-1,1), ylim = c(-1,1))
par(mar = c(4,0,1,8))
bgpc_hist <- hist(by_gene_pc, breaks = seq(-1,1,length.out = 50), plot = F)
barplot(bgpc_hist$counts, axes = TRUE, space = 0, horiz=TRUE, 
        xlab = "frequency")


dupemat_i <- split(dupemat, dupemat$individual)
by_indiv_ccc <- sapply(dupemat_i, function(x) ccc(x[,c("pq95p1", "pq95p2")]))
by_indiv_pc <- sapply(dupemat_i, function(x) cor(x[,c("pq95p1", "pq95p2")])[1,2])
par(mfrow = c(2,2), mar = c(0,4,8,1))
biccc_hist <- hist(by_indiv_ccc, breaks = seq(-1,1,length.out = 50), plot = F)
barplot(biccc_hist$counts, axes = TRUE, space = 0, 
        ylab = "frequency")
par(mar = c(0,0,8,8))
plot.new()
text(x = 0.5, y = 0.5, paste0("Within-individual distribution of paired correlations,\n(within-individual, multiple het sites),\n", dt_dupes$tissue[1]), xpd = NA)
par(mar = c(4,4,1,1))
plot(by_indiv_ccc, by_indiv_pc, xlab = "concordance correlation coefficient", ylab = "pearson's correlation coefficient",
     col = adjustcolor(1, 0.8), pch = 19, xlim = c(-1,1), ylim = c(-1,1))
par(mar = c(4,0,1,8))
bipc_hist <- hist(by_indiv_pc, breaks = seq(-1,1,length.out = 50), plot = F)
barplot(bipc_hist$counts, axes = TRUE, space = 0, horiz=TRUE, 
        xlab = "frequency")

#### Het-site annotation ####

# Load or download the GRCh38 genome annotations
#this all uses GTEx v8, which used GENCODE 26 https://www.gtexportal.org/home/methods, which matches ENSEMBL 98 https://www.gencodegenes.org/human/releases.html
txdb <- makeTxDbFromEnsembl(release=88, organism="Homo sapiens")
all_genes <- genes(txdb)
gene_info <- subset(all_genes, mcols(all_genes)$gene_id %in% targetsHG38$ensembl)
targetsHG38$ensembl <- do.call(rbind, strsplit(targetsHG38$ensg, "\\."))[,1]

# Loop through each row in targetsHG38
exons <- exonsBy(txdb, by="gene")
exons_by_transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)
transcripts <- transcriptsBy(txdb, by="gene")
GTEx_metadata <- fread("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
GTEx_isoforms <- fread("filtered_transcript_isopct.gct")
GTEx_isoforms$transcript_id <- gsub("\\..*", "", GTEx_isoforms$transcript_id)
GTEx_isoforms$gene_id <- gsub("\\..*", "", GTEx_isoforms$gene_id)
GTEx_isoforms <- melt(GTEx_isoforms, id.vars = c("transcript_id", "gene_id"), variable.name = "sample_id", value.name = "expression")
GTEx_isoforms_indivs <- sapply(lapply(strsplit(x = colnames(GTEx_isoforms)[-(1:2)], "-"), head, n=2), paste0, collapse = "-")
GTEx_isoforms_tissues <-  GTEx_metadata$SMTSD[match(colnames(GTEx_isoforms)[-(1:2)], GTEx_metadata$SAMPID)]
# split(GTEx_isoforms, GTEx_isoforms)

locus_info <- parallel::mclapply(1:nrow(targetsHG38), function(i) {
  
  mcprint(i)
  
  # Extract the current Ensembl ID and position
  current_ensembl <- targetsHG38$ensembl[i]
  current_pos <- targetsHG38$pos[i]
  
  # Find the gene in all_genes
  gene <- subset(all_genes, mcols(all_genes)$gene_id %in% current_ensembl)
  gene_id <- gene$gene_id
  
  # Calculate start, end, and proportion to 3' end
  gene_start <- start(gene)
  if(length(gene_start) == 0){return(NULL)}
  gene_end <- end(gene)
  gene_length <- gene_end - gene_start
  proportion_to_3_end <- ifelse(strand(gene) == "+", 
                                (current_pos - gene_start) / gene_length, 
                                (gene_end - current_pos) / gene_length)
  nloci_to_3_end <- ifelse(strand(gene) == "+", 
                                (current_pos - gene_start), 
                                (gene_end - current_pos))
  
  #TODO look at hetsites in the same isoform for consistency, not in the same gene
  GTEx_isoforms_xGene <- (GTEx_isoforms[GTEx_isoforms$gene_id == gene_id,])
  colnames(GTEx_isoforms_xGene)
  # Check if the position is within any exon of the gene
  target_range <- GRanges(seqnames=seqnames(gene), ranges=IRanges(start=current_pos, end=current_pos))
  in_exon <- any(subjectHits(findOverlaps(target_range, exons[[current_ensembl]])))
  
  # and if it is in a variable exon
  transcripts_for_gene <- transcripts[[current_ensembl]]
  mean(transcripts_for_gene$tx_name %in% names(exons_by_transcripts))
  
  # Extract transcript names for the current gene
  transcript_names <- transcripts[[current_ensembl]]$tx_name
  
  #TODO x-reference against expressed transcripts in GTEx?
  #https://gtexportal.org/home/downloads/adult-gtex/overview
  
  transcripts_with_site <- sapply(seq_along(transcript_names), function(j){
    # Extract exons for the current transcript
    current_exons <- exons_by_transcripts[[transcript_names[j]]]
    
    # Check for overlap between the target site and the current exons
    overlaps <- findOverlaps(target_range, current_exons)
    
    # Record if there is at least one overlap
    return(length(overlaps) > 0)
  })
  
  # Determine if the site is in a variable exon:
  # If the site is not present in the exons of all transcripts, it's in a variable exon
  in_variable_exon <- !all(transcripts_with_site)
  prop_transcripts <- mean(transcripts_with_site)
  
  # put the results in a df
  results_df <- data.frame(ensembl=current_ensembl, 
                           start=gene_start, 
                           end=gene_end, 
                           proportion_to_3_end=proportion_to_3_end,
                           nloci_to_3_end=nloci_to_3_end,
                           in_exon=in_exon,
                           in_variable_exon=in_variable_exon,
                           prop_transcripts=prop_transcripts,
                           id=targetsHG38$id[i],
                           chr=targetsHG38$chr[i],
                           pos=targetsHG38$pos[i],
                           symbol=strsplit(targetsHG38$gene[i], "_")[[1]][1])
  
  #and return
  return(results_df)
}, mc.cores = 15)
locus_info <- do.call(rbind, locus_info)
mean(locus_info$in_variable_exon)
hist(locus_info$proportion_to_3_end, breaks = -1:100/100, xlab = "relative distance to 3' end", main = "all het sites in study")
hist(log10(locus_info$nloci_to_3_end), xlab = latex2exp::TeX("log$_{10}$(absolute distance to 3' end in bp)"), 
     main = "all het sites in study", breaks = 100, axes = F)
axis(2); axis(1, at = 0:6, labels = 10^(0:6))
hist(locus_info$prop_transcripts, breaks = 0:100/100, xlab = "proportion presence in gene's transcripts", main = "all het sites in study")

# 3' Bias in Poly(A) Selection: Many RNA-seq library preparation methods include a poly(A) selection step, which targets and enriches mRNA by binding to its poly(A) tail at the 3' end. This can lead to a 3' bias, where reads closer to the 3' end of transcripts are more likely to be captured and sequenced. Consequently, these reads might be over-represented in the data.
# RNA Degradation: RNA molecules can degrade over time, with certain regions possibly being more susceptible to degradation than others. In general, RNA degradation tends to start from the ends of the molecules. However, because of the poly(A) tail protection and potentially more efficient capture of the 3' end regions during library preparation, reads near the 3' end might still be more stable or more likely to be detected in some contexts.

#### Monoallelic Hets GTEx Depth ####

rnaTargets = dir_ls("data/rnaseq-targets/",glob = "*")
names(rnaTargets) = gsub("data/rnaseq-targets/|.v8.wasp_corrected.ase_table.tsv","",names(rnaTargets))
rnaCounts = rnaTargets %>% map_dfr(fread, .id = "individual")
rnaCounts_x_tiss <- split(rnaCounts, rnaCounts$TISSUE_ID)
tiss_abbrev <- "ADPSBQ Adipose-Subcutaneous, ADPVSC Adipose-Visceral(Omentum), ADRNLG Adrenal Gland, ARTAORT Artery-Aorta, ARTCRN Artery-Coronary, ARTTBL Artery-Tibial, BLDDER Bladder, BRNAMY Brain-Amygdala, BRNACC Brain-Anterior cingulate cortex (BA24), BRNCDT Brain-Caudate(basal ganglia), BRNCHB Brain-Cerebellar Hemisphere, BRNCHA Brain-Cerebellum, BRNCTXA Brain-Cortex, BRNCTXB Brain-Frontal Cortex (BA9), BRNHPP Brain-Hippocampus, BRNHPT Brain-Hypothalamus, BRNNCC Brain Nucleus accumbens(basal ganglia), BRNPTM Brain-Putamen (basal ganglia), BRNSPC Brain-Spinal cord(cervical c-1), BRNSNG Brain-Substantia nigra, BREAST Breast-Mammary Tissue, LCL Cells-EBV-transformed lymphocytes, FIBRBLS Cells-Transformed fibroblasts, CVXECT Cervix-Ectocervix, CVSEND Cervix-Endocervix, CLNSGM Colon-Sigmoid, CLNTRN Colon-Transverse, ESPGEJ Esophagus-Gastroesophageal Junction, ESPMCS Esophagus-Mucosa, ESPMSL Esophagus-Muscularis, FLLPNT Fallopian Tube, HRTAA Heart-Atrial Appendage, HRTLV Heart-Left Ventricle, KDNCTX Kidney-Cortex, LIVER Liver, LUNG Lung, SLVRYG Minor Salivary Gland, MSCLSK Muscle-Skeletal, NERVET Nerve-Tibial, OVARY Ovary, PNCREAS Pancreas, PTTARY Pituitary, PRSTTE Prostate, SKINNS Skin-Not Sun Exposed (Suprapubic), SKINS Skin-Sun Exposed (Lower leg), SNTTRM Small Intestine-Terminal Ileum, SPLEEN Spleen, STMACH Stomach, TESTIS Testis, THYROID Thyroid, UTERUS Uterus, VAGINA Vagina, WHLBLD Whole Blood"
tiss_abbrev <- strsplit(tiss_abbrev, ", ")[[1]]
tiss_abbrev <- do.call(rbind, lapply(strsplit(tiss_abbrev, " "), function(x) c(x[1], paste0(x[-1], collapse = " "))))
tiss_abbrev <- setNames(tiss_abbrev[,2], tiss_abbrev[,1])
names(rnaCounts_x_tiss) <- gsub(" ", "", tiss_abbrev[names(rnaCounts_x_tiss)])
GTEx_metadata <- fread("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# compare coverage / depth from rnaseq vs mmPCRseq

#first do mmPCRseq
sample_prop_tail_to_remove <- 0
min_n <- 0
xlocs <- 0:100/100
r <- 50
window_width <- 0.05
mmPCRseq_geom_expm_n <- parallel::mclapply(tissue_codes, function(tissue_code){
  mcprint(paste0(tissue_code, " "))
  
  dt <- fread(file = paste0(stan_data_dir, tissue_code, "_", 
                            sample_prop_tail_to_remove, ".csv"))
  samplesize_expw <- sapply(xlocs, function(xloc){
    ase_coefs <- abs(dt$refCount / dt$totalCount - 0.5) * 2
    in_window <- ase_coefs < (xloc + window_width/2) &
      ase_coefs > (xloc - window_width/2) 
      
    
    if(!any(in_window)){
      return(NA)
    } else {
      w <- dexp(abs(xloc - ase_coefs[in_window]), rate = r)
      return(sum(w * log10(dt$totalCount[in_window])) / sum(w))
    }
    
  })
  10^samplesize_expw
})

par(mar = c(5,6,4,2))
plot(xlocs, mmPCRseq_geom_expm_n[[1]], type = "l", ylim = c(0, max(unlist(mmPCRseq_geom_expm_n), na.rm = T)), col = adjustcolor(1, 0.5), 
     xlab = "sample allelic balance coefficient", ylab = "tissue-specific sample size\n(laplace smoothed weighted geometric mean)", main = "mmPCRseq")
for(i in 2:length(tissue_codes)){lines(xlocs, mmPCRseq_geom_expm_n[[i]], col = adjustcolor(1, 0.5))}
text(x = par("usr")[1], y = par("usr")[4] - strheight("Ay") * 2, labels = "each line is a different GTEx tissue\nred line averages across tissues", pos = 4)
lines(xlocs, 10^apply(log10(do.call(rbind, mmPCRseq_geom_expm_n)), 2, mean, na.rm=T), col = 2, lwd = 4)

#now do rnaseq

rnaseq_geom_expm_n <- parallel::mclapply(gsub(x = names(tissue_codes), pattern = " ", replacement = ""), function(tissue_name){
  dt <- rnaCounts_x_tiss[[tissue_name]]
  
  mcprint(paste0(tissue_name, ": ", nrow(dt), "; "))
  samplesize_expw <- sapply(xlocs, function(xloc){
    ase_coefs <- abs(dt$REF_COUNT / dt$TOTAL_COUNT - 0.5) * 2
    in_window <- ase_coefs < (xloc + window_width/2) &
      ase_coefs > (xloc - window_width/2) 
    
    
    if(!any(in_window)){
      return(NA)
    } else {
      w <- dexp(abs(xloc - ase_coefs[in_window]), rate = r)
      return(sum(w * log10(dt$TOTAL_COUNT[in_window])) / sum(w))
    }
    
  })
  exp_n <- 10^samplesize_expw
  # if(any(exp_n > 1E3)){mcprint(tissue_name)}
  return(exp_n)
})

plot(xlocs, rnaseq_geom_expm_n[[1]], type = "l", 
     ylim = c(0, quantile(unlist(rnaseq_geom_expm_n), 0.999, na.rm = T)), col = adjustcolor(1, 0.5), 
     xlab = "sample allelic balance coefficient", ylab = "tissue-specific sample size\n(laplace smoothed weighted geometric mean)", main = "rnaseq")
for(i in 2:length(tissue_codes)){lines(xlocs, rnaseq_geom_expm_n[[i]], col = adjustcolor(1, 0.5))}
text(x = par("usr")[1], y = par("usr")[4] - strheight("Ay") * 2, labels = "each line is a different GTEx tissue\nred line averages across tissues", pos = 4)
lines(xlocs, 10^apply(log10(do.call(rbind, rnaseq_geom_expm_n)), 2, mean, na.rm=T), col = 2, lwd = 4)


#now compare the two directly, site by site + x-reference the axes
window_width <- 0.1
use_rnaseq_tot <- T
merge_rna_mmpcr_geom_expm_n <- lapply(setNames(names(tissue_codes), names(tissue_codes)), function(tissue_name){
  
  dt_rnaseq <- rnaCounts_x_tiss[[gsub(x = tissue_name, pattern = " ", replacement = "")]]
  dt_mmpcrseq <- fread(file = paste0(stan_data_dir, tissue_codes[tissue_name], "_", 
                                     sample_prop_tail_to_remove, ".csv"))
  setDT(dt_rnaseq)
  setDT(dt_mmpcrseq)
  
  if(nrow(dt_rnaseq) < 1 || nrow(dt_mmpcrseq) < 1){
    return(NA)
  }
  
  #merge and rename
  dt_merge <- dt_rnaseq[dt_mmpcrseq, .(
    individual = i.individual,
    gene = i.gene,
    tissue = i.tissue,
    variantID = i.variantID,
    refCount_rnaseq = x.REF_COUNT,
    altCount_rnaseq = x.ALT_COUNT,
    totalCount_rnaseq = x.TOTAL_COUNT,
    refCount_mmpcrseq = i.refCount,
    altCount_mmpcrseq = i.altCount,
    totalCount_mmpcrseq = i.totalCount
  ),
  on = .(VARIANT_ID = variantID, individual = individual)
  ]
  dt_merge <- dt_merge[!apply(apply(dt_merge[,c("refCount_mmpcrseq", 
                                                "altCount_rnaseq", 
                                                "altCount_mmpcrseq", 
                                                "refCount_rnaseq")], 1, is.na), 2, any),]
  
  mcprint(paste0(tissue_name, ": ", nrow(dt), "; "))
  samplesize_expw <- sapply(xlocs, function(xloc){
    if(use_rnaseq_tot){
      ase_coefs <- abs(dt_merge$refCount_mmpcrseq / dt_merge$totalCount_mmpcrseq - 0.5) * 2
    } else {
      ase_coefs <- abs(dt_merge$refCount_rnaseq / dt_merge$totalCount_rnaseq - 0.5) * 2
    }
    in_window <- ase_coefs < (xloc + window_width/2) &
      ase_coefs > (xloc - window_width/2) 
    if(!any(in_window)){
      return(NA)
    } else {
      w <- dexp(abs(xloc - ase_coefs[in_window]), rate = r)
      if(use_rnaseq_tot){
        return(sum(w * log10(dt_merge$totalCount_rnaseq[in_window])) / sum(w))
      } else {
        return(sum(w * log10(dt_merge$totalCount_mmpcrseq[in_window])) / sum(w))
      }
    }
    
  })
  exp_n <- 10^samplesize_expw
  # if(any(exp_n > 1E3)){mcprint(tissue_name)}
  return(exp_n)
})
merge_rna_mmpcr_geom_expm_n <- merge_rna_mmpcr_geom_expm_n[!sapply(sapply(merge_rna_mmpcr_geom_expm_n, is.na), all)]

plot(xlocs, merge_rna_mmpcr_geom_expm_n[[1]], type = "l", 
     ylim = c(0, quantile(unlist(merge_rna_mmpcr_geom_expm_n), 0.99, na.rm = T)), col = adjustcolor(1, 0.5), 
     xlab = "sample allelic balance coefficient", ylab = "tissue-specific sample size\n(laplace smoothed weighted geometric mean)", 
     main = ifelse(use_rnaseq_tot, "rnaseq total count, mmPCRseq ASE coef", "mmPCRseq total count, rnaseq ASE coef"))
for(i in 2:length(merge_rna_mmpcr_geom_expm_n)){lines(xlocs, merge_rna_mmpcr_geom_expm_n[[i]], col = adjustcolor(1, 0.5))}
text(x = par("usr")[1], y = par("usr")[4] - strheight("Ay") * 2, labels = "each line is a different GTEx tissue\nred line averages across tissues", pos = 4)
lines(xlocs, 10^apply(log10(do.call(rbind, merge_rna_mmpcr_geom_expm_n)), 2, mean, na.rm=T), col = 2, lwd = 4)

#now compare the two directly, but use an indicator of > or < the median to the observed val
window_width <- 0.1
use_rnaseq_tot <- T
merge_rna_mmpcr_geom_expm_n <- lapply(setNames(names(tissue_codes), names(tissue_codes)), function(tissue_name){
  
  mcprint(paste0(tissue_name, ": "))

  dt_rnaseq <- rnaCounts_x_tiss[[gsub(x = tissue_name, pattern = " ", replacement = "")]]
  dt_mmpcrseq <- fread(file = paste0(stan_data_dir, tissue_codes[tissue_name], "_", 
                                     sample_prop_tail_to_remove, ".csv"))
  setDT(dt_rnaseq)
  setDT(dt_mmpcrseq)
  
  if(nrow(dt_rnaseq) < 1 || nrow(dt_mmpcrseq) < 1){
    return(NA)
  }
  
  
  #merge and rename
  dt_merge <- dt_rnaseq[dt_mmpcrseq, .(
    individual = i.individual,
    gene = i.gene,
    tissue = i.tissue,
    variantID = i.variantID,
    refCount_rnaseq = x.REF_COUNT,
    altCount_rnaseq = x.ALT_COUNT,
    totalCount_rnaseq = x.TOTAL_COUNT,
    refCount_mmpcrseq = i.refCount,
    altCount_mmpcrseq = i.altCount,
    totalCount_mmpcrseq = i.totalCount
  ),
  on = .(VARIANT_ID = variantID, individual = individual)
  ]
  
  dt_merge <- dt_merge[!apply(apply(dt_merge[,c("refCount_mmpcrseq", 
                                                "altCount_rnaseq", 
                                                "altCount_mmpcrseq", 
                                                "refCount_rnaseq")], 1, is.na), 2, any),]
  
  
  if(nrow(dt_merge) == 0){
    return(NA)
  }
  
  total_counts_by_indiv <- lapply(split(dt_merge[,c("totalCount_rnaseq", "totalCount_mmpcrseq")], dt_merge$individual), function(x) apply(x, 2, sum))
  total_counts_by_indiv <- data.frame(individual = names(total_counts_by_indiv), do.call(rbind, total_counts_by_indiv), row.names = names(total_counts_by_indiv))
  dt_merge$totalCount_rnaseq_prop <- dt_merge$totalCount_rnaseq / total_counts_by_indiv[dt_merge$individual,]$totalCount_rnaseq
  dt_merge$totalCount_mmpcrseq_prop <- dt_merge$totalCount_mmpcrseq / total_counts_by_indiv[dt_merge$individual,]$totalCount_mmpcrseq
  
  prop_counts_by_site <- lapply(split(dt_merge, dt_merge$variantID), function(x) cbind(x, q_rnaseq_prop = rank(x$totalCount_rnaseq_prop) / nrow(x), q_mmpcrseq_prop = rank(x$totalCount_mmpcrseq_prop) / nrow(x)))
  dt_merge <- do.call(rbind, prop_counts_by_site)
  
  
  average_q <- sapply(xlocs, function(xloc){
    if(use_rnaseq_tot){
      ase_coefs <- abs(dt_merge$refCount_mmpcrseq / dt_merge$totalCount_mmpcrseq - 0.5) * 2
    } else {
      ase_coefs <- abs(dt_merge$refCount_rnaseq / dt_merge$totalCount_rnaseq - 0.5) * 2
    }
    in_window <- ase_coefs < (xloc + window_width/2) &
      ase_coefs > (xloc - window_width/2) 
    if(!any(in_window)){
      return(NA)
    } else {
      w <- dexp(abs(xloc - ase_coefs[in_window]), rate = r)
      if(use_rnaseq_tot){
        return(sum(w * dt_merge$q_rnaseq_prop[in_window]) / sum(w))
      } else {
        return(sum(w * dt_merge$q_mmpcrseq_prop[in_window]) / sum(w))
      }
    }
    
  })
  return(average_q)
})
merge_rna_mmpcr_geom_expm_n <- merge_rna_mmpcr_geom_expm_n[!sapply(sapply(merge_rna_mmpcr_geom_expm_n, is.na), all)]

plot(xlocs, merge_rna_mmpcr_geom_expm_n[[1]], type = "l", 
     ylim = c(0, 1), col = adjustcolor(1, 0.5), 
     xlab = "sample allelic balance coefficient", ylab = "tissue-specific quantile\n(laplace smoothed weighted geometric mean)", 
     main = ifelse(use_rnaseq_tot, "rnaseq total count, mmPCRseq ASE coef", "mmPCRseq total count, rnaseq ASE coef"))
for(i in 2:length(merge_rna_mmpcr_geom_expm_n)){lines(xlocs, merge_rna_mmpcr_geom_expm_n[[i]], col = adjustcolor(1, 0.5))}
text(x = par("usr")[1], y = par("usr")[3] + strheight("Ay") * 2, labels = "each line is a different GTEx tissue\nred line averages across tissues", pos = 4)
lines(xlocs, 10^apply(log10(do.call(rbind, merge_rna_mmpcr_geom_expm_n)), 2, mean, na.rm=T), col = 2, lwd = 4)



#### Het-site Consistency, All Tissues ####

badhets <- ""
sample_prop_tail_to_remove <- c("0", "0.1")[1]
min_n <- 0
multitiss_by_gene_hetsite <- parallel::mclapply(tissue_codes, function(tissue_code){
  mcprint(tissue_code)
  
  dt <- fread(file = paste0(stan_data_dir, tissue_code, "_", 
                            sample_prop_tail_to_remove, ".csv"))
  dt[abs(0.5-dt$refCount / dt$totalCount) == 0.5,]})
  dt_removed <- fread(file = paste0(stan_data_dir, tissue_code, "_", sample_prop_tail_to_remove, "_removed.csv"))
  dt_dupes <- rbind(dt[dt$gene_indiv %in% dt_removed$gene_indiv,], dt_removed)
  dt_dupes <- dt_dupes[dt_dupes$totalCount > min_n,]
  dt_dupes <- dt_dupes[dt_dupes$gene_indiv %in% names(which(table(dt_dupes$gene_indiv) >= 2)),]
  
  
  dt_dupes$count_hdev <- abs(dt_dupes$refCount - dt_dupes$totalCount / 2)
  dt_dupes$prop_hdev <- dt_dupes$count_hdev / dt_dupes$totalCount * 2
  dt_dupes$pmean_prop_hdev <- abs((dt_dupes$refCount+1) / (dt_dupes$totalCount+2) - 0.5) / 0.5
  dt_dupes$q95_prop_hdev <- qbeta_dev(p = 0.95, shape1 = dt_dupes$refCount+1, shape2 = dt_dupes$altCount+1) / 0.5
  
  dt_dupesplit <- split(dt_dupes, dt_dupes$gene_indiv)
  dupemat <- do.call(rbind, lapply(dt_dupesplit, function(x){
    y <- as.data.frame(x[,c("gene", "individual", "totalCount", "prop_hdev", "q95_prop_hdev", "pmean_prop_hdev", "variantID")])
    combins <- combn(1:nrow(y), m = 2)
    out <- data.frame(individual = y$individual[1], 
                      gene = y$gene[1],
                      s1 = y[combins[1,],"variantID"],
                      s2 = y[combins[2,],"variantID"],
                      p1 = y[combins[1,],"prop_hdev"],
                      p2 = y[combins[2,],"prop_hdev"],
                      pmp1 = y[combins[1,],"pmean_prop_hdev"],
                      pmp2 = y[combins[2,],"pmean_prop_hdev"],
                      pq95p1 = y[combins[1,],"q95_prop_hdev"],
                      pq95p2 = y[combins[2,],"q95_prop_hdev"],
                      n1 = y[combins[1,],"totalCount"],
                      n2 = y[combins[2,],"totalCount"])
    return(out)
    
  }))
  
  #remove "bad" heterozygous sites
  dupemat <- dupemat[!(dupemat$s1 %in% badhets | dupemat$s2 %in% badhets),]
  
  dupemat_swap <- dupemat
  colnames(dupemat_swap) <- gsub(1, 3, colnames(dupemat_swap))
  colnames(dupemat_swap) <- gsub(2, 1, colnames(dupemat_swap))
  colnames(dupemat_swap) <- gsub(3, 2, colnames(dupemat_swap))
  dupemat_swap <- dupemat_swap[colnames(dupemat)]
  dupemat_comb <- rbind(dupemat, dupemat_swap)
  
  dupemat_g <- split(dupemat, dupemat$gene)
  dupemat_h <- split(dupemat_comb, dupemat_comb$s1)
  
  by_gene_ccc <- sapply(dupemat_g, function(x) ccc(logit(x[,c("pq95p1", "pq95p2")])))
  by_hetsite_ccc <- sapply(dupemat_h, function(x) ccc(logit(x[,c("pq95p1", "pq95p2")])))
  
  by_gene_MAD <- sapply(dupemat_g, function(x) MAD(logit(x[,c("pq95p1", "pq95p2")])))
  by_hetsite_MAD <- sapply(dupemat_h, function(x) MAD(logit(x[,c("pq95p1", "pq95p2")])))
  
  by_gene_MAD_raw <- sapply(dupemat_g, function(x) MAD((x[,c("pq95p1", "pq95p2")])))
  by_hetsite_MAD_raw <- sapply(dupemat_h, function(x) MAD((x[,c("pq95p1", "pq95p2")])))
  
  return(list(by_gene_ccc = by_gene_ccc, by_hetsite_ccc = by_hetsite_ccc,
              by_gene_MAD = by_gene_MAD, by_hetsite_MAD = by_hetsite_MAD,
              by_gene_MAD_raw = by_gene_MAD_raw, by_hetsite_MAD_raw = by_hetsite_MAD_raw))

}, mc.cores = 16)
inv_tissue_codes <- setNames(names(tissue_codes), tissue_codes)

#now aggregate to genes or hetsites for ccc
aggregate_by_gene_ccc <- lapply(inv_tissue_codes, function(tiss){
  tissmat <- multitiss_by_gene_hetsite[[tiss]][["by_gene_ccc"]]
  tissmat <- data.frame(gene = names(tissmat), 
             ccc = as.numeric(tissmat), 
             tissue = tissue_codes[tiss])
  tissmat
  
})
aggregate_by_gene_ccc <- do.call(rbind, aggregate_by_gene_ccc)
rownames(aggregate_by_gene_ccc) <- NULL

aggregate_by_hetsite_ccc <- lapply(inv_tissue_codes, function(tiss){
  tissmat <- multitiss_by_gene_hetsite[[tiss]][["by_hetsite_ccc"]]
  tissmat <- data.frame(gene = names(tissmat), 
                        ccc = as.numeric(tissmat), 
                        tissue = tissue_codes[tiss])
  tissmat
  
})
aggregate_by_hetsite_ccc <- do.call(rbind, aggregate_by_hetsite_ccc)
rownames(aggregate_by_hetsite_ccc) <- NULL

#do some plotting
layout(mat = rbind(c(1,1,2,2,5,5), c(1,1,2,2,6,6), c(3,3,4,4,7,7)))
hist(aggregate_by_gene_ccc$ccc, breaks = seq(-1,1,length.out = 50), 
     xlab = "within-gene concordance correlation coefficient", main = "all gene x tissue pairs")

hist(aggregate_by_hetsite_ccc$ccc, breaks = seq(-1,1,length.out = 50), 
     xlab = "within-het-site concordance correlation coefficient", main = "all het-site x tissue pairs")

#look at means within genes, hetsites, and tissues
gsplit <- split(aggregate_by_gene_ccc$ccc, aggregate_by_gene_ccc$gene)
hsplit <- split(aggregate_by_hetsite_ccc$ccc, aggregate_by_hetsite_ccc$gene)
tgsplit <- split(aggregate_by_gene_ccc$ccc, aggregate_by_gene_ccc$tissue)
thsplit <- split(aggregate_by_hetsite_ccc$ccc, aggregate_by_hetsite_ccc$tissue)

gmeans <- sapply(gsplit, mean, na.rm = T)
hmeans <- sapply(hsplit, mean, na.rm = T)
tgmeans <- sapply(tgsplit, mean, na.rm = T)
thmeans <- sapply(thsplit, mean, na.rm = T)

hist(gmeans, breaks = seq(-1,1,length.out = 50), 
     xlab = "within-gene concordance correlation coefficient", main = "average within genes")
hist(tgmeans, breaks = seq(-1,1,length.out = 50), 
     xlab = "within-gene concordance correlation coefficient", main = "average within tissues")
hist(hmeans, breaks = seq(-1,1,length.out = 50), 
     xlab = "within-het-site concordance correlation coefficient", main = "average within het-sites")
hist(thmeans, breaks = seq(-1,1,length.out = 50), 
     xlab = "within-het-site concordance correlation coefficient", main = "average within tissues")
plot.new()

#compare hmeans to earlier het-site annotation
hetsite_comparison <- data.frame(prop_3end = setNames(locus_info$proportion_to_3_end, locus_info$id)[names(hmeans)], 
                                 dist_3end = setNames(locus_info$nloci_to_3_end, locus_info$id)[names(hmeans)], 
                                 mean_ccc = hmeans,
                                 prop_trans = setNames(locus_info$prop_transcripts, locus_info$id)[names(hmeans)])
hetsite_comparison <- hetsite_comparison[-which(is.na(hetsite_comparison), arr.ind = T)[,1],]

plot(hetsite_comparison$prop_3end, hetsite_comparison$mean_ccc,
     xlab = "relative distance to 3' end", ylab = "mean concordance correlation within het site")
cor(hetsite_comparison, use = "comp")
xlocs <- 0:100/100
r <- 100
ccc_expw <- sapply(xlocs, function(xloc){
  w <- dexp(abs(xloc - hetsite_comparison$prop_3end), rate = r)
  sum(w * hetsite_comparison$mean_ccc) / sum(w)
})
lines(xlocs, ccc_expw, col = 2, lwd = 2)

plot(log10(hetsite_comparison$dist_3end), hetsite_comparison$mean_ccc, 
     xlab = latex2exp::TeX("log$_{10}$(absolute distance to 3' end in bp)"), 
     ylab = "mean concordance correlation within het site", axes=F)
axis(2); axis(1, at = 0:6, labels = 10^(0:6)); box(which = "plot")
cor(hetsite_comparison, use = "comp")
xlocs <- 0:60/10
r <- 20
ccc_expw <- sapply(xlocs, function(xloc){
  w <- dexp(abs(xloc - log10(hetsite_comparison$dist_3end)), rate = r)
  sum(w * hetsite_comparison$mean_ccc, na.rm = T) / sum(w, na.rm = T)
})
lines(xlocs, ccc_expw, col = 2, lwd = 2)

plot(hetsite_comparison$prop_trans, hetsite_comparison$mean_ccc,
     xlab = "proportion transcripts", ylab = "mean concordance correlation within het site")
xlocs <- 0:100/100
r <- 100
ccc_expw <- sapply(xlocs, function(xloc){
  w <- dexp(abs(xloc - hetsite_comparison$prop_trans), rate = r)
  sum(w * hetsite_comparison$mean_ccc) / sum(w)
})
lines(xlocs, ccc_expw, col = 2, lwd = 2)


# badhets <- names(hmeans)[hmeans < 0]
# badhets <- badhets[!is.na(badhets)]

#now aggregate to genes or hetsites for MAD
aggregate_by_gene_MAD <- lapply(inv_tissue_codes, function(tiss){
  tissmat <- multitiss_by_gene_hetsite[[tiss]][["by_gene_MAD"]]
  tissmat <- data.frame(gene = names(tissmat), 
                        MAD = as.numeric(tissmat), 
                        tissue = tissue_codes[tiss])
  tissmat
  
})
aggregate_by_gene_MAD <- do.call(rbind, aggregate_by_gene_MAD)
rownames(aggregate_by_gene_MAD) <- NULL

aggregate_by_hetsite_MAD <- lapply(inv_tissue_codes, function(tiss){
  tissmat <- multitiss_by_gene_hetsite[[tiss]][["by_hetsite_MAD"]]
  tissmat <- data.frame(gene = names(tissmat), 
                        MAD = as.numeric(tissmat), 
                        tissue = tissue_codes[tiss])
  tissmat
  
})
aggregate_by_hetsite_MAD <- do.call(rbind, aggregate_by_hetsite_MAD)
rownames(aggregate_by_hetsite_MAD) <- NULL

#do some plotting
xlims <- c(0,7)
hist(aggregate_by_gene_MAD$MAD, breaks = seq(0,100,length.out=500), xlim = xlims,
     xlab = "within-gene MAD", main = "all gene x tissue pairs")
text(labels = ifelse(any(aggregate_by_gene_MAD$MAD > xlims[2]), paste0(round(mean(aggregate_by_gene_MAD$MAD > xlims[2]) * 100, 3), "% out-of-bounds"), ""),
     x = par("usr")[2], y = par("usr")[4], pos = 2, col = 2, cex = 0.75, xpd = NA)

hist(aggregate_by_hetsite_MAD$MAD, breaks = seq(0,100,length.out=500), xlim = xlims,
     xlab = "within-het-site MAD", main = "all het-site x tissue pairs")
text(labels = ifelse(any(aggregate_by_hetsite_MAD$MAD > xlims[2]), paste0(round(mean(aggregate_by_hetsite_MAD$MAD > xlims[2]) * 100, 3), "% out-of-bounds"), ""),
     x = par("usr")[2], y = par("usr")[4], pos = 2, col = 2, cex = 0.75, xpd = NA)

gsplit <- split(aggregate_by_gene_MAD$MAD, aggregate_by_gene_MAD$gene)
hsplit <- split(aggregate_by_hetsite_MAD$MAD, aggregate_by_hetsite_MAD$gene)
tgsplit <- split(aggregate_by_gene_MAD$MAD, aggregate_by_gene_MAD$tissue)
thsplit <- split(aggregate_by_hetsite_MAD$MAD, aggregate_by_hetsite_MAD$tissue)

gmeans <- sapply(gsplit, mean, na.rm = T)
hmeans <- sapply(hsplit, mean, na.rm = T)
tgmeans <- sapply(tgsplit, mean, na.rm = T)
thmeans <- sapply(thsplit, mean, na.rm = T)

hist(gmeans, breaks = seq(0,100,length.out=500), xlim = xlims,
     xlab = "within-gene MAD", main = "average within genes")
text(labels = ifelse(any(gmeans > xlims[2]), paste0(round(mean(gmeans > xlims[2]) * 100, 3), "% out-of-bounds"), ""),
     x = par("usr")[2], y = par("usr")[4], pos = 2, col = 2, cex = 0.75, xpd = NA)

hist(tgmeans, breaks = seq(0,100,length.out=500), xlim = xlims,
     xlab = "within-gene MAD", main = "average within tissues")
text(labels = ifelse(any(tgmeans > xlims[2]), paste0(round(mean(tgmeans > xlims[2]) * 100, 3), "% out-of-bounds"), ""),
     x = par("usr")[2], y = par("usr")[4], pos = 2, col = 2, cex = 0.75, xpd = NA)

hist(hmeans, breaks = seq(0,100,length.out=500), xlim = xlims,
     xlab = "within-het-site MAD", main = "average within het-sites")
text(labels = ifelse(any(hmeans > xlims[2]), paste0(round(mean(hmeans > xlims[2]) * 100, 3), "% out-of-bounds"), ""),
     x = par("usr")[2], y = par("usr")[4], pos = 2, col = 2, cex = 0.75, xpd = NA)

hist(thmeans, breaks = seq(0,100,length.out=500), xlim = xlims,
     xlab = "within-het-site MAD", main = "average within tissues")
text(labels = ifelse(any(thmeans > xlims[2]), paste0(round(mean(thmeans > xlims[2]) * 100, 3), "% out-of-bounds"), ""),
     x = par("usr")[2], y = par("usr")[4], pos = 2, col = 2, cex = 0.75, xpd = NA)
plot.new()

#compare hmeans to earlier het-site annotation
hetsite_comparison <- data.frame(prop_3end = setNames(locus_info$proportion_to_3_end, locus_info$id)[names(hmeans)], 
                                 dist_3end = setNames(locus_info$nloci_to_3_end, locus_info$id)[names(hmeans)], 
                                 mean_MAD = hmeans,
                                 prop_trans = setNames(locus_info$prop_transcripts, locus_info$id)[names(hmeans)])
hetsite_comparison <- hetsite_comparison[-which(is.na(hetsite_comparison), arr.ind = T)[,1],]

plot(hetsite_comparison$prop_3end, hetsite_comparison$mean_MAD,
     xlab = "relative distance to 3' end", ylab = "mean MAD within het site")
cor(hetsite_comparison, use = "comp")
xlocs <- 0:100/100
r <- 100
ccc_expw <- sapply(xlocs, function(xloc){
  w <- dexp(abs(xloc - hetsite_comparison$prop_3end), rate = r)
  sum(w * hetsite_comparison$mean_MAD) / sum(w)
})
lines(xlocs, ccc_expw, col = 2, lwd = 2)

plot(log10(hetsite_comparison$dist_3end), hetsite_comparison$mean_MAD, 
     xlab = latex2exp::TeX("log$_{10}$(absolute distance to 3' end in bp)"), 
     ylab = "mean MAD within het site", axes=F)
axis(2); axis(1, at = 0:6, labels = 10^(0:6)); box(which = "plot")
cor(hetsite_comparison, use = "comp")
xlocs <- 0:60/10
r <- 20
ccc_expw <- sapply(xlocs, function(xloc){
  w <- dexp(abs(xloc - log10(hetsite_comparison$dist_3end)), rate = r)
  sum(w * hetsite_comparison$mean_MAD, na.rm = T) / sum(w, na.rm = T)
})
lines(xlocs, ccc_expw, col = 2, lwd = 2)

plot(hetsite_comparison$prop_trans, hetsite_comparison$mean_MAD,
     xlab = "proportion transcripts", ylab = "mean MAD within het site")
xlocs <- 0:100/100
r <- 100
ccc_expw <- sapply(xlocs, function(xloc){
  w <- dexp(abs(xloc - hetsite_comparison$prop_trans), rate = r)
  sum(w * hetsite_comparison$mean_MAD) / sum(w)
})
lines(xlocs, ccc_expw, col = 2, lwd = 2)

#what about MAD_raw?
aggregate_by_gene_MAD_raw <- lapply(inv_tissue_codes, function(tiss){
  tissmat <- multitiss_by_gene_hetsite[[tiss]][["by_gene_MAD_raw"]]
  tissmat <- data.frame(gene = names(tissmat), 
                        MAD_raw = as.numeric(tissmat), 
                        tissue = tissue_codes[tiss])
  tissmat
  
})
aggregate_by_gene_MAD_raw <- do.call(rbind, aggregate_by_gene_MAD_raw)
rownames(aggregate_by_gene_MAD_raw) <- NULL

aggregate_by_hetsite_MAD_raw <- lapply(inv_tissue_codes, function(tiss){
  tissmat <- multitiss_by_gene_hetsite[[tiss]][["by_hetsite_MAD_raw"]]
  tissmat <- data.frame(gene = names(tissmat), 
                        MAD_raw = as.numeric(tissmat), 
                        tissue = tissue_codes[tiss])
  tissmat
  
})
aggregate_by_hetsite_MAD_raw <- do.call(rbind, aggregate_by_hetsite_MAD_raw)
rownames(aggregate_by_hetsite_MAD_raw) <- NULL

#do some plotting
hist(aggregate_by_gene_MAD_raw$MAD_raw,  breaks = seq(-1,1,length.out = 50),  xlim = c(0,1),
     xlab = "within-gene MAD (raw scale)", main = "all gene x tissue pairs")
hist(aggregate_by_hetsite_MAD_raw$MAD_raw,  breaks = seq(-1,1,length.out = 50),  xlim = c(0,1),
     xlab = "within-het-site MAD (raw scale)", main = "all het-site x tissue pairs")

gsplit <- split(aggregate_by_gene_MAD_raw$MAD_raw, aggregate_by_gene_MAD_raw$gene)
hsplit <- split(aggregate_by_hetsite_MAD_raw$MAD_raw, aggregate_by_hetsite_MAD_raw$gene)
tgsplit <- split(aggregate_by_gene_MAD_raw$MAD_raw, aggregate_by_gene_MAD_raw$tissue)
thsplit <- split(aggregate_by_hetsite_MAD_raw$MAD_raw, aggregate_by_hetsite_MAD_raw$tissue)

gmeans <- sapply(gsplit, mean, na.rm = T)
hmeans <- sapply(hsplit, mean, na.rm = T)
tgmeans <- sapply(tgsplit, mean, na.rm = T)
thmeans <- sapply(thsplit, mean, na.rm = T)

hist(gmeans,  breaks = seq(-1,1,length.out = 50),  xlim = c(0,1),
     xlab = "within-gene MAD (raw scale)", main = "average within genes")
hist(tgmeans,  breaks = seq(-1,1,length.out = 50),  xlim = c(0,1),
     xlab = "within-gene MAD (raw scale)", main = "average within tissues")
hist(hmeans,  breaks = seq(-1,1,length.out = 50),  xlim = c(0,1),
     xlab = "within-het-site MAD (raw scale)", main = "average within het-sites")
hist(thmeans,  breaks = seq(-1,1,length.out = 50),  xlim = c(0,1),
     xlab = "within-het-site MAD (raw scale)", main = "average within tissues")
plot.new()

#subset het sites to just the self-consistent ones?



#### Load and Plot Data ####
tissueColors = fread("gtex_tissue_colors.txt", select = c(1:5))

targetsHG38 = fread("targeted-hg38.bed", col.names = c("chr","pos","pos2","gene","pool","id","ensg"))
targetsHG38$id = paste0("chr",targetsHG38$id)

counts <- fread("eGTEx-mmPCR-ASE-aggr.txt.gz")
names(counts)[1:5] = c("tissue","sample","run","chr","pos")

counts$variantID = paste0(counts$chr, "_",
                          counts$pos, "_",
                          counts$refAllele, "_", 
                          counts$altAllele, "_b38")

ovcids = fread("OVCids.txt", col.names = c("stage","grade","individual"))
idCorrection = fread("OVC-sampleIDs-correction.txt", col.names = c("individual","sample"))
ovcids = left_join(ovcids, idCorrection, by="individual")

tmpData = subset(counts, 
              altCount/totalCount <= (1-sample_prop_tail_to_remove) & 
              altCount/totalCount >= sample_prop_tail_to_remove &
              !is.na(tissue))

options(repr.plot.width = 6, repr.plot.height = 6)
ggplot() +
    geom_density(data = subset(tmpData, tissue != "OVC"),
                 aes(x = altCount/totalCount), fill = "black", alpha = 0.5) +
    geom_density(data = subset(tmpData, tissue == "OVC"),
                 aes(x = altCount/totalCount), fill = "red", alpha = 0.5) +
    theme_pubr(base_size = 16) +
    theme(legend.position = "none") +
    xlab("Allelic Ratio")

countsOVC = subset(counts, tissue == "OVC")
countsOVC = left_join(countsOVC, ovcids, by = "sample")

countsOVC$grade[which(is.na(countsOVC$grade) & countsOVC$stage == "Benign")] = "Benign,\nUnstaged"
countsOVC$grade[which(is.na(countsOVC$grade) & is.na(countsOVC$stage))] = "Benign,\nUnstaged"

countsOVC$grade[which(countsOVC$grade == "GB")] = "B"

options(repr.plot.width = 4, repr.plot.height = 6)
ggplot(subset(countsOVC,altCount/totalCount <= (1-sample_prop_tail_to_remove) & 
              altCount/totalCount >= sample_prop_tail_to_remove &
              !is.na(stage))) + 
    geom_violin(aes(x = stage,
                    y = abs(0.5-altCount/totalCount),
                    fill = stage), alpha = 0.5) +
    geom_boxplot(aes(x = stage,
                     y = abs(0.5-altCount/totalCount),
                     fill = stage), width = 0.15) +
    theme_pubr(base_size = 16) +
    xlab("") + ylab("Allelic Deviation") +
    scale_fill_manual(values = c("Benign" = "gray", "Malignant" = "red")) +
    theme(legend.position = "none")

options(repr.plot.width = 9, repr.plot.height = 6)
ggplot(subset(countsOVC, altCount/totalCount < 0.9 & 
              altCount/totalCount > 0.1 &
              !is.na(grade))) + 
    geom_violin(aes(x = fct_relevel(grade, c("Benign,\nUnstaged","X","B",
                                             "I","II","III")),
                    y = abs(0.5-altCount/totalCount),
                    fill = grade), alpha = 0.5) +
    geom_boxplot(aes(x = fct_relevel(grade, c("Benign,\nUnstaged","X","B",
                                             "I","II","III")),
                     y = abs(0.5-altCount/totalCount),
                     fill = grade), width = 0.15) +
    theme_pubr(base_size = 16) +
    xlab("") + ylab("Allelic Deviation")  +
    theme(legend.position = "none")

subset(countsOVC, altCount/totalCount <= (1-sample_prop_tail_to_remove) & 
              altCount/totalCount >= sample_prop_tail_to_remove &
              !is.na(grade)) %>% 
        group_by(grade) %>%
        summarize(count = length(unique(individual.y))) %>%
        ggplot() +
            geom_bar(aes(x = fct_relevel(grade, c("Benign,\nUnstaged","X","B","I","II","III")),
                         y = count,
                         fill = grade),
                     stat = "identity") +
            theme_pubr(base_size = 16) +
            theme(legend.position = "none") +
            xlab("") + ylab("Number of Tumor Samples")

#### Comparing with RNA-Seq ASE ####
countsRNA = fread("RNAseq_OVBM_ase_counts.tsv")

rnaCompare = countsRNA[,c(1,2,3,4,5,6,13,11,9)]
mmCompare = countsOVC[,c(4,5,7,8,9,10,17)]

mmRNACompare = left_join(mmCompare, rnaCompare, by = c("chr" = "CHROM",
                                                       "pos" = "POS",
                                                       "refAllele" = "REF",
                                                       "altAllele" = "ALT",
                                                       "individual.y" = "OTB_Ids"))

mmRNACompare = mmRNACompare %>%
                    filter(SNP_het == 1 & individual.y != "A00947102") %>%
                    mutate(mmRatio = altCount/(refCount + altCount),
                           rnaRatio = AD2/(AD1+AD2)) %>%
                    select(individual.y, Sample_Type, mmRatio, rnaRatio)

benign = ggplot(subset(mmRNACompare, Sample_Type != "Malignant")) +
    geom_point(aes(x = rnaRatio, y = mmRatio)) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    theme_pubr(base_size = 22.5, x.text.angle = 45) +
    facet_wrap(~individual.y) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    xlab("RNA-Seq Allelic Ratio") + ylab("mmPCR-Seq Allelic Ratio")

malignant = ggplot(subset(mmRNACompare, Sample_Type == "Malignant")) +
    geom_point(aes(x = rnaRatio, y = mmRatio),color = "red") +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    theme_pubr(base_size = 22.5, x.text.angle = 45) +
    facet_wrap(~individual.y) +
    theme(strip.text.x = element_blank())+
    xlab("RNA-Seq Allelic Ratio") + ylab("")

options(repr.plot.width = 30, repr.plot.height = 15)
benign + malignant

#### Statistical inference and comparison with normal ovary ####
countsOVC = subset(counts, tissue == "OVC")
countsOVC = left_join(countsOVC, ovcids, by = "sample")
countsOVC$grade[which(is.na(countsOVC$grade) & countsOVC$stage == "Benign")] = "Benign,\nUnstaged"
countsOVC$grade[which(is.na(countsOVC$grade) & is.na(countsOVC$stage))] = "Benign,\nUnstaged"
countsOVC$grade[which(countsOVC$grade == "GB")] = "B"

#generate table for 255 genes
ocvTestCount = counts %>% filter(variantID %in% targetsHG38$id & 
                                 tissue == "OVC" & 
                                 altCount/totalCount <= (1-sample_prop_tail_to_remove) & altCount/totalCount >= sample_prop_tail_to_remove)
ocvTestCount = left_join(ocvTestCount, targetsHG38[,c(4,6)], by=c("variantID" = "id"))
ocvTestCount = separate(ocvTestCount, gene, c("gene","num"),"_")
ocvTestCount = left_join(ocvTestCount, ovcids, by = "sample")
ocvTestCount$grade[which(is.na(ocvTestCount$grade) & ocvTestCount$stage == "Benign")] = "Benign,\nUnstaged"
ocvTestCount$grade[which(is.na(ocvTestCount$grade) & is.na(ocvTestCount$stage))] = "Benign,\nUnstaged"
ocvTestCount$grade[which(ocvTestCount$grade == "GB")] = "B"


fit = vglm(cbind(altCount, totalCount - altCount) ~ log(totalCount),
                      family = betabinomial(lmu = "identitylink", nsimEIM=100),
                      data = ocvTestCount)
    
mu_fitted = fitted(fit)
sigma_fitted = exp(fit@coefficients[2])
    
ocvTestCount$alpha0 = mu_fitted / sigma_fitted
ocvTestCount$beta0 = (1 - mu_fitted) / sigma_fitted
ocvTestCount$alpha1 = ocvTestCount$alpha0 + ocvTestCount$altCount
ocvTestCount$beta1 = ocvTestCount$beta0 + ocvTestCount$totalCount - ocvTestCount$altCount

ocvTestCount$pointEstimates = ocvTestCount$alpha1 / (ocvTestCount$alpha1 + ocvTestCount$beta1)

ocvTestCount$lowEstimate  = qbeta(.025, ocvTestCount$alpha1, ocvTestCount$beta1)
ocvTestCount$highEstimate = qbeta(.975, ocvTestCount$alpha1, ocvTestCount$beta1)

#test overlap with 0.5
ocvTestCount$posteriorProb = 2*rowMins(cbind(pbeta(.5, ocvTestCount$alpha1, ocvTestCount$beta1),
                                  1-pbeta(.5, ocvTestCount$alpha1, ocvTestCount$beta1)))

ocvTestCount = ocvTestCount %>%
                    arrange(posteriorProb) %>%
                    mutate(qvalue = cummean(posteriorProb))

geneTest = ocvTestCount %>%
            group_by(gene) %>%
            summarize(geneAlpha = mean(alpha1),
                      geneBeta = mean(beta1)) %>%
            mutate(pointEstimates = geneAlpha/(geneAlpha + geneBeta),
                   lowEstimate = qbeta(.025, geneAlpha, geneBeta),
                   highEstimate = qbeta(.975, geneAlpha, geneBeta))

#average posterior shape params over genes, 
geneTest$posteriorProb = 2*rowMins(cbind(pbeta(.5, geneTest$geneAlpha, geneTest$geneBeta),
                                  1-pbeta(.5, geneTest$geneAlpha, geneTest$geneBeta)))

qvaluesGTEx = dir_ls(path = "tissue-sumstats/", glob = "*") %>% map_dfr(fread)

geneTestGTEx = qvaluesGTEx %>%
                    group_by(tissue, tissue_color_hex, gene) %>%
                    summarize(geneAlpha = mean(alpha1),
                              geneBeta = mean(beta1)) %>%
                    mutate(pointEstimates = geneAlpha/(geneAlpha + geneBeta),
                           lowEstimate = qbeta(.025, geneAlpha, geneBeta),
                           highEstimate = qbeta(.975, geneAlpha, geneBeta))

geneTestGTEx$posteriorProb = 2*rowMins(cbind(pbeta(.5, geneTestGTEx$geneAlpha, geneTestGTEx$geneBeta),
                                  1-pbeta(.5, geneTestGTEx$geneAlpha, geneTestGTEx$geneBeta)))

geneTestGTExOvary = subset(geneTestGTEx, tissue == "Ovary")

geneTestGTExAll = qvaluesGTEx %>%
                    group_by(gene) %>%
                    summarize(geneAlpha = mean(alpha1),
                              geneBeta = mean(beta1)) %>%
                    mutate(pointEstimates = geneAlpha/(geneAlpha + geneBeta),
                           lowEstimate = qbeta(.025, geneAlpha, geneBeta),
                           highEstimate = qbeta(.975, geneAlpha, geneBeta))

geneTestGTExAll$posteriorProb = 2*rowMins(cbind(pbeta(.5, geneTestGTExAll$geneAlpha, geneTestGTExAll$geneBeta),
                                  1-pbeta(.5, geneTestGTExAll$geneAlpha, geneTestGTExAll$geneBeta)))

head(geneTestGTExAll)

ovcOvary = left_join(geneTest[,c(1,4)], geneTestGTExOvary[,c(3,6)], by = "gene")
names(ovcOvary) = c("gene","ovarian cancer","normal ovary")

ovcOvary_all = left_join(geneTest[,c(1:4)], geneTestGTExOvary[,c(3:6)], by = "gene")
names(ovcOvary_all) = c("gene",
                        "ovarian cancer alpha",
                        "ovarian cancer beta",
                        "ovarian cancer est",
                    "normal ovary alpha", 
                    "normal ovary beta", 
                    "normal ovary est")
ovcOvary_all <- ovcOvary_all[!apply(apply(ovcOvary_all, 1, is.na), 2, any),]

#### Stephen's dataframe additions ####

tissueColors = fread("gtex_tissue_colors.txt", select = c(1:5))
targetsHG38 = fread("targeted-hg38.bed", col.names = c("chr","pos","pos2","gene","pool","id","ensg"))
targetsHG38$id = paste0("chr",targetsHG38$id)

counts = fread("eGTEx-mmPCR-ASE-aggr.txt.gz")
names(counts)[1:5] = c("tissue","sample","run","chr","pos")

counts$variantID = paste0(counts$chr, "_",
                          counts$pos, "_",
                          counts$refAllele, "_", 
                          counts$altAllele, "_b38")

countsGtex = counts %>% filter(variantID %in% targetsHG38$id & tissue != "OVC")
countsGtex = left_join(countsGtex, tissueColors[,1:2], by=c("tissue" = "tissue_site_detail"))
countsGtex = left_join(countsGtex, targetsHG38[,c(4,6)], by=c("variantID" = "id"))
countsGtex <- countsGtex[countsGtex$tissue == "Ovary",]
countsGtex = separate(countsGtex, gene, c("gene","num"),"_")
rm(counts); gc()

#subset and combine data
data_GTEx <- subset(countsGtex, tissue_site_detail_abbr == "OVARY" & 
                      altCount/totalCount <= (1-sample_prop_tail_to_remove) & altCount/totalCount >= sample_prop_tail_to_remove)
data_GTEx <- data_GTEx[,c("tissue", "gene", "chr", "pos", "variantID", "refAllele", "altAllele", 
                          "refCount", "totalCount", "individual")]
data_GTEx$grade <- "GTEx"
data_OVC <- ocvTestCount[,c("tissue", "gene", "chr", "pos", "variantID", "refAllele", "altAllele", 
                            "refCount", "totalCount", "individual.x", "grade")]
colnames(data_OVC)[colnames(data_OVC) == "individual.x"] <- "individual"
genes <- intersect(data_GTEx$gene, data_OVC$gene)

data_OVC <- data_OVC[data_OVC$gene %in% genes,]
data_GTEx <- data_GTEx[data_GTEx$gene %in% genes,]
d <- rbind(data_OVC, data_GTEx)
shared_variants <- intersect(data_OVC$variantID, data_GTEx$variantID)
data_OVC_sv <- data_OVC[data_OVC$variantID %in% shared_variants,]
data_GTEx_sv <- data_GTEx[data_GTEx$variantID %in% shared_variants,]

all(table(paste0(d$variantID, "-", d$gene)) == table(d$variantID)) #confirm matchup of many-to-one variant-to-gene

#perform binomial test
# hist(d$refCount / d$totalCount)
# plot(d$refCount / d$totalCount, y = d$totalCount, pch = 19, col = adjustcolor((d$grade == "GTEx") + 1, 0.1), cex = 2)
mirror0 <- function(x) c(x, -x)
mirrorP <- function(x, p) c(x, 2 * p - x)
cutoff <- 0.5 - sample_prop_tail_to_remove
gtex_dens <- density(mirrorP(mirror0(abs(0.5 - data_GTEx$refCount / data_GTEx$totalCount)), cutoff), from = 0, to = cutoff, adjust = 1/5)
ovc_dens <- density(mirrorP(mirror0(abs(0.5 - data_OVC$refCount / data_OVC$totalCount)), cutoff), from = 0, to = cutoff, adjust = 1/5)

# gtex_dens <- density(mirrorP(mirror0(abs(0.5 - data_GTEx_sv$refCount / data_GTEx_sv$totalCount)), cutoff), from = 0, to = cutoff, adjust = 1/5)
# ovc_dens <- density(mirrorP(mirror0(abs(0.5 - data_OVC_sv$refCount / data_OVC_sv$totalCount)), cutoff), from = 0, to = cutoff, adjust = 1/5)


#plot distributions of distances
plot(NULL, NULL, xlim = range(gtex_dens$x), ylim = c(0, max(c(gtex_dens$y, ovc_dens$y))), 
     xlab = "Distance of Sample ASE proportion from 0.5", ylab = "Density", 
     main = "Comparing Sample Proportions (ASE) of GTEx and Ovarian Cancer Samples\nNaive Pooling Across All Samples (Genes, Individuals, & Heterozygous Sites)")
ovary_cols <- c(GTEx = "blue", OVC = "red", diff = "orange", sigdiff = "purple")
ks_ase_coefs <- ks.test(abs(0.5 - data_OVC$refCount / data_OVC$totalCount), abs(0.5 - data_GTEx$refCount / data_GTEx$totalCount))
polygon(x = c(gtex_dens$x, rev(gtex_dens$x)),
        y = c(gtex_dens$y, rep(0, length(gtex_dens$y))), col = adjustcolor(ovary_cols["GTEx"], 0.2))
polygon(x = c(ovc_dens$x, rev(ovc_dens$x)),
        y = c(ovc_dens$y, rep(0, length(ovc_dens$y))), col = adjustcolor(ovary_cols["OVC"], 0.2))
legend("topright", legend = c("Healthy GTEx Ovary", "Ovarian Cancer"), col = adjustcolor(ovary_cols, 0.2), pch = 15, pt.cex = 2, cex = 1.25)

#MC simulation hypothesis testing of this density curve
perform_HS_of_curves <- F
if(perform_HS_of_curves){
  n_GTEx <- nrow(data_GTEx)
  total_obs_GTEx <- sum(data_GTEx$totalCount)
  data_GTEx$prop_row <- data_GTEx$totalCount / total_obs_GTEx
  n_OVC <- nrow(data_OVC)
  total_obs_OVC <- sum(data_OVC$totalCount)
  data_OVC$prop_row <- data_OVC$totalCount / total_obs_OVC
  
  mc_data_GTEx <- data_GTEx
  mc_data_OVC <- data_OVC
  
  #simulate bernoulli draws and compute densities
  n_sim <- 5E3
  mcprint <- function(...){
    system(sprintf('echo "%s"', paste0(..., collapse="")))
  }
  
  GTEx_dens <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_GTEx_row_inds <- table(sample(1:nrow(data_GTEx), size = total_obs_GTEx, replace = T, prob = data_GTEx$prop_row))
    mc_data_GTEx$mc_totalCount <- 0
    mc_data_GTEx$mc_totalCount[as.integer(names(mc_GTEx_row_inds))] <- mc_GTEx_row_inds
    mc_data_GTEx$mc_refCount <- sapply(1:n_GTEx, function(i) rbinom(1, mc_data_GTEx$mc_totalCount[i], 
                                                                    prob = mc_data_GTEx$refCount[i] / mc_data_GTEx$totalCount[i]))
    mc_gtex_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_GTEx$mc_refCount / mc_data_GTEx$mc_totalCount)), cutoff), 
                            from = 0, to = cutoff, adjust = 1/5)$y
    mc_gtex_dens
  }, mc.cores = 12)
  GTEx_dens <- GTEx_dens[sapply(GTEx_dens, class) != "try-error"]
  
  OVC_dens <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_OVC_row_inds <- table(sample(1:nrow(data_OVC), size = total_obs_OVC, replace = T, prob = data_OVC$prop_row))
    mc_data_OVC$mc_totalCount <- 0
    mc_data_OVC$mc_totalCount[as.integer(names(mc_OVC_row_inds))] <- mc_OVC_row_inds
    mc_data_OVC$mc_refCount <- sapply(1:n_OVC, function(i) rbinom(1, mc_data_OVC$mc_totalCount[i], 
                                                                  prob = mc_data_OVC$refCount[i] / mc_data_OVC$totalCount[i]))
    mc_ovc_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_OVC$mc_refCount / mc_data_OVC$mc_totalCount)), cutoff), 
                           from = 0, to = cutoff, adjust = 1/5)$y
  }, mc.cores = 15)
  OVC_dens <- OVC_dens[sapply(OVC_dens, class) != "try-error"]
  
  #can also treat rows as observations?
  GTEx_dens_rowobs <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_GTEx_row_inds <- sample(1:nrow(data_GTEx), size = n_GTEx, replace = T, prob = data_GTEx$prop_row)
    mc_gtex_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_GTEx$refCount[mc_GTEx_row_inds] / mc_data_GTEx$totalCount[mc_GTEx_row_inds])), cutoff), 
                            from = 0, to = cutoff, adjust = 1/5)$y
    mc_gtex_dens
  }, mc.cores = 4)
  
  OVC_dens_rowobs <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_OVC_row_inds <- sample(1:nrow(data_OVC), size = n_OVC, replace = T, prob = data_OVC$prop_row)
    mc_ovc_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_OVC$refCount[mc_OVC_row_inds] / mc_data_OVC$totalCount[mc_OVC_row_inds])), cutoff), 
                           from = 0, to = cutoff, adjust = 1/5)$y
    mc_ovc_dens
  }, mc.cores = 4)
  
  #and sample then uniformly?
  GTEx_dens_rowobs_unif <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_GTEx_row_inds <- sample(1:nrow(data_GTEx), size = n_GTEx, replace = T)
    mc_gtex_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_GTEx$refCount[mc_GTEx_row_inds] / mc_data_GTEx$totalCount[mc_GTEx_row_inds])), cutoff), 
                            from = 0, to = cutoff, adjust = 1/5)$y
    mc_gtex_dens
  }, mc.cores = 4)
  
  OVC_dens_rowobs_unif <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_OVC_row_inds <- sample(1:nrow(data_OVC), size = n_OVC, replace = T)
    mc_ovc_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_OVC$refCount[mc_OVC_row_inds] / mc_data_OVC$totalCount[mc_OVC_row_inds])), cutoff), 
                           from = 0, to = cutoff,  adjust = 1/5)$y
    mc_ovc_dens
  }, mc.cores = 4)
  
  #or people as observations?
  indivs_data_GTEx <- split(data_GTEx, data_GTEx$individual)
  GTEx_dens_indiv <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_GTEx_inds <- sample(1:length(indivs_data_GTEx), 
                                     size = length(indivs_data_GTEx), 
                                     replace = T)
    mc_data_GTEx_indivs <- do.call(rbind, indivs_data_GTEx[mc_GTEx_inds])
    mc_gtex_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_GTEx_indivs$refCount / mc_data_GTEx_indivs$totalCount)), cutoff), 
                            from = 0, to = cutoff, adjust = 1/5)$y
    mc_gtex_dens
  }, mc.cores = 12)
  
  indivs_data_OVC <- split(data_OVC, data_OVC$individual)
  OVC_dens_indiv <- parallel::mclapply(1:n_sim, function(sim_i){
    mcprint(paste0(sim_i, " "))
    mc_OVC_inds <- sample(1:length(indivs_data_OVC), 
                           size = length(indivs_data_OVC), 
                           replace = T)
    mc_data_OVC_indivs <- do.call(rbind, indivs_data_OVC[mc_OVC_inds])
    mc_OVC_dens <- density(mirrorP(mirror0(abs(0.5 - mc_data_OVC_indivs$refCount / mc_data_OVC_indivs$totalCount)), cutoff), 
                            from = 0, to = cutoff, adjust = 1/5)$y
    mc_OVC_dens
  }, mc.cores = 12)
  
  #now generate the plots
  dists_from_0.5 <- seq(from = 0, to = cutoff, length.out = 512)
  
  GTEx_dens_95QI <- t(apply(do.call(rbind, GTEx_dens), 2, quantile, probs = c(0.05, 0.95)))
  OVC_dens_95QI <- t(apply(do.call(rbind, OVC_dens), 2, quantile, probs = c(0.05, 0.95)))
  minvalrep <- min(c(length(GTEx_dens), length(OVC_dens)))
  diff_dens_95QI <- t(apply(do.call(rbind, OVC_dens[1:minvalrep]) - 
                              do.call(rbind, GTEx_dens[1:minvalrep]), 
                            2, quantile, probs = c(0.05, 0.95)))
  
  GTEx_dens_95QI <- t(apply(do.call(rbind, GTEx_dens_indiv), 2, quantile, probs = c(0.05, 0.95)))
  OVC_dens_95QI <- t(apply(do.call(rbind, OVC_dens_indiv), 2, quantile, probs = c(0.05, 0.95)))
  minvalrep <- min(c(length(GTEx_dens_indiv), length(OVC_dens_indiv)))
  diff_dens_95QI <- t(apply(do.call(rbind, OVC_dens_indiv[1:minvalrep]) - 
                              do.call(rbind, GTEx_dens_indiv[1:minvalrep]), 
                            2, quantile, probs = c(0.05, 0.95)))
  
  overlaps0 <- apply(diff_dens_95QI < 0, 1, sum) == 1
  signif_diff <- !overlaps0
  rlesig <- rle(signif_diff)
  breaks <- cumsum(c(0, rlesig$lengths))
  breaks[1] <- 1
  signif_diff_segms <- do.call(rbind, lapply(seq_along(rlesig$values), function(tfi){
    if(rlesig$values[tfi]) return(c(breaks[tfi],breaks[tfi+1]))
  }))
  signif_diff_segms[,2] <- signif_diff_segms[,2] - 1
  # GTEx_dens_95QI <- t(apply(do.call(rbind, GTEx_dens_rowobs), 2, quantile, probs = c(0.05, 0.95)))
  # OVC_dens_95QI <- t(apply(do.call(rbind, OVC_dens_rowobs), 2, quantile, probs = c(0.05, 0.95)))
  # 
  # GTEx_dens_95QI <- t(apply(do.call(rbind, GTEx_dens_rowobs_unif), 2, quantile, probs = c(0.05, 0.95)))
  # OVC_dens_95QI <- t(apply(do.call(rbind, OVC_dens_rowobs_unif), 2, quantile, probs = c(0.05, 0.95)))
  
  ylims <- c(min(c(0,diff_dens_95QI)), max(c(GTEx_dens_95QI, OVC_dens_95QI)))
  ylims[1] <- ylims[1] - diff(ylims) / 10
  plot(NULL, NULL, xlim = range(dists_from_0.5), ylim = ylims, 
       xlab = "Distance of Sample ASE proportion from 0.5", ylab = "Density", 
       main = "Comparing Sample Proportions (ASE) of GTEx and Ovarian Cancer Samples\nNaive Pooling Across All Samples (Genes, Individuals, & Heterozygous Sites)")
  # ks_ase_coefs <- ks.test(abs(0.5 - data_OVC$refCount / data_OVC$totalCount), abs(0.5 - data_GTEx$refCount / data_GTEx$totalCount))
  for(i in 1:nrow(signif_diff_segms)){
    segments(x0 = dists_from_0.5[signif_diff_segms[i,1]], 
             x1 = dists_from_0.5[signif_diff_segms[i,2]],
             y0 = ylims[1], y1 = ylims[1], xpd = NA, col = ovary_cols["sigdiff"], lwd = 2)
    segxvals <- dists_from_0.5[signif_diff_segms[i,1]:signif_diff_segms[i,2]]
    polygon(x = c(segxvals, rev(segxvals)),
            y = c(rep(ylims[2], length(segxvals)), rep(ylims[1], length(segxvals))), 
            col = adjustcolor(ovary_cols["sigdiff"], 0.1), border = NA)
  }
  polygon(x = c(dists_from_0.5, rev(dists_from_0.5)),
          y = c(GTEx_dens_95QI[,2], rev(GTEx_dens_95QI[,1])), 
          col = "white")
  polygon(x = c(dists_from_0.5, rev(dists_from_0.5)),
          y = c(OVC_dens_95QI[,2], rev(OVC_dens_95QI[,1])), 
          col = "white")
  polygon(x = c(dists_from_0.5, rev(dists_from_0.5)),
          y = c(diff_dens_95QI[,2], rev(diff_dens_95QI[,1])), 
          col = "white")
  
  polygon(x = c(dists_from_0.5, rev(dists_from_0.5)),
          y = c(GTEx_dens_95QI[,2], rev(GTEx_dens_95QI[,1])), 
          col = adjustcolor(ovary_cols["GTEx"], 0.2))
  polygon(x = c(dists_from_0.5, rev(dists_from_0.5)),
          y = c(OVC_dens_95QI[,2], rev(OVC_dens_95QI[,1])), 
          col = adjustcolor(ovary_cols["OVC"], 0.2))
  polygon(x = c(dists_from_0.5, rev(dists_from_0.5)),
          y = c(diff_dens_95QI[,2], rev(diff_dens_95QI[,1])), 
          col = adjustcolor(ovary_cols["diff"], 0.2))
  abline(h = 0, col = adjustcolor(1, 0.5), lwd = 2, lty = 2)
  
  legend("topright", legend = c("Healthy GTEx Ovary", "Ovarian Cancer", "Difference", "95% Diff Excludes 0"), 
         col = c(adjustcolor(ovary_cols[1:3], 0.2), ovary_cols[4]), pch = c(15, 15, 15, NA), lwd = c(NA, NA, NA, 2),
         border = 1, pt.cex = 2, cex = 1.25)
  
  
}

#plot distributions of p-values
d$binom_p <- sapply(1:nrow(d), function(i) binom.test(d$refCount[i], d$totalCount[i], alternative = "two.sided", p = 0.5)$p.value)
# d$binom_p_adj <- p.adjust(d$binom_p, "BH")
d$binom_p_adj <- IHW::ihw(d$binom_p, factor(d$grade), alpha = 0.1)@df$adj_pvalue
pval_breaks <- 0:20/20
ks_pvals <- ks.test(d$binom_p_adj[d$grade == "GTEx"], d$binom_p_adj[d$grade != "GTEx"])
hist(d$binom_p_adj, breaks = pval_breaks, col = adjustcolor(ovary_cols["OVC"], 0.2),
     main = paste0("Two-sample KS-test p-value = ", format(ks_pvals$p.value, scientific = T, digits = 3)),
     xlab = "Binomial Test IHW-adjusted (~ grade) P-Values (testing deviation of ASE prop from 0.5)")
hist(d$binom_p_adj[d$grade == "GTEx"], breaks = pval_breaks, col = "white", add = T)
hist(d$binom_p_adj[d$grade == "GTEx"], breaks = pval_breaks, col = adjustcolor(ovary_cols["GTEx"], 0.2), add = T)
legend("topright", legend = c("Healthy GTEx Ovary", "Ovarian Cancer"), col = adjustcolor(ovary_cols, 0.2), pch = 15, pt.cex = 2, cex = 1.25)

gtex_p_dens <- density(mirror0(log10(d$binom_p_adj[d$grade == "GTEx"])), from = 0, to = min(log10(d$binom_p_adj[d$binom_p_adj != 0])))
# gtex_p_dens$y <- log10(gtex_p_dens$y)
ovc_p_dens <- density(mirror0(log10(d$binom_p_adj[d$grade != "GTEx"])), from = 0, to = min(log10(d$binom_p_adj[d$binom_p_adj != 0])))
# ovc_p_dens$y <- log10(ovc_p_dens$y)

plot(NULL, NULL, xlim = c(-50,0), ylim = c(0, max(c(gtex_p_dens$y, ovc_p_dens$y), na.rm = T)), 
     xlab = latex2exp::TeX("log$_{10}$(p-value) from Binomial Test"), ylab = "Density", 
     main = "Comparing log10(P-Value) Distributions\nfrom Binomial Tests of ASE Across GTEx and Ovarian Cancer Samples\nNaive Pooling Across All Samples (Genes, Individuals, & Heterozygous Sites)", cex.main = 0.9)
lines(gtex_p_dens, col = ovary_cols["GTEx"])
lines(ovc_p_dens, col = ovary_cols["OVC"])
legend("topleft", legend = c("Healthy GTEx Ovary", "Ovarian Cancer"), col = ovary_cols, lwd = c(2,2), cex = 1.25)

#match individuals according to ASE variant
dsh <- d[d$variantID %in% shared_variants,]
adj_alpha <- 0.05
dsh$bh_sig <- dsh$binom_p_adj <= adj_alpha
dsh_split <- split(dsh, dsh$variantID)

dsig <- do.call(rbind, lapply(dsh_split, function(x){
  out <- x[1,1:7]
  gtex_sig <- x$bh_sig[x$grade == "GTEx"]
  ovc_sig <- x$bh_sig[x$grade != "GTEx"]
  table(ovc_sig)[c("TRUE", "FALSE")]
  out <- cbind(out, data.frame(GTEx_ASE = sum(gtex_sig), 
                               GTEx_total = length(gtex_sig),
                               OVC_ASE = sum(ovc_sig), 
                               OVC_total = length(ovc_sig)))
  two_x_two <- cbind(GTEx = c(ASE = out$GTEx_ASE, not_ASE = out$GTEx_total - out$GTEx_ASE),
        OVC = c(ASE = out$OVC_ASE, not_ASE = out$OVC_total - out$OVC_ASE))
  chisq_res <- suppressWarnings(chisq.test(two_x_two, correct=T))
  out$chi2p <- chisq_res$p.value
  exactTest_res <- fisher.test(two_x_two)
  out$fisher_exact_p <- exactTest_res$p.value
  out$OCV_more_ASE <-  (out$OVC_ASE / out$OVC_total) > (out$GTEx_ASE / out$GTEx_total)
  out
}))

hist(dsig$chi2p, breaks = 100, xlab = latex2exp::TeX(paste0("Nominal P-value from $\\Chi^{2}$ Test of 2 x 2 Contingency Table")),
     main = "")
title(latex2exp::TeX("Testing Each Het Site Individually for Differences in Proportion"), line = 2)
title(latex2exp::TeX("\"Significant\" ($\\alpha$ = 0.05) Individuals w/ Detectable ASE (GTEx Sample vs. Ovarian Cancer)"), line = 1)

hist(apply(cbind(dsig$GTEx_total, dsig$OVC_total), 1, min), xlab = "Number of Individuals", 
     main = "Minimum # of Total Individuals Across \nGTEx & Cancer Samples for Each Het Site")

#### Integrating COSMIC annotation ####
cosmic <- read.csv("cancer_annotations/COSMIC_Census_allThu Nov  2 19_33_39 2023.csv")
cosmic <- cosmic[cosmic$Tier == 1,]
cosmic_OVC <- cosmic[grepl("ovar", cosmic$Tumour.Types.Somatic.) | grepl("ovar", cosmic$Tumour.Types.Germline.),]
cosmic_OVC$Gene.Symbol
sum(unique(data_OVC$gene) %in% cosmic_OVC$Gene.Symbol)


data_OVC_cosm_sub <- data_OVC[(data_OVC$gene) %in% cosmic_OVC$Gene.Symbol,]
hist(data_OVC_cosm_sub$refCount / data_OVC_cosm_sub$totalCount, breaks = 1:9/10)
hist(data_OVC$refCount / data_OVC$totalCount, breaks = 1:9/10)

hist(d$binom_p_adj[d$gene %in% cosmic_OVC$Gene.Symbol], breaks = 0:20/20, freq = F)
hist(d$binom_p_adj[!(d$gene %in% cosmic_OVC$Gene.Symbol)], breaks = 0:20/20, freq = F)

hist(d$binom_p_adj[d$gene %in% cosmic_OVC$Gene.Symbol & d$grade != "GTEx"], 
     breaks = 0:20/20, freq = F, col = adjustcolor("blue", 0.2))

#### making some graphs and follow-up analyses ####

#kde plots for ovarian cancer vs healty tissue
ovcOvary <- ovcOvary[!apply(apply(ovcOvary, 1, is.na), 2, any),]
plot(ovcOvary$`normal ovary`, ovcOvary$`ovarian cancer`)
ndens <- density(ovcOvary$`normal ovary`, from = 0, to = 1, n = 2^10)
ndens <- density(geneTestGTEx$pointEstimates, from = 0, to = 1, n = 2^10) #aggregates over tissues
cdens <- density(ovcOvary$`ovarian cancer`, from = 0, to = 1, n = 2^10)
cdens <- density(ocvTestCount$pointEstimates, from = 0, to = 1, n = 2^10) #not unique per gene
plot(NULL, xlim = c(0,1), ylim = range(c(ndens$y, cdens$y)), xlab = "allelic expression coefficient",
     ylab = "density")
polygon(c(ndens$x, rev(ndens$x)), c(ndens$y, rep(0, length(ndens$y))), col = adjustcolor(1, 0.5))
polygon(c(cdens$x, rev(cdens$x)), c(cdens$y, rep(0, length(cdens$y))), col = adjustcolor(2, 0.5))

#q-q plot for this comparison

#raw 'data'
plot(quantile(ovcOvary$`normal ovary`, 0:100/100),
     quantile(ovcOvary$`ovarian cancer`, 0:100/100), type = "l",
     xlab = "normal ovary", ylab = "ovarian cancer", xlim = c(0,1), ylim = c(0,1))
abline(0,1,lty = 2, col =adjustcolor(1,0.5), lwd = 3)

#kdes
locs <- ndens$x
ps <- 0:100/100
nqs <- cumsum(ndens$y) / 2^10 #normal quantiles for above pts
cqs <- cumsum(cdens$y) / 2^10 #cancer quantiles for above pts
nxs <- ndens$x[sapply(ps, function(qi) max(1, findInterval(qi, nqs)))]
cxs <- cdens$x[sapply(ps, function(qi) max(1, findInterval(qi, cqs)))]
# plot(nqs, cqs[sapply(nqs, function(qi) max(1, findInterval(qi, cqs)))]); abline(0,1,col=2)
plot(nxs, cxs, type = "l", xlab = "normal ovary", ylab = "ovarian cancer", main = "q-q plot",
     xlim = c(0,1), ylim = c(0,1))
abline(0,1,lty = 2, col =adjustcolor(1,0.5), lwd = 3)

#ftest
var.test(ovcOvary$`normal ovary`, ovcOvary$`ovarian cancer`, alternative = "two.sided")

#gene-wise comparison
n_samp <- 5E3
posterior_probs_cancer_extreme <- do.call(rbind, lapply(1:nrow(ovcOvary_all), function(i){
  x <- ovcOvary_all[i,]
  snormal <- abs(rbeta(n_samp, 
                   x$`normal ovary alpha`,
                   x$`normal ovary beta`) - 0.5)
  scancer <- abs(rbeta(n_samp, 
                   x$`ovarian cancer alpha`,
                   x$`ovarian cancer beta`) - 0.5)
  sdiff <- snormal - scancer
  return(data.frame(gene = x$gene, 
           prob_cancer_more_extreme = mean(sdiff < 0),
           expectation_of_difference = mean(sdiff)))
})
)
ovc <- cbind(posterior_probs_cancer_extreme, ovcOvary)
hist(main = "", posterior_probs_cancer_extreme$prob_cancer_more_extreme)
plot(abs(ovc$`ovarian cancer` - 0.5) - abs(ovc$`normal ovary` - 0.5),
     ovc$prob_cancer_more_extreme)


#### Plotting Results ####
options(repr.plot.height = 7.5, repr.plot.width = 7.5)
ggplot() +
    geom_point(data = ovcOvary,
               aes(y = `ovarian cancer`, x = `normal ovary`),size = 3) +
    theme_pubr(base_size = 16) +
    xlim(c(0.2, 0.8)) + ylim(c(0.2, 0.8)) +
    geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5)

ovcOvary %>% arrange(-abs(`ovarian cancer` - `normal ovary`)) %>% head()

p1 = as.data.frame(subset(geneTest, gene == "PRDM1"))
p2 = as.data.frame(subset(geneTestGTEx, gene == "PRDM1", select = c(names(p1),"tissue","tissue_color_hex")))

p1$tissue = "OVC"
p1$tissue_color_hex = "#000000"

combined = rbind(p1,p2)

options(repr.plot.width = 15, repr.plot.height = 7.5)
    ggplot(combined) +
        geom_hline(yintercept = 0.5) +
        geom_pointrange(aes(x = reorder(tissue, pointEstimates), y = pointEstimates,
                            ymin = lowEstimate, ymax = highEstimate, color = tissue_color_hex), size = 1.5) +
        theme_pubr(base_size = 16, x.text.angle = 45) +
        scale_color_identity() +
        theme(plot.margin = margin(l=75)) + xlab("") + ylab("Estimated PRDM1 Allelic Ratio")

sessionInfo()


