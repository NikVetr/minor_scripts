library(bigsnpr)
library(data.table)

# gwas <- fread("~/data/f.20016.0.0_res.EUR.sumstats.MACfilt.txt.gz")
gwas <- fread("~/data/GWAS_EA_excl23andMe.txt")
genotypes <- fread("~/data/nik_AncestryDNA.txt")
mean(genotypes$rsid %in% gwas$MarkerName)
mean(gwas$MarkerName %in% genotypes$rsid)
hist(gwas$Pval)

#ld scores
info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"))

#munge sumstats
sumstats <- gwas[gwas$MarkerName %in% info$rsid,]
sumstats$na <- 0.5
sumstats_munged <- sumstats[,c("CHR", "POS", "MarkerName", "A1", "A2", "na", "SE", "Pval", "Beta", "na", "na")]
names(sumstats_munged) <- c("chr", "pos", "rsid", "a1", "a0", "n_eff", "beta_se", "p", "beta", "INFO", "MAF")
sumstats_munged$n_eff <- 766345

#constrain to matches in my genotype file
sumstats_munged <- sumstats_munged[sumstats_munged$rsid %in% genotypes$rsid,]

# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file
fam.order <- NULL
# preprocess the bed file (only need to do once for each data set)
snp_readBed("~/data/post-qc/EUR.QC.bed")
# now attach the genotype object
obj.bigSNP <- snp_attach("~/data/post-qc/EUR.QC.rds")
obj.bigSNP$genotypes

#add in my own data
library(data.table)
genotypes <- fread("~/data/nik_AncestryDNA.txt")
aligned_genotypes <- merge(obj.bigSNP$map, genotypes, by.x = "marker.ID", by.y = "rsid", all.x = TRUE)

# Assuming 'allele1' and 'allele2' are columns from your `new_genotypes` data.frame
code_genotype <- function(ref_allele, alt_allele, allele1, allele2) {
  if (is.na(allele1) || is.na(allele2)) {
    return(NA_integer_)  # missing genotype
  }

  alleles <- c(allele1, allele2)
  if (all(alleles == ref_allele)) {
    return(as.integer(0))  # homozygous reference
  } else if (all(alleles == alt_allele)) {
    return(as.integer(2))  # homozygous alternative
  } else {
    return(as.integer(1))  # heterozygous
  }
}

# Apply the conversion for each SNP in aligned_genotypes
# You need to ensure that you're using the correct reference and alternative alleles from your map.
genotype_vector <- mapply(code_genotype,
                          ref_allele = aligned_genotypes$allele1.x,
                          alt_allele = aligned_genotypes$allele2.x,
                          allele1 = aligned_genotypes$allele1.y,
                          allele2 = aligned_genotypes$allele2.y)

obj.bigSNP$genotypes[1,] <- as.numeric(genotype_vector)
# new_genotype_matrix <- matrix(genotype_vector, nrow = 1, ncol = length(genotype_vector))
# all_genotypes <- obj.bigSNP$genotypes  # This is a big.matrix
# 
# # Add a row to the big.matrix
# num_variants <- ncol(obj.bigSNP$genotypes)
# num_samples <- nrow(obj.bigSNP$genotypes)
# 
# # Define the path for the new backing file
# backingfile_path <- "~/data/post-qc/combined_genotypes.bk"  # Update this path as necessary
# 
# # Create the FBM with an extra row for the new genotype
# new_fbm <- bigstatsr::FBM(nrow = num_samples + 1, ncol = num_variants, backingfile = backingfile_path)
# new_fbm[1:num_samples, ] <- obj.bigSNP$genotypes[]
# new_fbm[num_samples + 1, ] <- genotype_vector
# 
# # Now 'genotypes' includes your new data, and you can reassign it to your bigSNP object
# # Path to the .bed file (update the path as necessary)
# bedfile_path <- "~/data/post-qc/combined_genotypes.bed"
# 
# # Use the snp_writeBed function to write the FBM object to a bed file
# snp_writeBed(new_fbm, bedfile_path)

# obj.bigSNP$genotypes <- new_fbm
# Save the updated object if needed
# obj.bigSNP$genotypes$save(file.path("~/data/post-qc/", "EUR.QC.bk"))



# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats_munged, map, match.min.prop = 0)
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = "~/repos/1000-genomes-genetic-maps/interpolated_OMNI/")
# calculate LD
for (chr in 1:22) {
  print(chr)
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = NCORES,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld,
                    length(ld),
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff,
                    blocks = NULL)
h2_est <- ldsc[["h2"]]

# Get adjusted beta from the auto model
multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
  ncores = NCORES
)

beta_auto <- sapply(multi_auto, function(auto)
  auto$beta_est)

genotype <- obj.bigSNP$genotypes
# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
pred_auto <-
  big_prodMat(genotype,
              beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)
# scale the PRS generated from AUTO
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-
  rowMeans(beta_auto[,
                     abs(pred_scaled -
                           median(pred_scaled)) <
                       3 * mad(pred_scaled)])
pred_auto <-
  big_prodVec(genotype,
              final_beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)


scores <- pred_auto / sd(pred_auto)
scores <- scores - mean(scores)
hist(scores, breaks = 50, main = "EA Scores")
abline(v = scores[1], col = 2, lwd = 5)
mean(scores[1] > scores)
