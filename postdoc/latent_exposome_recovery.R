# Load necessary libraries
library(MASS)   # For multivariate normal simulations
library(ggplot2) # For plotting
library(reshape2) # For data manipulation
library(parallel)

# Set seed for reproducibility
set.seed(123)

# Simulation parameters
n_individuals <- 5000   # Number of individuals
n_snps <- 2000         # Number of SNPs
n_traits <- 100          # Number of traits
n_env_factors <- 20      # Number of latent environmental factors

# Simulate genotype matrix (assuming SNPs are unlinked and biallelic)
genotype_matrix <- matrix(rbinom(n_individuals * n_snps, 2, 0.5),
                          nrow = n_individuals, ncol = n_snps)

# Simulate true SNP effects for each trait
beta_matrix <- matrix(rnorm(n_snps * n_traits, mean = 0, sd = 0.1),
                      nrow = n_snps, ncol = n_traits)

# Simulate latent environmental factors
env_factors <- matrix(rnorm(n_individuals * n_env_factors),
                      nrow = n_individuals, ncol = n_env_factors)

# Simulate environmental loadings for each trait
env_loadings <- matrix(rnorm(n_env_factors * n_traits, mean = 0, sd = 1),
                       nrow = n_env_factors, ncol = n_traits)

# Simulate error term
residuals <- matrix(rnorm(n_individuals * n_traits, mean = 0, sd = 0),
                    nrow = n_individuals, ncol = n_traits)

# Generate phenotypes
genetic_effects <- genotype_matrix %*% beta_matrix
environmental_effects <- env_factors %*% env_loadings

apply(genetic_effects, 2, var)
apply(environmental_effects, 2, var)
phenotypes <- genetic_effects + environmental_effects + residuals

# Split data into training and testing sets for GWAS
train_indices <- sample(1:n_individuals, size = n_individuals / 4 * 3)
test_indices <- setdiff(1:n_individuals, train_indices)

# Perform GWAS for each trait on training data
mcprint <- function(...){system(sprintf('printf "%s"', paste0(..., collapse="")))}
gwas_results <- mclapply(1:n_traits, function(t) {
  mcprint(paste0(t, " "))
  phenotype_train <- phenotypes[train_indices, t]
  genotype_train <- genotype_matrix[train_indices, ]
  # Fit linear model for each SNP (univariate GWAS)
  betas <- apply(genotype_train, 2, function(snp) {
    coef(lm(phenotype_train ~ snp))[2]
  })
  return(betas)
}, mc.cores = 12)

# Compute polygenic scores for test data
polygenic_scores_test <- matrix(0, nrow = length(test_indices), ncol = n_traits)
for (t in 1:n_traits) {
  genotype_test <- genotype_matrix[test_indices, ]
  betas <- gwas_results[[t]]
  polygenic_scores_test[, t] <- genotype_test %*% betas
}

polygenic_scores_train <- matrix(0, nrow = length(train_indices), ncol = n_traits)
for (t in 1:n_traits) {
  genotype_train <- genotype_matrix[train_indices, ]
  betas <- gwas_results[[t]]
  polygenic_scores_train[, t] <- genotype_train %*% betas
}

# Get observed phenotypes
phenotypes_test <- phenotypes[test_indices, ]
phenotypes_train <- phenotypes[train_indices, ]

# Calculate residuals
residuals_train <- phenotypes_train - polygenic_scores_train

# Perform PCA on residuals
pca_residuals_train <- prcomp(residuals_train, center = TRUE, scale. = TRUE)

# Extract the top PCs (as many as the number of true env factors)
var_prop <- 0.5
n_PCs <- sum(cumsum(pca_residuals_train$sdev^2) < var_prop * n_traits)
estimated_env_factors <- pca_residuals_train$x[, 1:n_PCs]

# Compare estimated env factors with true env factors
# For comparison, we'll compute the correlation between estimated PCs and true env factors

true_env_factors_train <- env_factors[train_indices, ]

# Compute correlation matrix
cor_matrix <- cor(estimated_env_factors, true_env_factors_train)
apply(cor_matrix, 2, which.max)
apply(cor_matrix, 2, max)
apply(cor_matrix, 2, sort, T)

# Refit GWAS including estimated environmental PCs as covariates
adjusted_gwas_results <- mclapply(1:n_traits, function(t) {
  mcprint(paste0(t, " "))
  phenotype_train <- phenotypes[train_indices, t]
  genotype_train <- genotype_matrix[train_indices, ]
  # Fit linear model for each SNP including PCs
  betas <- apply(genotype_train, 2, function(snp) {
    model <- lm(phenotype_train ~ snp + estimated_env_factors)
    coef(model)[2]
  })
  return(betas)
}, mc.cores = 12)

# Evaluate whether including estimated env factors reduces confounding
# Here, we could compare the variance explained or the accuracy of SNP effect estimates

# For demonstration, let's plot SNP effect estimates before and after adjustment for one trait
trait_to_plot <- 1
betas_unadjusted <- gwas_results[[trait_to_plot]]
betas_adjusted <- adjusted_gwas_results[[trait_to_plot]]
betas_true <- beta_matrix[, trait_to_plot]

df_betas <- data.frame(
  SNP = 1:n_snps,
  True_Effect = betas_true,
  Unadjusted_Effect = betas_unadjusted,
  Adjusted_Effect = betas_adjusted
)

# Plotting true vs estimated SNP effects
library(ggplot2)
ggplot(df_betas, aes(x = True_Effect, y = Unadjusted_Effect)) +
  geom_point(alpha = 0.5) +
  labs(title = "Unadjusted SNP Effects vs True Effects", x = "True SNP Effect", y = "Estimated SNP Effect (Unadjusted)") +
  geom_abline(slope = 1, intercept = 0, color = "red")

ggplot(df_betas, aes(x = True_Effect, y = Adjusted_Effect)) +
  geom_point(alpha = 0.5) +
  labs(title = "Adjusted SNP Effects vs True Effects", x = "True SNP Effect", y = "Estimated SNP Effect (Adjusted)") +
  geom_abline(slope = 1, intercept = 0, color = "red")

# Compute correlation between true and estimated SNP effects
cor_unadjusted <- cor(betas_true, betas_unadjusted)
cor_adjusted <- cor(betas_true, betas_adjusted)

print(paste("Correlation between true and unadjusted SNP effects:", round(cor_unadjusted, 3)))
print(paste("Correlation between true and adjusted SNP effects:", round(cor_adjusted, 3)))

#compute polygenic scores
# Compute polygenic scores for test data
adj_polygenic_scores_test <- matrix(0, nrow = length(test_indices), ncol = n_traits)
for (t in 1:n_traits) {
  genotype_test <- genotype_matrix[test_indices, ]
  betas <- adjusted_gwas_results[[t]]
  adj_polygenic_scores_test[, t] <- genotype_test %*% betas
}

genetic_effects_test <- genetic_effects[test_indices,]
plot(x = sapply(1:n_traits, function(t) cor(genetic_effects_test[,t], polygenic_scores_test[,t])),
     y = sapply(1:n_traits, function(t) cor(genetic_effects_test[,t], adj_polygenic_scores_test[,t])),
     xlab = "unadjusted correlations", ylab = "adjusted correlations", main = "Pearson Correlation of Polygenic Score\nwith Genetic Liability")
abline(0,1, col=2, lty = 2, lwd = 3)
