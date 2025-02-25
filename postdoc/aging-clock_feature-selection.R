# -----------------------------
# Simulation Setup
# -----------------------------
set.seed(123)

# Number of individuals, tissues, and genes
n <- 300
k <- 3      # e.g., heart, muscle, blood
p <- 100     # genes per tissue (for demonstration)

# True effect sizes
prop_age_genes <- 0.1   # Proportion of genes truly affected by age
prop_dis_genes <- 0.1   # Proportion of genes truly affected by disease

# Correlation structure
# We'll create a block correlation among (k * p) "variables"
# to mimic "stronger" correlation across tissues, "moderate" across genes.

# 1) Create correlation blocks by tissue
#    Each tissue: p genes with moderate correlation among them
# 2) Then induce correlation across tissues to make it "strong"
#    so that the same gene across different tissues is more highly correlated.

library(MASS)  # for mvrnorm

# Build a block correlation matrix
rho_within_tissue <- 0.4  # moderate correlation within a tissue's genes
rho_across_tissue <- 0.7  # stronger correlation across the same gene in different tissues

# Initialize a (k*p x k*p) correlation matrix
Sigma <- matrix(0, nrow = k*p, ncol = k*p)
rownames(Sigma) <- colnames(Sigma) <- paste0("T", rep(1:k, each=p), "_G", rep(1:p, k))

# We'll define the diagonal blocks (genes within a single tissue)
for(t in seq_len(k)) {
  block_idx <- ((t-1)*p + 1):(t*p)
  # Fill diagonal block with moderate correlation
  Sigma[block_idx, block_idx] <- rho_within_tissue
  diag(Sigma[block_idx, block_idx]) <- 1
}

# Now induce cross-tissue correlation for the same gene
for(g in seq_len(p)) {
  # For each gene g, we find the positions across tissues
  gene_positions <- seq(g, k*p, by=p)
  # Increase correlation among these positions
  for(i in gene_positions) {
    for(j in gene_positions) {
      if(i != j) {
        Sigma[i, j] <- rho_across_tissue
        Sigma[j, i] <- rho_across_tissue
      }
    }
  }
}

# Ensure Sigma is positive-definite (it should be, but just in case)
# A small eigenvalue fix or nearPD step could be applied if needed
# library(Matrix); Sigma <- as.matrix(nearPD(Sigma)$mat)

# -----------------------------
# Simulate Individual Covariates
# -----------------------------
# Age ~ Uniform(40, 80)
age <- runif(n, 20, 80)

# Disease risk depends on age: logit(prob) = beta0 + beta1 * age
beta0 <- -12  # baseline intercept
beta1 <- 0.15 # effect of age
logit_prob <- beta0 + beta1 * age
disease_prob <- 1 / (1 + exp(-logit_prob))
disease <- rbinom(n, 1, disease_prob)

# -----------------------------
# Determine which genes truly have aging/disease effects
# -----------------------------
# "Sparse" means only a fraction of genes (prop_age_genes, prop_dis_genes) are affected.
# We'll pick them randomly, though you can also fix them for reproducibility.

num_age_genes <- round(k*p * prop_age_genes)
num_dis_genes <- round(k*p * prop_dis_genes)

# Indices of genes with an age effect
age_gene_idx <- sample(k*p, num_age_genes)
# Indices of genes with a disease effect
dis_gene_idx <- c(sample(age_gene_idx, num_dis_genes / 2), sample(k*p, num_dis_genes / 2))
shared_idx <- intersect(age_gene_idx, dis_gene_idx)

# True effect sizes
# For demonstration, make them modest in magnitude
age_effects <- rep(0, k*p)
age_effects[age_gene_idx] <- runif(num_age_genes, 0.05, 0.15)

dis_effects <- rep(0, k*p)
dis_effects[dis_gene_idx] <- runif(num_dis_genes, 0.05, 0.15)

shared_effects <- rep(0, k*p)
shared_effects[shared_idx] <- runif(length(shared_idx), 0.05, 0.15)
age_effects <- age_effects + shared_effects
dis_effects <- dis_effects + shared_effects

# -----------------------------
# Simulation of "Baseline" Expression
# -----------------------------
# We'll sample a baseline multivariate normal vector for each individual (k*p dimension),
# then add the linear age/disease effects.

# The random baseline for all k*p genes in each individual
# We assume mean = 0 for demonstration
baseline_draws <- mvrnorm(n = n, mu = rep(0, k*p), Sigma = Sigma)

# Now create the final expression matrix: X[i, t, g]
# We'll store it as a matrix (n x (k*p)) or an array. Here, a matrix for simplicity:
expression_matrix <- matrix(NA, nrow = n, ncol = k*p)

for(i in seq_len(n)) {
  # linear predictor for each gene-tissue dimension
  linpred <- baseline_draws[i, ] +
    age[i] * age_effects +
    disease[i] * dis_effects
  
  expression_matrix[i, ] <- linpred
}

# Optionally add random measurement noise
# expression_matrix <- expression_matrix + matrix(rnorm(n*k*p, 0, 0.1), n, k*p)

# -----------------------------
# Assemble into a Tidy Structure (Optional)
# -----------------------------
# It's often convenient to put this in a long data.frame
library(dplyr)
library(tidyr)

# Let's define a helper that maps column index to tissue & gene
df <- as.data.frame(expression_matrix)
df$id <- seq_len(n)
long_df <- df %>%
  pivot_longer(
    -id,
    names_to = "var",
    values_to = "expr"
  ) %>%
  mutate(
    age = age[id],
    disease = disease[id]
  )

# Parse tissue & gene labels from var
long_df <- long_df %>%
  separate(var, into = c("tissue_gene"), sep = "_", remove = FALSE) %>%
  # 'tissue_gene' is like "T1" or "T2" etc. for the label we defined earlier
  separate(tissue_gene, into = c("tissue", "gene"), sep = "G", remove = FALSE) %>%
  mutate(
    tissue = gsub("T", "", tissue),
    tissue = as.numeric(tissue),
    gene   = as.numeric(gene)
  ) %>%
  select(-var, -tissue_gene)

# Now long_df has columns:
#   id, expr, age, disease, tissue, gene
# This is workable for fitting models.

# -----------------------------
# Demonstration / Analysis
# -----------------------------
# For instance, we might naïvely regress expression on age ignoring disease,
# or properly include disease, or compare simple vs. multi-level models, etc.

# Example: Fit a naive linear model ignoring disease for Tissue=1, Gene=1
naive_model <- lm(expr ~ age, data = subset(long_df, tissue == 1 & gene == 1))
summary(naive_model)

# Example: Fit a "true" model including disease
true_model <- lm(expr ~ age + disease, data = subset(long_df, tissue == 1 & gene == 1))
summary(true_model)

# You could compare bias in estimated age coefficients for all genes/tissues
# ignoring vs. including disease or random effects. This helps reveal confounding.


# -----------------------------
# Next Steps
# -----------------------------
# 1) Scale up to realistic k, p (maybe 10+ tissues, 1000+ genes).
# 2) Extend the correlation structure or effect definitions (e.g. factor models).
# 3) Fit multi-level or penalized models to test how well they recover
#    the true age effects when disease is a confounder.
# 4) Evaluate partial correlation, partial R^2, etc. to see how ignoring disease
#    or unobserved confounders biases the "aging" coefficients.

