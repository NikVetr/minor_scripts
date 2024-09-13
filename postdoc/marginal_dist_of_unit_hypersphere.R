# Load necessary library
library(ggplot2)

# Function to generate points on a hypersphere
generate_hypersphere_points <- function(dim, n_points) {
  # Generate random points from normal distribution
  points <- matrix(rnorm(n_points * dim), ncol = dim)
  # Normalize each point to lie on the unit hypersphere
  points <- points / sqrt(rowSums(points^2))
  return(points)
}

# Function to calculate the theoretical Beta-related PDF
beta_related_pdf <- function(x, k) {
  coefficient <- gamma(k / 2) / (sqrt(pi) * gamma((k - 1) / 2))
  return(coefficient * (1 - x^2) ^ ((k - 3) / 2))
}

# Parameters
dim <- 10  # Dimension of the hypersphere
n_points <- 1E6  # Number of points to sample

# Generate points on the hypersphere
points <- generate_hypersphere_points(dim, n_points)

# Consider one coordinate, say x1 (first column)
x1 <- points[, 1]

# Create a dataframe for ggplot
df <- data.frame(x1 = x1)

# Compute the theoretical PDF values
x_vals <- seq(-1, 1, length.out = 1000)
pdf_vals <- beta_related_pdf(x_vals, dim)
theoretical_df <- data.frame(x = x_vals, y = pdf_vals)

# Plot using ggplot2
ggplot(df, aes(x = x1)) +
  geom_histogram(aes(y = ..density..), bins = 100, color = "black", fill = "green", alpha = 0.6) +
  geom_line(data = theoretical_df, aes(x = x, y = y), color = "blue", size = 1) +
  ggtitle(paste("Empirical vs Theoretical Distribution of x1 from", dim, "D Hypersphere")) +
  xlab("x1") +
  ylab("Density") +
  theme_minimal()
