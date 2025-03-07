# 1) Define the transformation function
transform_points <- function(x, orig_span, orig_range, relative_range) {
  new_range <- orig_span[1] + relative_range * diff(orig_span)
  mapped_x <- new_range[1] + (x - orig_range[1]) * (diff(new_range) / diff(orig_range))
  return(mapped_x)
}

# 2) Example usage and plotting

# Define an original span of [0, 10]
orig_span <- c(0, 10)

# Suppose we have an original sub-range of [3, 6]
orig_range <- c(4, 9)

# Suppose we want the new sub-range to occupy relative positions 0.2 and 0.6
# within the full span [0, 1], i.e. 20% to 60% across the 0..10 segment
relative_range <- c(0.2, 0.6)

# Create a set of points x in the original span
x <- seq(orig_range[1], orig_range[2], length.out = 10)

# Transform them
x_mapped <- transform_points(x, orig_span, orig_range, relative_range)

# For reference, compute the actual new_range:
new_range <- orig_span[1] + relative_range * diff(orig_span)

# Now let's make a simple plot
plot(
  x, 
  rep(0, length(x)), 
  pch = 16, xlim = orig_span,
  xlab = "Value", 
  ylab = "Illustration Only", 
  main = "Transforming Points from orig_range to new_range"
)
abline(v = orig_span[1]:orig_span[2], lty = 3)

# Show original x points in blue
points(x, rep(0, length(x)), pch = 16, col = "blue")

# Draw a horizontal line marking orig_range
segments(orig_range[1], 0, orig_range[2], 0, lwd = 2, col = "blue")

# Show mapped x points in red (shift them vertically for clarity)
points(x_mapped, rep(1, length(x_mapped)), pch = 16, col = "red")

# Draw a horizontal line marking new_range
segments(new_range[1], 1, new_range[2], 1, lwd = 2, col = "red")

# Add a legend
legend(
  "bottom", 
  legend = c("Original points", "Mapped points"), 
  col = c("blue", "red"), 
  pch = 16, 
  horiz = TRUE
)
