# Load necessary libraries
library(grDevices)
library(grid)

# Read the input text file
font_family <- hershey_fonts <- c(
  "HersheyGothicEnglish",
  "HersheyGothicGerman",
  "HersheyScript",
  "HersheyPlain",
  "HersheySans",
  "HersheySerif",
  "HersheySerifItalic",
  "HersheySerifBold",
  "HersheySerifBoldItalic"
)[3]

input_file <- "~/therapist_ty.txt"  # Path to your input text file
txtlines <- trimws(readLines(input_file, warn = FALSE))
output_file <- "~/therapist_ty.svg"  # Path to your output SVG file
aspect_ratio <- 0.7  # Example aspect ratio

# Calculate width based on aspect ratio (assuming fixed height)

line_dims <- list(w = strwidth(txtlines, units = "inches"), h = strheight(txtlines, units = "inches"))
total_area <- sum(line_dims$w * line_dims$h) * 4
height <- sqrt(total_area / aspect_ratio)
width <- total_area / height
ncharwidth <- width / max(strwidth(LETTERS, units = "inches")) * 1.5
# Wrap the text to fit the calculated width
wrapped_text <- strwrap(txtlines, width = ncharwidth)  

# Create an SVG file
svg(output_file, width = width, height = height)

# Set up the graphics parameters to use Hershey font
par(family = font_family)

# Calculate the height of a line of text
line_height <- convertUnit(unit(1, "lines"), "npc", valueOnly = TRUE) * aspect_ratio

# Starting y position
y_start <- 1 - line_height / 2  # Start near the top, slightly adjusted

# Plot each line of the wrapped text individually
# grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, just = "center", gp = gpar(lty = "dashed"))
y_pos <- y_start - seq(0, length(wrapped_text) - 1) * line_height
for (i in seq_along(wrapped_text)) {
  grid.text(wrapped_text[i], x = 0.01, y = y_pos[i], 
            just = "left", gp = gpar(fontfamily = font_family))
}



# Close the SVG device
dev.off()
