library(png)
library(magick)
library(qrcode)
library(corrplot)

#### make logo ####
draw_NIK <- F

# svg("~/nik_name_logo.svg", width = 400/72, height = 400/72)
png("~/nik_name_logo.png", width = 400, height = 400)
par(lwd = 4, xpd = NA, mar = c(0,0,0,0))
# Define the parameters for the plot
n_strand <- 1000
t <- seq(-pi / 2, pi / 2, length.out = n_strand)   # Position along the strand for a single twist
amplitude <- 1                                 # Amplitude of the sinusoids
phase_shift <- pi                              # Phase shift for the second sinusoid
frequency <- 1                                 # Frequency of the sinusoids
base_pair_interval <- pi / 8                   # Base pair every pi/4 units
proportional_distance <- 0                  # Proportional distance from the opposite strand

# Creating strand_df
strand_df <- data.frame(
  y = t,
  x1 = amplitude * sin(frequency * t),         # First strand x-coordinates
  x2 = amplitude * sin(frequency * t + phase_shift)  # Second strand x-coordinates
)

# Maximum width between the strands and whitespace to be left over
max_width <- max(abs(strand_df$x2 - strand_df$x1))
ws_strand <- proportional_distance * max_width

# Base pair indices
bp_indices <- sapply(seq(-pi / 2, pi / 2, by = base_pair_interval), function(pos) which.min(abs(t - pos)))

# Creating bp_df
bp_df <- data.frame(
  y0 = strand_df$y[bp_indices],
  x0 = strand_df$x1[bp_indices],
  y1 = strand_df$y[bp_indices],
  x1 = rep(NA, length(bp_indices)),
  w = rep(NA, length(bp_indices)),
  dir = rep(NA, length(bp_indices))
)

# Fill in w and dir in bp_df
bp_df$w <- pmax(abs(strand_df$x2[bp_indices] - strand_df$x1[bp_indices]) - ws_strand, 0)
bp_df$dir <- ifelse(bp_df$y0 >= 0, -1, 1)  # -1 for left, 1 for right

# Correctly calculate x1 based on w and dir
bp_df$x1 <- bp_df$x0 + (bp_df$w * bp_df$dir)

# Ensure x1 values do not exceed the strand limits
bp_df$x1 <- pmin(pmax(bp_df$x1, min(strand_df$x1, strand_df$x2)), max(strand_df$x1, strand_df$x2))

# Plotting
if(draw_NIK){
  plot(NULL, NULL, xlim = c(-3.5, 3), ylim = c(-3.25, 3.25), frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")  
} else {
  plot(NULL, NULL, xlim = c(-1, 1), ylim = c(-1.5, 1.5), frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
}

with(strand_df, {
  lines(x1, y, col = "black")
  lines(x2, y, col = "black")
})
with(bp_df, segments(x0, y0, x1, y1, col = "black"))

#connect internal bps
intersection_index <- which.min(abs(strand_df$x1 - strand_df$x2 - ws_strand))
lines(strand_df$x2[intersection_index:n_strand] + ws_strand, strand2_y[intersection_index:n_strand], col = "black")
lines(strand_df$x2[1:(n_strand-intersection_index)] - ws_strand, strand2_y[1:(n_strand-intersection_index)], col = "black")

if(draw_NIK){
  
# draw a letter "N"

#sinusoid for N
lines(strand_df$x2-2.5, strand_df$y, col = "black")

#vertical bars of "N"
segments(x0 = min(strand_df$x2-2.5), y0 = min(strand_df$y), x1 = min(strand_df$x2-2.5), y1 = max(strand_df$y))
segments(x0 = max(strand_df$x2-2.5), y0 = min(strand_df$y), x1 = max(strand_df$x2-2.5), y1 = max(strand_df$y))

# draw a letter "K"

#vertical bar of "K"
segments(x0 = max(strand_df$x2+0.5), y0 = min(strand_df$y), x1 = max(strand_df$x2+0.5), y1 = max(strand_df$y))

#horizontal bar of "K"
strand_intersect_index <- which.min(abs(strand_df$x1 - strand_df$x2))
segments(x0 = max(strand_df$x2+0.5), y0 = strand_df$y[strand_intersect_index], x1 = max(strand_df$x2+1), y1 = strand_df$y[strand_intersect_index])

#strands of "K"
lines(strand_df$x2[1:strand_intersect_index] + 2, strand_df$y[1:strand_intersect_index], col = "black")
lines(strand_df$x1[strand_intersect_index:n_strand] + 2, strand_df$y[strand_intersect_index:n_strand], col = "black")

}

dev.off()
####

#logo with just the DNA strand, pixelated

# Read and preprocess the image
img <- readPNG("~/nik_name_logo.png")

# Convert the image to a matrix
target_dim <- 17
img_matrix <- t(img[,,1] + img[,,2] + img[,,3] < 0.1)
black_pix <- ceiling(which(img_matrix, arr.ind = T) / nrow(img_matrix) * target_dim)
black_pix <- black_pix[!duplicated(black_pix),]
img_matrix <- matrix(F, target_dim, target_dim)
img_matrix[black_pix] <- T

# Plotting
svg("~/nik_name_logo_pixelated.svg", width = 400 / 72, height = 400 / 72)

par(mar = c(0,0,0,0))
plot(NULL, NULL, xlim = c(1/2,ncol(img_matrix)+1/2), 
     ylim = c(1/2, nrow(img_matrix)+1/2), frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "", asp = 1, xpd = NA)

# Draw rectangles for each non-white pixel
rect(xleft = black_pix[,1] - 1/2 - 1E-2, 
     ybottom = black_pix[,2] - 1/2 - 1E-2, 
     xright = black_pix[,1] + 1/2 + 1E-2, 
     ytop = black_pix[,2] + 1/2 + 1E-2, col = 1, border = NA)

dev.off()


####
#make qc code
out <- qr_code("https://github.com/NikVetr/CV/blob/master/vetr_resume_1pg.pdf", ecl = "H")
out <- add_logo(out, logo = "~/nik_name_logo_pixelated.svg", ecl = "L")
qrcode::generate_svg(out, "~/resume_qr_code.svg")
cairo_pdf("~/resume_qr_code.pdf", width = 1000, height = 1000)
plot(out)
dev.off()
