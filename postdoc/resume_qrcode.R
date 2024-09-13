library(png)
library(magick)
library(qrcode)
library(corrplot)

source("~/repos/polylines/R/functions.R")

#### make logo ####
draw_NIK <- F
generate_QR <- F

# svg("~/nik_name_logo.svg", width = 400/72, height = 400/72)

# Define the parameters for the plot
lwd <- 16
n_strand <- 1000
dh_bounds <- c(-1.5*pi, 1.5*pi)
t <- seq(dh_bounds[1], dh_bounds[2], length.out = n_strand)   # Position along the strand for a single twist
amplitude <- 1                                 # Amplitude of the sinusoids
phase_shift <- pi                              # Phase shift for the second sinusoid
frequency <- 1                                 # Frequency of the sinusoids
base_pair_interval <- pi / 6                   # Base pair every pi/4 units
proportional_distance <- 0.2                  # Proportional distance from the opposite strand
draw_internal_strand <- F
alternate_sides <- F
partial_displacement <- T
chirality_prop <- 0.4
prop_square_wave <- 0.2
bp_thickness <- 1.25
strand_thickness <- 0.1
extend_straight_by <- pi*2
straight_y <- seq(dh_bounds[2], dh_bounds[2] + 
                    extend_straight_by, by = diff(dh_bounds)/n_strand)
nstraight <- length(straight_y)
rot <- 90

# Creating strand_df
strand_df <- data.frame(
  y = t,
  x1 = amplitude * sin(frequency * t),         # First strand x-coordinates
  x2 = amplitude * sin(frequency * t + phase_shift),  # Second strand x-coordinates
  z1 = amplitude * sin(frequency * t + pi/2),         # First strand x-coordinates
  z2 = amplitude * sin(frequency * t + phase_shift + pi/2)  # Second strand x-coordinates
)

strand_df <- rbind(strand_df, 
                   data.frame(
                     y = straight_y,
                     x1 = rep(strand_df$x1[n_strand], nstraight),
                     x2 = rep(strand_df$x2[n_strand], nstraight),
                     z1 = rep(strand_df$z1[n_strand], nstraight),
                     z2 = rep(strand_df$z2[n_strand], nstraight)
                   ))

bounds <- range(strand_df$y)
maxw_strand <- range(strand_df$x1)

#distort away from perfect sine wave
# strand_df$x1 <- strand_df$x1 * asinh(abs(strand_df$x1)^-prop_square_wave)
# strand_df$x2 <- strand_df$x2 * asinh(abs(strand_df$x2)^-prop_square_wave)
strand_df$x1 <- abs(strand_df$x1)^(1-prop_square_wave) * sign(strand_df$x1) * prop_square_wave + 
  (1-prop_square_wave) * strand_df$x1
strand_df$x2 <- abs(strand_df$x2)^(1-prop_square_wave) * sign(strand_df$x2) * prop_square_wave + 
  (1-prop_square_wave) * strand_df$x2

#visual width of the strand
strand_df$lwd1 <- (rescale01(strand_df$z1) * sqrt(strand_thickness))^2
strand_df$lwd2 <- (rescale01(strand_df$z2) * sqrt(strand_thickness))^2 

# Maximum width between the strands and whitespace to be left over
max_width <- max(abs(strand_df$x2 - strand_df$x1))
ws_strand <- proportional_distance * max_width

# Base pair indices
bp_indices <- sapply(seq(bounds[1], bounds[2], by = base_pair_interval), function(pos) which.min(abs(strand_df$y - pos)))

# Creating bp_df
bp_df <- data.frame(
  y0 = strand_df$y[bp_indices],
  x0 = strand_df$x1[bp_indices],
  y1 = strand_df$y[bp_indices],
  x1 = rep(NA, length(bp_indices)),
  w = rep(NA, length(bp_indices)),
  dir = rep(NA, length(bp_indices)),
  lwd0 = strand_df$lwd1[bp_indices],
  lwd1 = strand_df$lwd2[bp_indices]
)

# Fill in w and dir in bp_df
bp_df$w <- pmax(abs(strand_df$x2[bp_indices] - strand_df$x1[bp_indices]) - ws_strand, 0)
bp_df$dir <- ifelse(bp_df$x0 >= 0, -1, 1)  # -1 for left, 1 for right

# Correctly calculate x1 based on w and dir
bp_df$x1 <- bp_df$x0 + (bp_df$w * bp_df$dir)

# Ensure x1 values do not exceed the strand limits
bp_df$x1 <- pmin(pmax(bp_df$x1, min(strand_df$x1, strand_df$x2)), max(strand_df$x1, strand_df$x2))

#if alternating sides, displace x values by the ws remainging
if(alternate_sides){
  bp_df$x0 <- bp_df$x0 + 
    (ws_strand * as.numeric(bp_df$dir==1) - 
    ws_strand * as.numeric(bp_df$dir==-1)) * rep(0:1, length.out = nrow(bp_df)) / 
    ifelse(partial_displacement, 2, 1)
  bp_df$x1 <- bp_df$x1 + 
    (ws_strand * as.numeric(bp_df$dir==1) - 
    ws_strand * as.numeric(bp_df$dir==-1)) * rep(0:1, length.out = nrow(bp_df)) / 
    ifelse(partial_displacement, 2, 1)
}

#remove 0 width bases
bp_df <- bp_df[bp_df$w > 1E-6,] 

#adjust for chirality
strand_df$seen <- abs(strand_df$x1 - strand_df$x2) > chirality_prop
on_top_1 <- diff(strand_df$x1) > 0
strand_df$on_top_1 <- c(on_top_1[1], on_top_1)
strand_df$seen_1 <- strand_df$seen_2 <- T
strand_df$seen_1[!strand_df$seen & !strand_df$on_top_1] <- F
strand_df$seen_2[!strand_df$seen & strand_df$on_top_1] <- F

#rotate if requested
center <- c(mean(c(strand_df$x1, strand_df$x2)), mean(strand_df$y), mean(c(strand_df$z1, strand_df$z2)))
rotmat <- rotmat_00(rot)
strand_df$

# Plotting
png("~/nik_name_logo.png", width = 400 * diff(bounds) / pi, height = 400 * diff(bounds) / pi)
par(lwd = lwd, xpd = NA, mar = c(0,0,0,0))

if(draw_NIK){
  plot(NULL, NULL, xlim = c(-3.5, 3), ylim = c(-3.25, 3.25), frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")  
} else {
  plot(NULL, NULL, xlim = c(-diff(bounds)/2, diff(bounds)/2), ylim = bounds, frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
}

seen_sets_1 <- split(which(strand_df$seen_1), cumsum(!strand_df$seen_1)[strand_df$seen_1])
seen_sets_2 <- split(which(strand_df$seen_2), cumsum(!strand_df$seen_2)[strand_df$seen_2])

ss1_n <- sapply(seen_sets_1, length)
ss1_thresh_n <- mean(range(ss1_n))
ss2_n <- sapply(seen_sets_2, length)
ss2_thresh_n <- mean(range(ss2_n))

#strand 1
for(i in seq_along(seen_sets_1)){
  seen_set <- seen_sets_1[[i]]
  sdf <- strand_df[seen_set,]
  polylines(sdf$x1, sdf$y, lwd = sdf$lwd1, col = "black")
}

#strand 2
for(i in seq_along(seen_sets_2)){
  seen_set <- seen_sets_2[[i]]
  sdf <- strand_df[seen_set,]
  polylines(sdf$x2, sdf$y, lwd = sdf$lwd2, col = "black")
}

# with(bp_df, segments(x0, y0, x1, y1, col = "black"))
nbp <- n_strand / nrow(bp_df) 
for(i in 1:nrow(bp_df)){
  segments(bp_df$x0[i], bp_df$y0[i], bp_df$x1[i], bp_df$y1[i], col = "black")
  xbp <- seq(bp_df$x0[i], bp_df$x1[i], length.out = nbp)
  ybp <- seq(bp_df$y0[i], bp_df$y1[i], length.out = nbp)
  
  #length and direction tell us how wide the nucs are
  lwd <- seq(bp_df$lwd0[i], 
             bp_df$lwd1[i] - (bp_df$lwd1[i] - bp_df$lwd0[i]) * proportional_distance, 
             length.out = nbp)
  lwd <- rescale01(lwd)^2 * diff(range(lwd)) + min(lwd)
  lwd <- lwd * bp_thickness
  polylines(x = xbp, y = ybp, lwd = lwd)
}

#connect internal bps
if(draw_internal_strand){
  intersection_index <- which.min(abs(strand_df$x1 - strand_df$x2 - ws_strand))
  lines(strand_df$x2[intersection_index:n_strand] + ws_strand, strand_df$y[intersection_index:n_strand], col = "black")
  lines(strand_df$x2[1:(n_strand-intersection_index)] - ws_strand, strand_df$y[1:(n_strand-intersection_index)], col = "black")
}

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



#### logo with just the DNA strand, pixelated ####

# Read and preprocess the image
img <- readPNG("~/nik_name_logo.png")

# Convert the image to a matrix
target_dim <- 17
img_matrix <- t(img[,,1] + img[,,2] + img[,,3] < 0.1)
black_pix <- ceiling(which(img_matrix, arr.ind = T) / max(dim(img_matrix)) * target_dim)
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
if(generate_QR){
  out <- qr_code("https://github.com/NikVetr/CV/blob/master/vetr_resume_1pg.pdf", ecl = "H")
  out <- add_logo(out, logo = "~/nik_name_logo_pixelated.svg", ecl = "L")
  qrcode::generate_svg(out, "~/resume_qr_code.svg")
}
cairo_pdf("~/resume_qr_code.pdf", width = 1000, height = 1000)
plot(out)
dev.off()
