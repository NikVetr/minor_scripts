source("~/repos/polylines/R/functions.R")

#### function ####

draw_DNA <- function(
    n_strand = 1000,
    dh_bounds = c(-1.5 * pi, 3.5 * pi),
    amplitude = 1,
    phase_shift = pi,
    frequency = 1,
    base_pair_interval = pi / 6, #should really be 5.25 for 
    proportional_distance = 0.2,
    draw_internal_strand = FALSE,
    alternate_sides = FALSE,
    partial_displacement = TRUE,
    chirality_prop = 0.4,
    prop_square_wave = 0.2,
    bp_thickness = 1,
    strand_thickness = 0.1,
    extend_straight = "n",
    extend_straight_by = pi * 2,
    rot = 0,
    target_center = c(0,0),
    box_dim = c(10,10),
    col = 1,
    ...
) {
  
  # initializing the main data frame
  xyr <- xyrat()
  t <- seq(dh_bounds[1], dh_bounds[2], length.out = n_strand)   # Position along the strand for a single twist
  
  strand_df <- data.frame(
    y1 = t,
    y2 = t,
    x1 = amplitude * sin(frequency * t) * xyr,         # First strand x-coordinates
    x2 = amplitude * sin(frequency * t + phase_shift) * xyr,  # Second strand x-coordinates
    z1 = amplitude * sin(frequency * t + pi/2),         # First strand x-coordinates
    z2 = amplitude * sin(frequency * t + phase_shift + pi/2)  # Second strand x-coordinates
  )
  
  #add in straight segment to one end
  
  #calculate straight points
  straight_y <- seq(dh_bounds[2], dh_bounds[2] + 
                      extend_straight_by, by = diff(dh_bounds)/n_strand)
  nstraight <- length(straight_y)
  
  if(extend_straight != "n"){
    strand_df <- rbind(strand_df, 
                       data.frame(
                         y1 = straight_y,
                         y2 = straight_y,
                         x1 = rep(strand_df$x1[n_strand], nstraight),
                         x2 = rep(strand_df$x2[n_strand], nstraight),
                         z1 = rep(strand_df$z1[n_strand], nstraight),
                         z2 = rep(strand_df$z2[n_strand], nstraight)
                       ))
  }
  
  if(extend_straight == "1" || extend_straight == 1){
    strinds <- (n_strand+1):(n_strand+nstraight)
    strand_df$y2[strinds] <- strand_df$y2[n_strand]
  }
  
  if(extend_straight == "2" || extend_straight == 2){
    strinds <- (n_strand+1):(n_strand+nstraight)
    strand_df$y1[strinds] <- strand_df$y1[n_strand]
  }
  
  #compute high level parameters
  bounds <- range(c(strand_df$y1, strand_df$y2))
  maxw_strand <- range(strand_df$x1)
  
  #distort away from perfect sine wave
  # strand_df$x1 <- strand_df$x1 * asinh(abs(strand_df$x1)^-prop_square_wave)
  # strand_df$x2 <- strand_df$x2 * asinh(abs(strand_df$x2)^-prop_square_wave)
  strand_df$x1 <- abs(strand_df$x1)^(1-prop_square_wave) * sign(strand_df$x1) * prop_square_wave + 
    (1-prop_square_wave) * strand_df$x1
  strand_df$x2 <- abs(strand_df$x2)^(1-prop_square_wave) * sign(strand_df$x2) * prop_square_wave + 
    (1-prop_square_wave) * strand_df$x2
  
  #visual width of the strand
  strand_df$lwd1 <- (rescale01(strand_df$z1) * sqrt(strand_thickness))^2 + strand_thickness / 2
  strand_df$lwd2 <- (rescale01(strand_df$z2) * sqrt(strand_thickness))^2 + strand_thickness / 2
  
  # Maximum width between the strands and whitespace to be left over
  max_width <- max(abs(strand_df$x2 - strand_df$x1))
  ws_strand <- proportional_distance * max_width
  
  # Base pair indices
  if(extend_straight == "b" || extend_straight == "n"){
    bp_indices <- sapply(seq(bounds[1], bounds[2], by = base_pair_interval), 
                         function(pos) which.min(abs(strand_df$y1 - pos)))
  } else {
    bp_indices <- sapply(seq(bounds[1], bounds[2], by = base_pair_interval), 
                         function(pos) which.min(abs(t - pos)))
  }
  
  # Creating bp_df
  bp_df <- data.frame(
    y0 = strand_df$y1[bp_indices],
    x0 = strand_df$x1[bp_indices],
    y1 = strand_df$y2[bp_indices],
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
  center <- c(x = mean(c(strand_df$x1, strand_df$x2)), 
              y = mean(c(strand_df$y1, strand_df$y2)), 
              z = mean(c(strand_df$z1, strand_df$z2)))
  rotmat <- rotmat_00(rot)
  strand_df[,c("x1", "y1")] <- t(t(as.matrix(t(t(strand_df[,c("x1", "y1")]) - 
                                                 center[c("x", "y")])) %*% rotmat) + target_center)
  strand_df[,c("x2", "y2")] <- t(t(as.matrix(t(t(strand_df[,c("x2", "y2")]) - 
                                                 center[c("x", "y")])) %*% rotmat) + target_center)
  bp_df[,c("x0", "y0")] <- t(t(as.matrix(t(t(bp_df[,c("x0", "y0")]) - 
                                             center[c("x", "y")])) %*% rotmat) + target_center)
  bp_df[,c("x1", "y1")] <- t(t(as.matrix(t(t(bp_df[,c("x1", "y1")]) - 
                                             center[c("x", "y")])) %*% rotmat) + target_center)
  
  
  #rescale to size of bounding box
  min_pts <- c(x = min(c(strand_df$x1, strand_df$x2)), 
               y = min(c(strand_df$y1, strand_df$y2)), 
               z = min(c(strand_df$z1, strand_df$z2)))
  sdf_pts <- strand_df[,c("x1", "y1", "z1", "x2", "y2", "z2")]
  bdf_pts <- bp_df[,c("x0", "y0", "x1", "y1")]
  
  sdf_pts <- t(t(sdf_pts) - c(min_pts, min_pts))
  bdf_pts <- t(t(bdf_pts) - c(min_pts[c("x", "y")], min_pts[c("x", "y")]))
  max_pts <- apply(sdf_pts, 2, max)
  max_pts <- c(x = max(max_pts[c("x1", "x2")]), 
               y = max(max_pts[c("y1", "y2")]), 
               z = max(max_pts[c("z1", "z2")]))
  rescale_factor <- min(box_dim / max_pts[c("x", "y")])
  sdf_pts <- t(t(sdf_pts * rescale_factor))
  bdf_pts <- t(t(bdf_pts * rescale_factor))
  
  #recenter these
  center <- c(x = mean(sdf_pts[,c("x1", "x2")]), 
              y = mean(sdf_pts[,c("y1", "y2")]), 
              z = mean(sdf_pts[,c("z1", "z2")]))
  sdf_pts[,c("x1", "y1", "x2", "y2")] <- t(t(sdf_pts[,c("x1", "y1", "x2", "y2")]) - center[c("x", "y", "x", "y")] + 
                                            c(target_center, target_center))
  bdf_pts[,c("x0", "y0", "x1", "y1")] <- t(t(bdf_pts[,c("x0", "y0", "x1", "y1")]) - center[c("x", "y", "x", "y")] + 
                                            c(target_center, target_center))
  
  
  #and reassign
  strand_df[,c("x1", "y1", "z1", "x2", "y2", "z2")] <- sdf_pts
  bp_df[,c("x0", "y0", "x1", "y1")] <- bdf_pts
  
  
  #find where the strands go underneath each other
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
    polylines(sdf$x1, sdf$y1, lwd = sdf$lwd1, col = col, ...)
  }
  
  #strand 2
  for(i in seq_along(seen_sets_2)){
    seen_set <- seen_sets_2[[i]]
    sdf <- strand_df[seen_set,]
    polylines(sdf$x2, sdf$y2, lwd = sdf$lwd2, col = col, ...)
  }
  
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
    polylines(x = xbp, y = ybp, lwd = lwd, col = col, ...)
  }
  
}

png("~/dna_strand.png", width = 2000, height = 2000)
par(xpd = NA, mar = c(0,0,0,0))

plot(NULL, NULL, xlim = c(-8,8), ylim = bounds, frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
draw_DNA(rot = 135, extend_straight = "b", target_center = c(2,1), box_dim = c(8,14))

dev.off()


#### initial non-function-wrapped code ####

# Define the parameters for the plot
lwd <- 16
n_strand <- 1000
dh_bounds <- c(-1.5*pi, 3.5*pi)
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
extend_straight <- c("1", "2", "b", "n")[4]
extend_straight_by <- pi*2
rot <- 270

# initializing the main data frame
t <- seq(dh_bounds[1], dh_bounds[2], length.out = n_strand)   # Position along the strand for a single twist

strand_df <- data.frame(
  y1 = t,
  y2 = t,
  x1 = amplitude * sin(frequency * t),         # First strand x-coordinates
  x2 = amplitude * sin(frequency * t + phase_shift),  # Second strand x-coordinates
  z1 = amplitude * sin(frequency * t + pi/2),         # First strand x-coordinates
  z2 = amplitude * sin(frequency * t + phase_shift + pi/2)  # Second strand x-coordinates
)

#add in straight segment to one end

#calculate straight points
straight_y <- seq(dh_bounds[2], dh_bounds[2] + 
                    extend_straight_by, by = diff(dh_bounds)/n_strand)
nstraight <- length(straight_y)

if(extend_straight != "n"){
  strand_df <- rbind(strand_df, 
                     data.frame(
                       y1 = straight_y,
                       y2 = straight_y,
                       x1 = rep(strand_df$x1[n_strand], nstraight),
                       x2 = rep(strand_df$x2[n_strand], nstraight),
                       z1 = rep(strand_df$z1[n_strand], nstraight),
                       z2 = rep(strand_df$z2[n_strand], nstraight)
                     ))
}

if(extend_straight == "1" || extend_straight == 1){
  strinds <- (n_strand+1):(n_strand+nstraight)
  strand_df$y2[strinds] <- strand_df$y2[n_strand]
}

if(extend_straight == "2" || extend_straight == 2){
  strinds <- (n_strand+1):(n_strand+nstraight)
  strand_df$y1[strinds] <- strand_df$y1[n_strand]
}

#compute high level parameters
bounds <- range(c(strand_df$y1, strand_df$y2))
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
if(extend_straight == "b" || extend_straight == "n"){
    bp_indices <- sapply(seq(bounds[1], bounds[2], by = base_pair_interval), 
                     function(pos) which.min(abs(strand_df$y1 - pos)))
} else {
  bp_indices <- sapply(seq(bounds[1], bounds[2], by = base_pair_interval), 
                       function(pos) which.min(abs(t - pos)))
}

# Creating bp_df
bp_df <- data.frame(
  y0 = strand_df$y1[bp_indices],
  x0 = strand_df$x1[bp_indices],
  y1 = strand_df$y2[bp_indices],
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
center <- c(mean(c(strand_df$x1, strand_df$x2)), mean(c(strand_df$y1, strand_df$y2)), mean(c(strand_df$z1, strand_df$z2)))
rotmat <- rotmat_00(rot)
strand_df[,c("x1", "y1")] <- t(t(as.matrix(t(t(strand_df[,c("x1", "y1")]) - center[1:2])) %*% rotmat) + center[1:2])
strand_df[,c("x2", "y2")] <- t(t(as.matrix(t(t(strand_df[,c("x2", "y2")]) - center[1:2])) %*% rotmat) + center[1:2])
bp_df[,c("x0", "y0")] <- t(t(as.matrix(t(t(bp_df[,c("x0", "y0")]) - center[1:2])) %*% rotmat) + center[1:2])
bp_df[,c("x1", "y1")] <- t(t(as.matrix(t(t(bp_df[,c("x1", "y1")]) - center[1:2])) %*% rotmat) + center[1:2])

  
# Plotting
png("~/dna_strand.png", width = 400 * diff(bounds) / pi, height = 400 * diff(bounds) / pi)
par(lwd = lwd, xpd = NA, mar = c(0,0,0,0))

plot(NULL, NULL, xlim = c(-diff(bounds)/2, diff(bounds)/2), ylim = bounds, frame = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")

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
  polylines(sdf$x1, sdf$y1, lwd = sdf$lwd1, col = "black")
}

#strand 2
for(i in seq_along(seen_sets_2)){
  seen_set <- seen_sets_2[[i]]
  sdf <- strand_df[seen_set,]
  polylines(sdf$x2, sdf$y2, lwd = sdf$lwd2, col = "black")
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

dev.off()
