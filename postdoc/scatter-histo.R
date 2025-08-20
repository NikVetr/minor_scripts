## ===================================================================== ##
##   2‑D NORMAL SCATTER WITH MATCHED MARGINAL HISTOGRAMS  (REFINED)      ##
## ===================================================================== ##
## ---------- USER‑ADJUSTABLE SETTINGS --------------------------------- ##
n            <- 1000      # sample size
r            <- -0.8       # correlation  (‑1 < r < 1)
target.ticks <- 10        # desired number of axis ticks (≈)
pad.prop     <- 0.04      # padding around extreme point   (fraction of side)
gap.prop     <- 0.06      # gap between scatter & histograms
hist.prop    <- 0.25      # histogram thickness            (fraction of side)
alpha.pt     <- 0.2       # point transparency
pch.pt       <- 19        # point symbol
equal_axes <- F
plot_red_point <- T
plot_green_point <- T
red_point_approx_q <- c(0.995, 0.95)
draw_blue_circle <- T
blue_circle_thickness_propto_dens <- T

set.seed(123)               # reproducibility
## --------------------------------------------------------------------- ##

### 1 ─ DATA GENERATION -------------------------------------------------
x <- rnorm(n)
y <- r * x + sqrt(1 - r^2) * rnorm(n)      # correlated partner

### 2 ─ CORE GEOMETRY  (robust ‑‑ no point ever outside) -------------
x.range  <- diff(range(x))
y.range  <- diff(range(y))

side.data <- max(x.range, y.range)          # longest data span
pad.len   <- pad.prop * side.data           # absolute padding length
side      <- side.data + 2 * pad.len        # full square side (incl. pad)

half.side <- side / 2                       # handy helper

# centre of the data cloud
x.mid <- mean(range(x))
y.mid <- mean(range(y))

# square limits (scatter frame)
x0 <- x.mid - half.side;  x1 <- x.mid + half.side
y0 <- y.mid - half.side;  y1 <- y.mid + half.side

#exactly equal axes
if(equal_axes){
  x0 <- y0 <- min(c(x0, y0))
  x1 <- y1 <- max(c(x1, y1))
}

# surrounding layout
hist.side <- hist.prop * side
gap       <- gap.prop  * side
xlim      <- c(x0, x1 + gap + hist.side)
ylim      <- c(y0, y1 + gap + hist.side)

### 3 ─ DEVICE PREP -----------------------------------------------------
op <- par(no.readonly = TRUE)              # save user settings
par(mar = c(4, 4, 1, 1), xpd = NA)     # room for ticks/labels
plot.new()
plot.window(xlim, ylim, asp = 1)

### 4 ─ UNIFIED TICK POSITIONS -----------------------------------------
ticks <- ticks_orig <- pretty(c(x, y), n = target.ticks)           # common breaks
ticks <- ticks[ticks >= x0 & ticks <= x1]             # keep inside square
tick.len <- 0.02 * S
lab.step <- max(1, ceiling(length(ticks) / target.ticks * 2))  # label stride

### 5 ─ SCATTER ---------------------------------------------------------
points(x, y, pch = pch.pt, col = rgb(0, 0, 0, alpha.pt))
rect(x0, y0, x1, y1, lwd = 1)

### 6 ─ HISTOGRAMS (BREAKS = TICKS) ------------------------------------
hx <- hist(x, breaks = ticks_orig, plot = FALSE)
hy <- hist(y, breaks = ticks_orig, plot = FALSE)

if(plot_red_point){
  red_point_coords <- c(hx$breaks[sum(cumsum(hx$counts) / sum(hx$counts) < red_point_approx_q[1])],
                        hy$breaks[sum(cumsum(hy$counts) / sum(hy$counts) < red_point_approx_q[2])])
}


## top histogram
y.base <- y1 + gap
for (i in seq_along(hx$counts)) {
  bar.ht <- (hx$counts[i] / max(hx$counts)) * hist.side
  if(bar.ht == 0){next()}
  if(plot_red_point){
    if(red_point_coords[1] <= hx$breaks[i]){
      rect(hx$breaks[i], y.base,
           hx$breaks[i + 1], y.base + bar.ht,
           col = adjustcolor(2, 0.5), border = adjustcolor(2, 0.8))
      next()
    }
  }
  rect(hx$breaks[i], y.base,
       hx$breaks[i + 1], y.base + bar.ht,
       col = "grey80", border = "grey40")
}

## right histogram
x.base <- x1 + gap
for (i in seq_along(hy$counts)) {
  bar.wd <- (hy$counts[i] / max(hy$counts)) * hist.side
  if(bar.wd == 0){next()}
  if(plot_red_point){
    if(red_point_coords[2] <= hy$breaks[i]){
      rect(x.base, hy$breaks[i],
           x.base + bar.wd, hy$breaks[i + 1],
           col = adjustcolor(2, 0.5), border = adjustcolor(2, 0.8))
      next()
    }
  }
  rect(x.base, hy$breaks[i],
       x.base + bar.wd, hy$breaks[i + 1],
       col = "grey80", border = "grey40")
}

#border for histograms
min_breaks <- min(c(hx$breaks, hy$breaks))
segments(x0 = x.base, x1 = x.base, y0 = min_breaks, y1 = y.base, col = 1, lwd = 1.5)
segments(x0 = min_breaks, x1 = x.base, y0 = y.base, y1 = y.base, col = 1, lwd = 1.5)

### 7 ─ AXES, TICKS & LABELS -------------------------------------------
## bottom axis
xticks <- ticks[ticks >= x0 & ticks <= x1]
segments(xticks, y0, xticks, y0 - tick.len)                  # ticks
lab.idx <- seq(1, length(xticks), by = lab.step)            # skip labels if crowded
text(xticks[lab.idx], y0 - tick.len, labels = xticks[lab.idx], pos = 1, col = 1)  # labels
text((x0 + x1) / 2, y0 - tick.len - max(strheight(xticks[lab.idx])) * 2, labels = expression('X'[1]),
     font = 1, pos = 1, cex = 1.5, col = 1)      # axis title

## left axis
yticks <- ticks[ticks >= y0 & ticks <= y1]
segments(x0, y0, x0, y1)
segments(x0, yticks, x0 - tick.len, yticks)
text(x0 - tick.len, yticks[lab.idx], yticks[lab.idx], adj = 1, pos = 2, col = 1)
text(x0 - tick.len - max(strwidth(yticks[lab.idx])) * 1.25, (y0 + y1) / 2, labels = expression('X'[2]), font = 1,
     cex = 1.5, col = 1, pos = 2)

if(plot_green_point){
  green_point_coords <- sum(red_point_coords * c(1, -1)) / sum(c(1, -1) * c(1, -1)) * c(1, -1)
  points(x = green_point_coords[1], y = green_point_coords[2], pch = pch.pt, col = 3)
  aroff_len_prop <- 0.05
  aroff <- sqrt(sum((green_point_coords - red_point_coords)^2)) * aroff_len_prop
  arrows(x0 = green_point_coords[1] + aroff, y0 = green_point_coords[2] + aroff,
         x1 = red_point_coords[1] - aroff, y1 = red_point_coords[2] - aroff, length = 0.1, col = 1)
}

if(draw_blue_circle){
  nb_coords <- 1E2
  gr_dist <- sqrt(sum((green_point_coords - red_point_coords)^2))
  rads <- seq(0, 2 * pi, length.out = nb_coords)
  blue_pts <- t(rbind(sin(rads), cos(rads)) * gr_dist + green_point_coords)
  
  if(blue_circle_thickness_propto_dens){
    
    bp_ldens <- mvtnorm::dmvnorm(blue_pts, mean = c(0,0), sigma = diag(2) + r - diag(2) * r, log = T)
    bp_dens <- exp(bp_dens)
    bp_ldens_norm <- bp_ldens - min(bp_ldens) + 1
    source(file = "~/repos/polylines/R/functions.R")
    polylines(x = blue_pts[,1], y =  blue_pts[,2], lwd = bp_ldens_norm/30, 
              col = 4, border = 4)
  } else {
    lines(blue_pts, col = 4, lwd = 2)  
  }
  
}


#plot red point
if(plot_red_point){
  points(x = red_point_coords[1], y = red_point_coords[2], pch = pch.pt, col = 2)
  #vertical guiding line
  segments(x0 = red_point_coords[1], x1 = red_point_coords[1], 
           y0 = red_point_coords[2], y1 = y.base, 
           col = 2, lty = 3)
  #horizontal guiding line
  segments(x0 = red_point_coords[1], x1 = x.base, 
           y0 = red_point_coords[2], y1 = red_point_coords[2], 
           col = 2, lty = 3)
  #shade in histogram tails
}

### 8 ─ CLEAN‑UP --------------------------------------------------------
par(op)                             # restore original graphical settings
