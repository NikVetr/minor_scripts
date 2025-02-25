
#compute smooth shape
n <- 200
n_tot <- 200
t <- seq(0, 2 * pi, length.out = n)
f <- function(t, nridges = 4, rad_range = c(1,1.25)) {
  rad_range[1] + sin(nridges * t) * diff(rad_range)
}
r <- sqrt(f(t, 8, c(1,1.1)))
x <- sin(t) * r
y <- cos(t) * r 
d_smooth <- cbind(x=x,y=y)
plot(d_smooth, type = "l")

#compute spiked shape
regular_polygon <- function(k, r = 1) {
  angles <- seq(0, 2 * pi, length.out = k + 1)[-k - 1]
  x <- r * cos(angles)
  y <- r * sin(angles)
  polygon_coords <- cbind(x = x, y = y)
  return(rbind(polygon_coords, polygon_coords[1,]))
}
d_spiked <- regular_polygon(4, 1.125)
lines(d_spiked)

#interpolate the two
arc_length <- function(coords) {
  # coords is a Nx2 matrix
  diffs <- diff(coords)
  dist_steps <- sqrt(rowSums(diffs^2))
  cumsum(c(0, dist_steps))
}

resample_curve <- function(coords, npoints = 200) {
  # Resample curve coords to have npoints equally spaced along arc-length
  al <- arc_length(coords)
  total_len <- tail(al, 1)
  new_al <- seq(0, total_len, length.out = npoints)
  # interpolate x and y separately
  x_new <- approx(al, coords[,1], xout = new_al)$y
  y_new <- approx(al, coords[,2], xout = new_al)$y
  return(cbind(x = x_new, y = y_new))
}

# Resample both to have the same number of points (e.g., 200)
d_smooth_resampled <- resample_curve(d_smooth, npoints = n_tot)
d_spiked_resampled <- resample_curve(d_spiked, npoints = n_tot)

plot.new()
plot.window(xlim = range(d_smooth_resampled[,1]), ylim = range(d_smooth_resampled[,2]))
text_plot_subsample <- round(seq(1, n_tot, length.out = n_tot/2))
text(d_smooth_resampled[text_plot_subsample,1],
     d_smooth_resampled[text_plot_subsample,2],
     labels = text_plot_subsample)
text(d_spiked_resampled[text_plot_subsample,1],
     d_spiked_resampled[text_plot_subsample,2],
     labels = text_plot_subsample, col = 2, xpd = NA)


#rotate indices to align shapes
find_best_rotation <- function(d1, d2, ind){
  d1_pt <- d1[ind,]
  best_ind_away <- which.min(apply(t((t(d2) - d1_pt))^2, 1, sum)) - ind
  if(sign(best_ind_away) == -1){
    return(nrow(d1) + best_ind_away)
  } else {
    return(best_ind_away)  
  }
}

n_points_to_evaluate <- 100
points_to_evaluate <- round(seq(1, n_tot, length.out = n_points_to_evaluate + 1)[-1])
optimal_displacements <- sapply(points_to_evaluate, function(i) find_best_rotation(d_smooth_resampled, d_spiked_resampled, i))
opt_disp_prop <- optimal_displacements / n_tot
logit <- function(p) log(p / (1-p))
inv_logit <- function(x) exp(x)/(1+exp(x))
disp_dens <- density(logit(opt_disp_prop))
optimal_disp_prop <- disp_dens$x[which.max(disp_dens$y)]
optimal_displacement <- round(inv_logit(optimal_disp) * n_tot)
new_inds <- (1:n_tot + optimal_displacement) %% n_tot + 1
d_smooth_resampled <- d_smooth_resampled[new_inds,]

# Define a common parameter t in [0,1]
tvals <- seq(0,1,length.out=n_tot)

#------------------------------------
# 2. Compute curvature for the spiked shape
#------------------------------------

compute_curvature <- function(coords) {
  # coords: Nx2
  # Use finite differences to approximate first and second derivatives
  # Param: tvals from 0 to 1
  dt <- mean(diff(tvals))  # uniform parameter spacing
  
  x <- coords[,1]
  y <- coords[,2]
  
  # First derivatives
  dx <- diff(x)/dt
  dy <- diff(y)/dt
  
  # Second derivatives
  ddx <- diff(dx)/dt
  ddy <- diff(dy)/dt
  
  # Curvature formula at interior points:
  # κ = |x'y'' - y'x''| / ( (x'^2 + y'^2)^(3/2) )
  # We'll lose a couple points at ends, so pad with NA or handle indexing
  kappa <- rep(NA, length(x))
  for (i in 2:(length(x)-1)) {
    denom <- (dx[i-1]^2 + dy[i-1]^2)^(3/2)
    num <- abs(dx[i-1]*ddy[i-1] - dy[i-1]*ddx[i-1])
    kappa[i] <- num / denom
  }
  # Replace NAs at ends by nearest available values
  kappa[1] <- kappa[2]
  kappa[length(kappa)] <- kappa[length(kappa)-1]
  kappa
}

kappa_spiked <- compute_curvature(d_spiked_resampled)
kappa_smooth <- compute_curvature(d_smooth_resampled)

# Normalize curvature to [0,1]
kappa_min <- min(kappa_spiked, na.rm=TRUE)
kappa_max <- max(kappa_spiked, na.rm=TRUE)
w <- (kappa_spiked - kappa_min) / (kappa_max - kappa_min)  # weight: 0 for flat, 1 for spike

#------------------------------------
# 3. Blend the shapes using curvature-based weights
#------------------------------------
d_final <- cbind(
  x = w * d_spiked_resampled[,"x"] + (1 - w) * d_smooth_resampled[,"x"],
  y = w * d_spiked_resampled[,"y"] + (1 - w) * d_smooth_resampled[,"y"]
)

#------------------------------------
# 4. Plot the final result
#------------------------------------
plot(d_smooth, type="l", asp=1, main="Blended Shape")
lines(d_spiked, col="blue")
lines(d_final, col="red", lwd=2)
legend("topright", legend=c("Smooth","Spiked","Final"), col=c("black","blue","red"), lty=1)