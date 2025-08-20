# Load necessary libraries
library(grDevices)
library(grid)

#' Plot Idealized Simulated EKG Readout
#'
#' Generates and plots a simulated, idealized Electrocardiogram (EKG/ECG) waveform
#' without noise. Allows customization of the number of cycles, waveform
#' characteristics, and display options like a grid.
#'
#' @param n_cycles Integer. Number of EKG cycles (heartbeats) to simulate. Default: 3.
#' @param duration_cycle Numeric. Duration of a single EKG cycle in seconds.
#'   Determines the heart rate (HR = 60 / duration_cycle). Default: 0.8 (75 bpm).
#' @param sampling_rate Integer. Number of data points per second. Higher values
#'   result in a smoother curve but require more computation. Default: 250.
#' @param start_coords Numeric vector `c(x, y)`. Coordinates for the starting
#'   point of the EKG trace. Default: `c(0, 0)`.
#' @param p_amp Numeric. Amplitude (peak height) of the P wave. Default: 0.15.
#' @param q_amp Numeric. Amplitude (depth) of the Q wave (negative value). Default: -0.1.
#' @param r_amp Numeric. Amplitude (peak height) of the R wave. Default: 1.0.
#' @param s_amp Numeric. Amplitude (depth) of the S wave (negative value). Default: -0.25.
#' @param t_amp Numeric. Amplitude (peak height) of the T wave. Default: 0.3.
#' @param p_dur Numeric. Duration of the P wave in seconds. Default: 0.08.
#' @param qrs_dur Numeric. Duration of the QRS complex in seconds. Default: 0.10.
#' @param t_dur Numeric. Duration of the T wave in seconds. Default: 0.18.
#' @param pr_interval Numeric. Duration from the start of P wave to the start of QRS
#'   complex (includes P wave duration). Default: 0.16.
#' @param st_interval Numeric. Duration from the end of QRS complex to the start of T
#'   wave (the isoelectric segment). Default: 0.12.
#' @param baseline_level Numeric. The baseline voltage level (isoelectric line). Default: 0.
#' @param add_grid Logical. If TRUE, adds a standard EKG grid background. Default: TRUE.
#' @param grid_col_major Character or numeric. Color for major grid lines. Default: "pink".
#' @param grid_col_minor Character or numeric. Color for minor grid lines. Default: "lightpink".
#' @param line_col Character or numeric. Color of the EKG trace. Default: "black".
#' @param lwd Numeric. Line width for the EKG trace. Default: 1.5.
#' @param xlim Numeric vector `c(min, max)`. Optional override for x-axis limits.
#'   If NULL, limits are calculated automatically. Default: NULL.
#' @param ylim Numeric vector `c(min, max)`. Optional override for y-axis limits.
#'   If NULL, limits are calculated automatically. Default: NULL.
#' @param return_data Logical. If TRUE, returns a data frame with time and voltage
#'   invisibly, instead of plotting. Default: FALSE.
#' @param ... Additional arguments passed to the `plot` and `lines` functions.
#'
#' @return If `return_data` is FALSE (default), plots the idealized EKG and returns NULL invisibly.
#'   If `return_data` is TRUE, returns a data frame with columns 'time' and 'voltage' invisibly.
#' @export
#'
#' @examples
#' # Default idealized EKG plot
#' plot_ekg_ideal()
#'
#' # More cycles, faster heart rate (shorter duration)
#' plot_ekg_ideal(n_cycles = 5, duration_cycle = 0.6)
#'
#' # Different color
#' plot_ekg_ideal(line_col = "blue")
#'
#' # No grid, custom amplitudes
#' plot_ekg_ideal(add_grid = FALSE, r_amp = 1.5, t_amp = 0.5, n_cycles = 2)
#'
#' # Get the data instead of plotting
#' ekg_data <- plot_ekg_ideal(return_data = TRUE, n_cycles = 1)
#' head(ekg_data)
#' plot(ekg_data, type='l')
#'
plot_ekg_ideal <- function( # Renamed slightly for clarity, but you can keep the old name
  n_cycles = 3,
  duration_cycle = 0.8, # seconds (implies 60/0.8 = 75 bpm)
  sampling_rate = 250, # points per second (Hz)
  start_coords = c(0, 0), # c(x_start, y_start)
  p_amp = 0.15,
  q_amp = -0.1,
  r_amp = 1.0,
  s_amp = -0.25,
  t_amp = 0.3,
  p_dur = 0.08,
  qrs_dur = 0.10,
  t_dur = 0.18,
  pr_interval = 0.16, # Start of P to start of QRS
  st_interval = 0.12, # End of QRS to start of T
  baseline_level = 0,
  # noise_level parameter removed
  add_grid = F,
  grid_col_major = "pink",
  grid_col_minor = "lightpink",
  line_col = "black",
  lwd = 1.5,
  xlim = NULL,
  ylim = NULL,
  return_data = FALSE,
  ...
) {
  
  # --- Input Validation ---
  stopifnot(
    is.numeric(n_cycles) && n_cycles >= 0 && n_cycles == floor(n_cycles),
    is.numeric(duration_cycle) && duration_cycle > 0,
    is.numeric(sampling_rate) && sampling_rate > 0 && sampling_rate == floor(sampling_rate),
    is.numeric(start_coords) && length(start_coords) == 2,
    is.numeric(p_amp), is.numeric(q_amp), is.numeric(r_amp), is.numeric(s_amp), is.numeric(t_amp),
    is.numeric(p_dur) && p_dur > 0, is.numeric(qrs_dur) && qrs_dur > 0, is.numeric(t_dur) && t_dur > 0,
    is.numeric(pr_interval) && pr_interval >= p_dur,
    is.numeric(st_interval) && st_interval >= 0,
    is.numeric(baseline_level),
    is.logical(add_grid),
    is.logical(return_data),
    is.numeric(lwd) && lwd > 0
  )
  
  # Total number of points per cycle
  n_points_cycle <- floor(duration_cycle * sampling_rate)
  if (n_points_cycle < 10) {
    warning("Low sampling rate for cycle duration, waveform may be poorly defined.")
  }
  if (n_cycles == 0) {
    if (return_data) return(invisible(data.frame(time=numeric(0), voltage=numeric(0))))
    else return(invisible(NULL)) # Nothing to plot
  }
  
  # Time vector for a single cycle
  t_cycle <- seq(0, duration_cycle, length.out = n_points_cycle)
  
  # --- Generate Single Cycle Waveform ---
  y_cycle <- rep(baseline_level, n_points_cycle)
  
  # Helper function to generate a wave segment (approximated by Gaussian)
  generate_wave <- function(time_vec, center_time, duration, amplitude, baseline) {
    # Use standard deviation related to duration for shape control
    sd <- duration / 4 # Adjust this divisor for wider/narrower waves
    if (sd <= 0) return(rep(0, length(time_vec))) # Avoid errors with zero duration
    
    wave <- amplitude * exp(-(time_vec - center_time)^2 / (2 * sd^2))
    return(wave)
  }
  
  # Define timing within the cycle
  p_start <- 0
  p_center <- p_start + p_dur / 2
  qrs_start <- pr_interval
  qrs_mid <- qrs_start + qrs_dur / 2
  q_time <- qrs_start + qrs_dur * 0.15 # Approx timing
  r_time <- qrs_start + qrs_dur * 0.4  # Approx timing
  s_time <- qrs_start + qrs_dur * 0.7  # Approx timing
  qrs_end <- qrs_start + qrs_dur
  t_start <- qrs_end + st_interval
  t_center <- t_start + t_dur / 2
  t_end <- t_start + t_dur
  
  # Check if T wave fits within cycle duration
  if (t_end > duration_cycle) {
    warning("Calculated T wave end exceeds cycle duration. Adjusting T duration or intervals might be needed.")
    t_dur <- max(0, duration_cycle - t_start)
    t_center <- t_start + t_dur / 2
  }
  
  # Generate waves and add to baseline
  y_cycle <- y_cycle + generate_wave(t_cycle, p_center, p_dur, p_amp, baseline_level)
  
  # QRS Complex
  qrs_sd_factor <- 6 # Make QRS sharper than P/T waves
  y_cycle <- y_cycle + generate_wave(t_cycle, q_time, qrs_dur / 3, q_amp, baseline_level) # Q
  y_cycle <- y_cycle + generate_wave(t_cycle, r_time, qrs_dur / 2.5, r_amp, baseline_level) # R
  y_cycle <- y_cycle + generate_wave(t_cycle, s_time, qrs_dur / 3, s_amp, baseline_level) # S
  
  # T Wave
  y_cycle <- y_cycle + generate_wave(t_cycle, t_center, t_dur, t_amp, baseline_level)
  
  # --- Combine Cycles ---
  y_full <- rep(y_cycle, n_cycles)
  t_full <- seq(0, duration_cycle * n_cycles, length.out = n_points_cycle * n_cycles)
  
  # --- Noise Addition Section Removed ---
  # if (noise_level > 0) {
  #   noise <- rnorm(length(y_full), mean = 0, sd = noise_level)
  #   y_full <- y_full + noise
  # }
  
  # --- Adjust for Start Coordinates ---
  t_full <- t_full + start_coords[1]
  y_full <- y_full + start_coords[2]
  
  # --- Prepare Data Frame ---
  ekg_data <- data.frame(time = t_full, voltage = y_full)
  
  # --- Return Data if requested ---
  if (return_data) {
    return(invisible(ekg_data))
  }
  
  # --- Plotting ---
  
  # Determine plot limits if not provided
  if (is.null(xlim)) {
    xlim <- range(ekg_data$time)
    # Add a little padding
    xlim <- xlim + c(-0.05, 0.05) * diff(xlim)
    if(diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5) # Handle single point case
  }
  if (is.null(ylim)) {
    ylim <- range(ekg_data$voltage)
    # Add padding, ensuring reasonable vertical space
    padding <- max(0.2, 0.1 * diff(ylim))
    ylim <- ylim + c(-padding, padding)
    if(diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5) # Handle flat line case
  }
  
  # Set up plot area
  plot(NULL, type = "n", xlim = xlim, ylim = ylim,
       xlab = "Time (s)", ylab = "Voltage (mV)",
       main = "Simulated Idealized EKG Readout", ...) # Updated title
  
  # --- Add EKG Grid ---
  if (add_grid) {
    x_range <- diff(xlim)
    y_range <- diff(ylim)
    x_step_minor <- ifelse(x_range < 0.4, 0.1, 0.04)
    x_step_major <- ifelse(x_range < 2.0, 0.5, 0.2)
    y_step_minor <- ifelse(y_range < 1.0, 0.2, 0.1)
    y_step_major <- ifelse(y_range < 2.5, 1.0, 0.5)
    
    abline(v = seq(floor(xlim[1] / x_step_minor) * x_step_minor, xlim[2], by = x_step_minor),
           col = grid_col_minor, lty = "dotted")
    abline(h = seq(floor(ylim[1] / y_step_minor) * y_step_minor, ylim[2], by = y_step_minor),
           col = grid_col_minor, lty = "dotted")
    abline(v = seq(floor(xlim[1] / x_step_major) * x_step_major, xlim[2], by = x_step_major),
           col = grid_col_major, lty = "solid", lwd=1.2)
    abline(h = seq(floor(ylim[1] / y_step_major) * y_step_major, ylim[2], by = y_step_major),
           col = grid_col_major, lty = "solid", lwd=1.2)
  }
  
  # --- Plot EKG Trace ---
  lines(ekg_data$time, ekg_data$voltage, col = line_col, lwd = lwd, ...)
  
  invisible(NULL) # Return NULL invisibly when plotting
}

# --- Examples ---

#get ecg coords and info
n_cycles <- 1
r_amp <- 0.7
ecg_coords <- plot_ekg_ideal(n_cycles = n_cycles, duration_cycle = 1, 
                             r_amp = r_amp, return_data = T, sampling_rate = 800)
xr_ecg <- range(ecg_coords$time)
w_ecg <- diff(xr_ecg)
yr_ecg <- range(ecg_coords$voltage)
h_ecg <- diff(yr_ecg)

#get name coords and info

#get connecting tail from "a"
string <- "a lavender"
source("~/scripts/minor_scripts/postdoc/hershey_vec.R")
a_vec <- orig_smooth_vec[[1]]

#process rest of name
name_vec <- orig_smooth_vec[-1]
cat_name <- data.frame(do.call(rbind, name_vec))
nseg <- length(name_vec)
name_inds <- rep(1:nseg, sapply(name_vec, nrow))
xr_name <- range(cat_name[,1])
w_name <- diff(xr_name)
yr_name <- range(cat_name[,2])
h_name <- diff(yr_name)
name_scale <- h_ecg / h_name / 1.5
cat_name <- cat_name * name_scale

#add tail to ecg 
a_vec <- a_vec * name_scale
bottom_of_name <- min(a_vec[,2])
height_of_connection <- tail(a_vec[,2], 1)
a_tail_poss_inds <- floor(seq(0.75 * nrow(a_vec), nrow(a_vec), by = 1))
tail_start <- which.min(a_vec[a_tail_poss_inds,2] - bottom_of_name) + a_tail_poss_inds[1] - 1
a_tail <- a_vec[tail_start:nrow(a_vec),]
a_tail_disp <- c(x = tail(ecg_coords[,1], 1) -  a_tail[1,1],
                 y = tail(ecg_coords[,2], 1) -  a_tail[1,2])
tailed_ecg_coords <- rbind(as.matrix(ecg_coords), cbind(a_tail[-1,1] + a_tail_disp[1], 
                                                        a_tail[-1,2] + a_tail_disp[2]))

#connect name to tail
# name_disp <- a_tail_disp
# name_disp[1] <- name_disp[1] + tail(a_vec$x, 1) - cat_name$x[1]
name_disp <- as.numeric(tail(tailed_ecg_coords, 1) - head(cat_name, 1))
disp_cat_name <- data.frame(x = cat_name$x + name_disp[1], 
                  y = cat_name$y + name_disp[2])
name_vec <- split(disp_cat_name, name_inds)

#translate these so they all touch
# ──────────────────────────────────────────────────────────────────────────────
# Align every path to the previous one by a pure X‑shift
# (works on the list structure you showed: each element is a data.frame x / y)
# ──────────────────────────────────────────────────────────────────────────────
align_vecs_horiz <- function(vecs) {
  stopifnot(length(vecs) >= 2)
  
  # helper: find the closest point in B to every point in A, get global minimum
  closest_shift <- function(A, B) {
    # matrices: rows = points, cols = x,y   (RANN likes matrices)
    Am <- as.matrix(A)
    Bm <- as.matrix(B)
    
    # fast nearest‑neighbour search
    res <- RANN::nn2(data = Bm, query = Am, k = 1)
    
    k     <- which.min(res$nn.dists)       # index of globally closest pair
    idx_A <- k                             # point in A
    idx_B <- res$nn.idx[k, 1]              # its mate in B
    
    # shift needed so the X‑coords touch
    shift_x <- A$x[idx_A] - B$x[idx_B]
    
    list(shift = shift_x,
         idx_A = idx_A,
         idx_B = idx_B,
         dist  = res$nn.dists[k])
  }
  
  out <- vecs            # copy to modify
  for (i in 2:length(vecs)) {
    sh <- closest_shift(out[[i - 1]], out[[i]])
    out[[i]]$x <- out[[i]]$x + sh$shift     # translate current vector
  }
  out
}

vecs <- c(list(data.frame(x = tailed_ecg_coords[,1], y = tailed_ecg_coords[,2])), name_vec)
# vecs <- align_vecs_horiz(vecs)

#plot output
source("~/repos/polylines/R/functions.R")
lwd <- 1.5
plot.new()
plot.window(xlim = range(do.call(rbind, vecs)[,1]), 
            ylim = range(do.call(rbind, vecs)[,2]), asp = 1)
for(i in 1:length(vecs)){
  lines(x = vecs[[i]]$x, 
        y = vecs[[i]]$y, col = 1, lwd = lwd)
}


#### try polylines version ####

#' Dynamic thickness for pen‑style strokes
#' x, y        : numeric vectors of equal length (≥3)
#' w_dir       : weight for downward‑stroke cue
#' w_curv      : weight for curvature cue
#' w_speed     : weight for speed cue
#' tmin, tmax  : output range
calc_thickness <- function(x, y,
                           w_dir   = 0.5,
                           w_curv  = 0.3,
                           w_speed = 0.2,
                           tmin    = 0.5,
                           tmax    = 2.0,
                           eps = 1e-6) {
  
  n <- length(x)
  if (length(y) != n || n < 3)
    stop("x and y must be the same length (≥3).")
  
  ## segment vectors ----------------------------------------------------------
  dx <- diff(x)
  dy <- diff(y)
  seg_len <- sqrt(dx^2 + dy^2) + eps            # speed proxy
  
  ## 1. Direction component (dot with downward unit vector) -------------------
  dir_comp <- pmax(0, (-dy) / seg_len)          # 0 (not down) … 1 (pure down)
  dir_comp <- c(dir_comp, tail(dir_comp, 1))    # pad to length n
  
  ## 2. Curvature component (angle between consecutive segments) --------------
  v1x <- dx[-length(dx)]; v1y <- dy[-length(dy)]
  v2x <- dx[-1];           v2y <- dy[-1]
  cosang <- (v1x*v2x + v1y*v2y) /
    (sqrt(v1x^2+v1y^2) * sqrt(v2x^2+v2y^2) + eps)
  angle   <- acos(pmin(1, pmax(-1, cosang)))    # 0 … π
  curv_comp <- c(0, angle / pi, 0)              # normalize, pad
  
  ## 3. Speed component (inverse length, normalised) --------------------------
  speed_comp <- 1 / seg_len
  speed_comp <- speed_comp / max(speed_comp)
  speed_comp <- c(speed_comp, tail(speed_comp, 1))
  
  ## weighted blend -----------------------------------------------------------
  wtot <- w_dir + w_curv + w_speed
  score <- (w_dir   * dir_comp +
              w_curv  * curv_comp +
              w_speed * speed_comp) / wtot
  
  ## rescale to [tmin, tmax] --------------------------------------------------
  tmin + (tmax - tmin) * score
}

redistribute <- function(d, n, oversample_factor = 10) {
  # Ensure input is matrix
  d <- as.matrix(d)
  stopifnot(ncol(d) == 2)
  
  # Oversample between each pair of points
  interpolate_segment <- function(p1, p2, m) {
    t <- seq(0, 1, length.out = m + 1)[-1]  # remove first point to avoid duplication
    outer(1 - t, p1) + outer(t, p2)
  }
  
  segments <- lapply(1:(nrow(d) - 1), function(i) {
    interpolate_segment(d[i, ], d[i + 1, ], oversample_factor)
  })
  path <- rbind(d[1, ], do.call(rbind, segments))
  
  # Compute cumulative distance
  deltas <- diff(path)
  dists <- sqrt(rowSums(deltas^2))
  cumdist <- c(0, cumsum(dists))
  
  # Interpolate n equally spaced points along the path
  target_dists <- seq(0, max(cumdist), length.out = n)
  x_interp <- approx(cumdist, path[,1], xout = target_dists)$y
  y_interp <- approx(cumdist, path[,2], xout = target_dists)$y
  
  data.frame(x = x_interp, y = y_interp)
}

blending_simplex_path <- function(k, n, shape = 2, pow_ends = 2) {
  a <- 10^shape
  t_grid <- seq(0, 1, length.out = n)
  t_breaks <- seq(0, 1, length.out = k + 1)
  sigmoid <- function(x) 1 / (1 + exp(-x))
  weights <- matrix(0, nrow = n, ncol = k)
  for (i in seq_len(k)) {
    lower <- sigmoid(a * (t_grid - t_breaks[i]))
    upper <- sigmoid(a * (t_grid - t_breaks[i + 1]))
    weights[, i] <- lower - upper
  }
  weights_pow <- abs(t_grid * 2 - 1)^pow_ends * (max(0, pow_ends-1)) + 1
  weights <- weights^weights_pow
  row_sums <- rowSums(weights)
  weights <- weights / row_sums
  return(weights)
}

smooth_exp <- function(d, rate = 0.5, buffer = NULL) {
  d <- as.matrix(d)
  stopifnot(ncol(d) == 2)
  n <- nrow(d)
  
  if (is.null(buffer)) {
    # reasonable default window size
    buffer <- ceiling(10 / rate)
  }
  
  # Pad the data to handle edges (symmetric extension)
  d_padded <- rbind(
    d[rep(1, buffer), ],
    d,
    d[rep(n, buffer), ]
  )
  
  # Precompute weights
  k <- 0:buffer
  weights <- dexp(k, rate)
  weights_full <- c(rev(weights[-1]), weights)
  weights_full <- weights_full / sum(weights_full)
  
  # Smoothed output
  smoothed <- t(sapply((buffer + 1):(buffer + n), function(i) {
    window <- d_padded[(i - buffer):(i + buffer), ]
    colSums(window * weights_full)
  }))
  
  colnames(smoothed) <- colnames(d)
  data.frame(smoothed)
}

all_coords <- list()
tail_pow <- 2
tail_prop_floor <- 0.5
base_lwd <- 0.01
lwd_pow <- 2
plot.new()
plot.window(xlim = range(do.call(rbind, vecs)[,1]), 
            ylim = range(do.call(rbind, vecs)[,2]), asp = 1)
for(i in 1:length(vecs)){
  coords <- vecs[[i]]
  coords <- redistribute(coords, nrow(coords))
  coords <- smooth_exp(coords)
  ncoords <- nrow(coords)
  tail_weights <- c(1 - blending_simplex_path(2, floor(ncoords/2), 1, pow_ends = tail_pow)[,1],
                    blending_simplex_path(2, ceiling(ncoords/2), 1, pow_ends = tail_pow)[,1]) * (1-tail_prop_floor) + tail_prop_floor
  lwd <- rep(base_lwd, ncoords) * calc_thickness(x = coords$x, y = coords$y) ^ lwd_pow * tail_weights
  all_coords[[length(all_coords) + 1]] <- data.frame(coords, lwd = lwd)
}

#now transform all_coords into polygons
polys <- list()
for(i in 1:length(all_coords)){
  polys[[i]] <- polylines(all_coords[[i]]$x, y = all_coords[[i]]$y, 
                          lwd = all_coords[[i]]$lwd, complex = T, return_info = T)
}

#now merge these
fuse_polys <- function(polys, buffer_dist = 0.01) {
  # Convert list of data.frames to sf POLYGONs
  sf_polys <- lapply(polys, function(p) {
    if (!all(p[1, ] == p[nrow(p), ])) {
      p <- rbind(p, p[1, ])  # Ensure closed
    }
    sf::st_polygon(list(as.matrix(p)))
  })
  
  # Combine into one sf object
  sf_union_input <- sf::st_union(do.call(c, lapply(sf_polys, sf::st_buffer, dist = buffer_dist)))
  
  # Contract back by the same distance
  result <- sf::st_buffer(sf_union_input, dist = -buffer_dist*9/10)
  
  return(result)
}

fused_poly <- fuse_polys(polys, buffer_dist = 0.0)
plot(fused_poly, col = 1, border = 1)


#try to add first name in too?
string <- "Katya"
source("~/scripts/minor_scripts/postdoc/hershey_vec.R")


#####
name_vec <- orig_smooth_vec

#get rid of tail from last letter
nname <- length(name_vec)
bottom_end_ind <- max(order(abs(name_vec[[nname]]$y - min(name_vec[[nname]]$y)))[1:5])
name_vec[[nname]] <- name_vec[[nname]][1:bottom_end_ind,]

#do the earlier processing
cat_name <- data.frame(do.call(rbind, name_vec))
nseg <- length(name_vec)
name_inds <- rep(1:nseg, sapply(name_vec, nrow))
xr_name <- range(cat_name[,1])
w_name <- diff(xr_name)
yr_name <- range(cat_name[,2])
h_name <- diff(yr_name)
name_scale <- h_ecg / h_name / 1.5
cat_name <- cat_name * name_scale

#extend lhs of ecg
ecg_coords_lext <- data.frame(rbind(cbind(seq(-0.5, 0, by = diff(tail(ecg_coords$time, 2))), ecg_coords$voltage[1]), 
                         as.matrix(head(ecg_coords, nrow(ecg_coords) * 0.9))))
right_ext <- apply(do.call(rbind, lapply(1:200, function(x) apply(tail(ecg_coords_lext, 2), 2, diff))), 2, cumsum)
right_ext[,1] <- right_ext[,1] + as.numeric(tail(ecg_coords_lext, 1)[1])
right_ext[,2] <- right_ext[,2] + as.numeric(tail(ecg_coords_lext, 1)[2])
ecg_coords_lext <- rbind(ecg_coords_lext, right_ext)
ecg_coords_lext <- ecg_coords_lext * 0.8
plot(ecg_coords_lext, type = "l")

#combine with ecg
name_disp <- as.numeric(head(ecg_coords_lext, 1) - tail(cat_name, 1))
disp_cat_name <- data.frame(x = cat_name$x + name_disp[1], 
                            y = cat_name$y + name_disp[2])
name_vec <- split(disp_cat_name, name_inds)

vecs <- c(name_vec, list(data.frame(x = ecg_coords_lext[,1], y = ecg_coords_lext[,2])))
vecs <- c(head(name_vec, length(name_vec) - 1), 
          list(rbind(tail(name_vec, 1)[[1]], 
                     data.frame(x = ecg_coords_lext[,1], y = ecg_coords_lext[,2]))))

#get the polygons
all_coords <- list()
tail_pow <- 2
tail_prop_floor <- 0.5
base_lwd <- 0.0075
lwd_pow <- 2
plot.new()
plot.window(xlim = range(do.call(rbind, vecs)[,1]), 
            ylim = range(do.call(rbind, vecs)[,2]), asp = 1)
for(i in 1:length(vecs)){
  coords <- vecs[[i]]
  coords <- redistribute(coords, nrow(coords))
  coords <- smooth_exp(coords)
  ncoords <- nrow(coords)
  tail_weights <- c(1 - blending_simplex_path(2, floor(ncoords/2), 1, pow_ends = tail_pow)[,1],
                    blending_simplex_path(2, ceiling(ncoords/2), 1, pow_ends = tail_pow)[,1]) * (1-tail_prop_floor) + tail_prop_floor
  lwd <- rep(base_lwd, ncoords) * calc_thickness(x = coords$x, y = coords$y) ^ lwd_pow * tail_weights
  all_coords[[length(all_coords) + 1]] <- data.frame(coords, lwd = lwd)
}
polys <- list()
for(i in 1:length(all_coords)){
  polys[[i]] <- polylines(all_coords[[i]]$x, y = all_coords[[i]]$y, 
                          lwd = all_coords[[i]]$lwd, complex = T, return_info = T)
}

fused_poly_fname <- fuse_polys(polys, buffer_dist = 0.0)
plot(fused_poly_fname, col = 1, border = 1)

#### now plot them together ###
plot(fused_poly, col = 1, border = 1, ylim = c(0.25,1), xpd = NA)
plot(fused_poly_fname * 1.85 + c(2.65,0.75), col = 1, border = 1, add = T, xpd = NA)
