#libraries
library(foreach)
library(doParallel)
library(parallel)

#functions
line_lengths <- function(x, segs){
  segn0 <- table(segs)
  nsegs <- length(unique(segs))
  
  seg_dists <- lapply(1:nsegs, function(segi){
    diffs <- t(sapply(which(segs == segi)[-1], function(i) x[i,] - x[i-1,]))
    dists <- apply(diffs, 1, function(i) sqrt(sum(unlist(i)^2)))
  })
  
  sapply(seg_dists, sum)
}

interpolate_to_length <- function(x, n, equal_increments = F){
  if(equal_increments){
    diffs <- diff(x)
    dists <- abs(diffs)
    sumd <- sum(dists)
    cumsumd <- c(0, cumsum(dists))
    at_len <- seq(0, sumd, length.out = n)
    floor_at_len <- sapply(at_len, function(i) sum(i >= cumsumd))
    prop_from_floor <- (at_len - cumsumd[floor_at_len]) / c(dists, 1)[floor_at_len]
    out <- x[floor_at_len] + 
      (x[sapply(floor_at_len + 1, function(i) 
        min(c(i, length(x))))] - x[floor_at_len]) * 
      prop_from_floor
    
    
  } else {
    n_at <- seq(1, length(x), length.out = n)
    n_at_leftover <- n_at %% 1
    n_at_floor <- floor(n_at)
    n_at_ceil <- ceiling(n_at)
    out <- x[n_at_floor] + (x[n_at_ceil] - x[n_at_floor]) * n_at_leftover  
  }
  
  out
}

# x <- c(1:3, 2, 5)
# plot(x, type = "l")
# points(seq(1, length(x), length.out = 30), interpolate_to_length(x, 30, T), col = 2)

break_integer <- function(x, proportions) {
  # Multiply the integer with the proportions and round the result
  chunks <- round(x * proportions)
  
  # Calculate the difference between the original integer and the sum of the rounded chunks
  diff <- x - sum(chunks)
  
  # Find the indices with the largest absolute error due to rounding
  abs_error <- abs(chunks - x * proportions)
  idx <- order(abs_error, decreasing = TRUE)
  
  # Distribute the difference to the elements with the largest absolute error
  if (diff > 0) {
    chunks[idx[1:diff]] <- chunks[idx[1:diff]] + 1
  } else if (diff < 0) {
    chunks[idx[1:(-diff)]] <- chunks[idx[1:(-diff)]] - 1
  }
  
  return(chunks)
}

interpolate_to_length_kd <- function(x, n, equal_increments = F, segs = NA){
  if(equal_increments){
    
    if(is.na(segs[1])){
      segs <- rep(1, nrow(x))
    }
    segn0 <- table(segs)
    nsegs <- length(unique(segs))
    
    seg_dists <- lapply(1:nsegs, function(segi){
      diffs <- t(sapply(which(segs == segi)[-1], function(i) x[i,] - x[i-1,]))
      dists <- apply(diffs, 1, function(i) sqrt(sum(unlist(i)^2)))
    })
    
    seg_cumsumds <- lapply(1:nsegs, function(segi){
        c(0, cumsum(seg_dists[[segi]]))
    })
    
    sumd <- sapply(seg_cumsumds, tail, 1)
    segn <- break_integer(n, sumd / sum(sumd))
    seg_start_index <- cumsum(c(0,segn0[-nsegs]))
      
    out <- lapply(1:nsegs, function(segi){
      at_len <- seq(0, sumd[segi], length.out = segn[segi])
      floor_at_len <- sapply(at_len, function(i) sum(i >= seg_cumsumds[[segi]]))
      prop_from_floor <- (at_len - seg_cumsumds[[segi]][floor_at_len]) / c(seg_dists[[segi]], 1)[floor_at_len]
      x[seg_start_index[segi] + floor_at_len,] + 
        (x[seg_start_index[segi] + sapply(floor_at_len + 1, function(i) min(c(i, segn0[[segi]]))),] - 
           x[seg_start_index[segi] + floor_at_len,]) * prop_from_floor
      
    })
    
    out <- do.call(rbind, out)
    out <- cbind(out, i = rep(1:nsegs, segn))
    
  } else {
  
    out <- apply(x, 2, interpolate_to_length, n = n, equal_increments = equal_increments)
  
  }
  
  out
  
}

#note, only consistent length along original shape, not guaranteed for new shape

# x <- cbind(sample(1:5), sample(1:5))
# interpolate_to_length_kd(x, 10, T)
# plot(x, type = "l")
# points(interpolate_to_length_kd(x, 10, F))

wavelet_smooth_noise_thresh <- function(d, wavelet_filter = "mb4", threshold_multiplier = 3){
  
  # select and apply transformation
  # wavelet_filter <- c("haar", "la8", "d4", "d6", "mb4")[5]
  if(class(d) %in% c("matrix", "data.frame")){
    d <- unlist(d)
  }
  
  n_levels <- max(1, floor(log2(length(d)) / 2))
  waveout <- waveslim::modwt(d, wf = wavelet_filter, n.levels = n_levels)
  
  # Estimate sd of noise and set threshold from noise standard deviation
  noise_sigma <- sd(waveout$d1)
  threshold <- threshold_multiplier * noise_sigma
  
  # Apply soft thresholding to the detail coefficients
  thresholded_coeffs <- waveout
  
  for (i in 1:n_levels) {
    
    level_name <- paste0("d", i)
    
    thresholded_coeffs[[level_name]] <- sign(thresholded_coeffs[[level_name]]) * 
      (abs(thresholded_coeffs[[level_name]]) - threshold) * 
      (abs(thresholded_coeffs[[level_name]]) > threshold)
    
  }
  
  denoised_signal <- waveslim::imodwt(thresholded_coeffs)
  
  # Combine the denoised signals into a 2D time series
  data.frame(denoised_signal)
  
}


wavelet_smooth <- function(d, wavelet_filter = "mb4", var_thresh = 0.98){
  
  # select and apply transformation
  # wavelet_filter <- c("haar", "la8", "d4", "d6", "mb4")[5]
  if(class(d) %in% c("matrix", "data.frame")){
    d <- unlist(d)
  }
  
  n_levels <- max(1, floor(log2(length(d)) / 2))
  waveout <- waveslim::modwt(d, wf = wavelet_filter, n.levels = n_levels)
  
  # Estimate sd of noise and set threshold from noise standard deviation
  squares <- unlist(waveout)^2
  total_ss <- sum(squares)
  sorted_squares <- sort(squares, T)
  threshold <- sqrt(sorted_squares[sum(cumsum(sorted_squares) / total_ss <= max(c(var_thresh, sorted_squares[1] / total_ss * 1.01)))])
  
  # Apply soft thresholding to the detail coefficients
  thresholded_coeffs <- waveout
  
  levels <- c(paste0("d", 1:n_levels), paste0("s", n_levels))
  for (level_name in levels) {
    thresholded_coeffs[[level_name]] <- sign(thresholded_coeffs[[level_name]]) * 
      (abs(thresholded_coeffs[[level_name]]) - threshold) * 
      (abs(thresholded_coeffs[[level_name]]) > threshold)
    
  }
  
  denoised_signal <- waveslim::imodwt(thresholded_coeffs)
  
  # Combine the denoised signals into a 2D time series
  data.frame(denoised_signal)
  
}

logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (1+ exp(x))
dinvlogit <- function(x) exp(x) / (1+ exp(x))^2
ifelse2 <- function(test, yes, no) if(test){return(yes)}else{return(no)}
invlogit_bounds <- function(z, n = 100){
  zvals <- seq(-z, z, length.out = n)
  weights <- seq(1, 0, length.out = floor(n/2))
  weights <- c(weights, ifelse2(n%%2, 0, numeric(0)), rev(weights))
  vals_in_01 <- invlogit(zvals) * (1-weights) + seq(0, 1, length.out = n) * weights
  vals_in_01
}

sigmoid_01 <- function(logc = 2, loc = 0.5, n = 100){
  x <- seq(0, 1, length.out = n)
  conc <- exp(logc) * 2
  shapes <- c(conc * loc, conc * (1-loc))
  weights <- pbeta(x, shapes[1], shapes[2])
  weights
}

wavelet_smooth_interpolate <- function(d1, d2, wavelet_filter = "mb4", threshold_multipliers = c(3,3), props){
  
  # select and apply transformation
  # wavelet_filter <- c("haar", "la8", "d4", "d6", "mb4")[5]
  if(class(d1) %in% c("matrix", "data.frame")){
    d1 <- unlist(d1)
  }
  
  if(class(d2) %in% c("matrix", "data.frame")){
    d2 <- unlist(d2)
  }
  
  #process first signal
  n_levels <- max(1, floor(log2(length(d1)) / 2))
  waveout1 <- waveslim::modwt(d1, wf = wavelet_filter, n.levels = n_levels)
  
  noise_sigma1 <- sd(waveout1$d1)
  threshold1 <- threshold_multipliers[1] * noise_sigma1
  
  thresholded_coeffs1 <- waveout1
  for (i in 1:n_levels) {
    level_name <- paste0("d", i)
    thresholded_coeffs1[[level_name]] <- sign(thresholded_coeffs1[[level_name]]) * 
      (abs(thresholded_coeffs1[[level_name]]) - threshold1) * 
      (abs(thresholded_coeffs1[[level_name]]) > threshold1)
  }
  
  #process second signal
  waveout2 <- waveslim::modwt(d2, wf = wavelet_filter, n.levels = n_levels)
  
  noise_sigma2 <- sd(waveout2$d1)
  threshold2 <- threshold_multipliers[2] * noise_sigma2
  
  thresholded_coeffs2 <- waveout2
  for (i in 1:n_levels) {
    level_name <- paste0("d", i)
    thresholded_coeffs2[[level_name]] <- sign(thresholded_coeffs2[[level_name]]) * 
      (abs(thresholded_coeffs2[[level_name]]) - threshold2) * 
      (abs(thresholded_coeffs2[[level_name]]) > threshold2)
  }
  
  #combine signals
  thresholded_coeffs <- lapply(setNames(props, props), function(prop) thresholded_coeffs2)
  for(prop in props){
    for (i in 1:length(thresholded_coeffs2)) {
      level_name <- i
      thresholded_coeffs[[as.character(prop)]][[level_name]] <- 
        thresholded_coeffs1[[level_name]] * prop + 
        thresholded_coeffs2[[level_name]] * (1 - prop)
    }
  }
  
  # Combine the denoised signals into a time series
  denoised_signal <- lapply(setNames(props, props), function(prop) 
    data.frame(waveslim::imodwt(thresholded_coeffs[[as.character(prop)]])))
  
  
  denoised_signal
  
}


fourier_smooth <- function(d){
  
  if(!(class(d) %in% c("matrix", "data.frame"))){
    d <- data.frame(d)
  }
  
  fftout <- fft(as.matrix(d), F)
  mag <- Mod(fftout)
  
  #magnitude-based filtering
  #can also filter by variance ie w/ Parseval's theorem
  noise_sigma <- median(mag[(nrow(fftout) %/% 2 + 1):nrow(fftout), (ncol(fftout) %/% 2 + 1):ncol(fftout)], na.rm = TRUE) / 0.6745
  threshmag <- noise_sigma * sqrt(2 * log(nrow(d))) #universal threshold
  threshmag <- noise_sigma * sqrt(2 * log(nrow(d))^2 / log(2)) #visushrink threshold
  to_discard <- lapply(1:ncol(d), function(i){
    x <- mag[,i]
    which(x < threshmag)
  })
  to_discard <- Reduce(intersect, to_discard)
  
  #recompose
  fftout[to_discard,] <- 0 + 0i
  filtered_inds <- Re(fft(fftout, T)) / nrow(d) / ncol(d)
  data.frame(filtered_inds)
  
}



linear_smooth_interpolate <- function(d1 = sx1, d2 = sx2, props = props){
  
  lapply(setNames(props, props), function(prop){
    d1 * prop + d2 * (1-prop)
  })
  
}

#maybe also pass it through a more "loop"-y intermediate stage? ie increaase intermediate frequencies?

#for different segment #s, find shortest segments in larger sequence 
#and create fictitious 0-length segments in longer sequence
#at wherever the intersection points end up

#alternatively, just lump short adjacent segments together?

#maybe handle different numbers of line segments by growing and shrinking them?
#and then match the number of line segments in 1-to-1 object sets?
#as in, create fictitious shapes for the extra segments in the shape that is lacking

#can also try to trace an existing signature's path -- 
#collapse to thickness of 1, and then if you reach a junction / cross pick the direction that minimizes
#change in direction

#can also apply a thickening transform -- smoothly increase line thickness propto how much time was spent at that point


#### multi-object matching ####

#functions
# Function to determine if lines intersect
# Each line is defined by two points, row1 = (x1, y1) and (x2, y2)
intersect_func <- function(line1, line2, return_loc = F) {
  
  #find slopes
  s1 <- line1[2,] - line1[1,]
  s1 <- s1$y / s1$x
  s2 <- line2[2,] - line2[1,]
  s2 <- s2$y / s2$x
  
  if(s1 == s2){
    return(F)
  }
  
  if(s2 == Inf){
    return(sum(line1$x > line2$x[1]) == 1)
  } else if(s1 == Inf){
    return(sum(line2$x > line1$x[1]) == 1)
  }
  
  #find intercepts
  i1 <- line1[1,2] - s1 * line1[1,1]
  i2 <- line2[1,2] - s2 * line2[1,1]
  
  #find intersection x_coord  
  x_intersect <- (i1 - i2) / (s2 - s1)
  
  #check it is within bounds of either line
  intersects_TF <- sum(line1$x > x_intersect) == 1
  
  if(!intersects_TF){return(intersects_TF)}
  
  if(return_loc){
    y_intersect <- i1 + s1 * intersects_TF
    return(c(x = x_intersect, y = y_intersect))
  } else {
    return(intersects_TF)  
  }
  
  
  
}

check_overlap <- function(rect1, rect2) {
  if (all(max(rect1$x) < rect2$x) || all(min(rect1$x) > rect2$x)) {
    return(FALSE)
  }
  if (all(max(rect1$y) < rect2$y) || all(min(rect1$y) > rect2$y)) {
    return(FALSE)
  }
  return(TRUE)
}

compress_redundancy_xy <- function(obj){
  obj[!apply(cbind(unlist(lapply(rle(obj$x)$lengths, function(n) rep(c(FALSE, TRUE), times = c(1, n-1)))),
                   unlist(lapply(rle(obj$y)$lengths, function(n) rep(c(FALSE, TRUE), times = c(1, n-1))))), 1, all),]
}


plot_example <- F
if(plot_example){
  
  #### read in data ####
  nv <- read.csv("~/data/nikvetr_sig_ts.csv")
  nl <- read.csv("~/data/nl_sig_ts.csv")
  
  #get block inds
  nv$i <- c(1, cumsum(diff(nv$t) > 50) + 1)
  nl$i <- c(1, cumsum(diff(nl$t) > 50) + 1)
  
  #subset to just the first block for now
  nv <- nv[nv$i == 1,]
  nl <- nl[nl$i == 1,]
  
  #transform to equal lengths
  orig_res <- c(nl = nrow(nl), nv = nrow(nv))
  orig_res_rat <- orig_res / max(orig_res)
  max_len <- max(c(nrow(nl), nrow(nv)))
  nl2 <- interpolate_to_length_kd(x = nl[,c("x", "y")], n = max_len, 
                                  segs = nl$i, equal_increments = T)
  nv2 <- interpolate_to_length_kd(x = nv[,c("x", "y")], n = max_len, 
                                  segs = nv$i, equal_increments = T)
  plot(nl2$x, -nl2$y, col = nl2$i)
  plot(nv2$x, -nv2$y, col = nv2$i)
  
  #center text on origin
  nvw <- diff(range(nv2$x))
  nvh <- diff(range(nv2$y))
  nv2$x <- nv2$x - min(nv2$x) - nvw / 2
  nv2$y <- nv2$y - min(nv2$y) - nvh / 2
  
  nlw <- diff(range(nl2$x))
  nlh <- diff(range(nl2$y))
  nl2$x <- nl2$x - min(nl2$x) - nlw / 2
  nl2$y <- nl2$y - min(nl2$y) - nlh / 2
  
  #get equal length text lines
  sig_lengths <- c(nv = sum(line_lengths(nv2[,c("x", "y")], nv2$i)), 
                   nl = sum(line_lengths(nl2[,c("x", "y")], nl2$i)))
  sig_props <- sig_lengths / max(sig_lengths)
  nv2[,c("x", "y")] <- nv2[,c("x", "y")] * sig_props["nv"]
  nl2[,c("x", "y")] <- nl2[,c("x", "y")] * sig_props["nl"]
  
  #replace w/ new vars
  nv <- nv2
  nl <- nl2
  
  #invert vertical axis
  nv$y <- -nv$y
  nl$y <- -nl$y
  
  plot(nl$x, nl$y, type = "l", col = 0)
  for(i in 1:max(nl$i)){
    lines(nl$x[nl$i == i], nl$y[nl$i == i], type = "l", col = 2, lwd = 2)
  }
  
  plot(nv$x, nv$y, type = "l", col = 0)
  for(i in 1:max(nv$i)){
    lines(nv$x[nv$i == i], nv$y[nv$i == i], type = "l", col = 2, lwd = 2)
  }
  
  #### smoothed signature processing ####
  base_threshold_multiplier <- 50
  prop_var <- 0.999
  use_var_thresh <- T
  nls <- lapply(1:max(nl$i), function(i){
    
    inds <- which(nl$i == i)
    
    sx <- nl$x[inds]
    sx <- c(sx, rev(sx))
    if(use_var_thresh){
      sx <- unlist(wavelet_smooth(sx, var_thresh = prop_var))  
    } else {
      sx <- unlist(wavelet_smooth_noise_thresh(sx, threshold_multiplier = base_threshold_multiplier / orig_res_rat["nl"]^1.5))  
    }
    
    
    sx <- sx[1:length(inds)]
    
    sy <- nl$y[inds]
    sy <- c(sy, rev(sy))
    if(use_var_thresh){
      sy <- unlist(wavelet_smooth(sy, var_thresh = prop_var))
    } else {
      sy <- unlist(wavelet_smooth_noise_thresh(sy, threshold_multiplier = base_threshold_multiplier / orig_res_rat["nl"]^1.5))
    }
    sy <- sy[1:length(inds)]
    
    return(data.frame(x = sx, y = sy, i = i))
  })
  
  nvs <- lapply(1:max(nv$i), function(i){
    
    inds <- which(nv$i == i)
    
    sx <- nv$x[inds]
    sx <- c(sx, rev(sx))
    if(use_var_thresh){
      sx <- unlist(wavelet_smooth(sx, var_thresh = prop_var))
    } else {
      sx <- unlist(wavelet_smooth_noise_thresh(sx, threshold_multiplier = base_threshold_multiplier / orig_res_rat["nv"]^1.5))
    }
    sx <- sx[1:length(inds)]
    
    sy <- nv$y[inds]
    sy <- c(sy, rev(sy))
    if(use_var_thresh){
      sy <- unlist(wavelet_smooth(sy, var_thresh = prop_var))
    } else {
      sy <- unlist(wavelet_smooth_noise_thresh(sy, threshold_multiplier = base_threshold_multiplier / orig_res_rat["nv"]^1.5))
    }
    sy <- sy[1:length(inds)]
    
    return(data.frame(x = sx, y = sy, i = i))
  })
  
  
  #basic plotting
  par(mfrow = c(2,1))
  plot(nl$x, nl$y, type = "l", col = 0)
  for(i in 1:max(nl$i)){
    lines(nls[[i]]$x, nls[[i]]$y, type = "l", col = 1, lwd = 2)
  }
  
  plot(nv$x, nv$y, type = "l", col = 0)
  for(i in 1:max(nv$i)){
    lines(nvs[[i]]$x, nvs[[i]]$y, type = "l", col = 1, lwd = 2)
  }
  
  #### interpolation ####
  
  #wavelet interpolation?
  props <- 0:4/4
  inds <- which(nl$i == 1)
  
  sx1 <- nl$x[inds]
  sx1 <- c(sx1, rev(sx1))
  sx2 <- nv$x[inds]
  sx2 <- c(sx2, rev(sx2))
  
  sxs <- wavelet_smooth_interpolate(d1 = sx1, d2 = sx2, 
                                    threshold_multipliers = c(base_threshold_multiplier, base_threshold_multiplier) / orig_res_rat, 
                                    props = props)
  
  sy1 <- nl$y[inds]
  sy1 <- c(sy1, rev(sy1))
  sy2 <- nv$y[inds]
  sy2 <- c(sy2, rev(sy2))
  
  sys <- wavelet_smooth_interpolate(d1 = sy1, d2 = sy2, 
                                    threshold_multipliers = c(base_threshold_multiplier, base_threshold_multiplier) / orig_res_rat, 
                                    props = props)
  
  #plot it
  i <- 5
  ranges <- rbind(range(sxs[[i]]), range(sys[[i]]))
  plot(1, 1, xlim = ranges[1,], ylim = ranges[2,], type = "l", col = 0, 
       xaxt = "n", yaxt = "n", frame = F, xlab = "", ylab = "")
  lines(unlist(sxs[[i]]), unlist(sys[[i]]), type = "l", col = 1, lwd = 2)
  
  #linear interpolation
  sxs_l <- linear_smooth_interpolate(d1 = wavelet_smooth(sx1, threshold_multiplier = 100), 
                                     d2 = wavelet_smooth(sx2, threshold_multiplier = 100), props = props)
  sys_l <- linear_smooth_interpolate(d1 = wavelet_smooth(sy1, threshold_multiplier = 100), 
                                     d2 = wavelet_smooth(sy2, threshold_multiplier = 100), props = props)
  
  #now write to disk and create animation
  anim_dir <- "~/Pictures/sig_interp/"
  if(!dir.exists(anim_dir)){dir.create(anim_dir)}
  frames_dir <- paste0(anim_dir, "frames/")
  if(!dir.exists(frames_dir)){dir.create(frames_dir)}
  if(length(list.files(frames_dir)) > 0){file.remove(paste0(frames_dir, list.files(frames_dir)))}
  
  fps <- 60
  thin <- 1
  fps <- 60
  ns <- 2.5
  thin <- 1
  frame_indices <- 1:(fps*ns)
  plotting_indices <- seq(1, max(frame_indices), by = thin)
  props <- 1:length(frame_indices) / length(frame_indices)
  
  inds <- which(nl$i == 1)
  
  sxs <- wavelet_smooth_interpolate(d1 = sx1, d2 = sx2, 
                                    threshold_multipliers = c(base_threshold_multiplier, base_threshold_multiplier) / orig_res_rat, 
                                    props = props)
  sys <- wavelet_smooth_interpolate(d1 = sy1, d2 = sy2, 
                                    threshold_multipliers = c(base_threshold_multiplier, base_threshold_multiplier) / orig_res_rat, 
                                    props = props)
  
  sxs_l <- linear_smooth_interpolate(d1 = wavelet_smooth(sx1, threshold_multiplier = 100), 
                                     d2 = wavelet_smooth(sx2, threshold_multiplier = 100), props = props)
  sys_l <- linear_smooth_interpolate(d1 = wavelet_smooth(sy1, threshold_multiplier = 100), 
                                     d2 = wavelet_smooth(sy2, threshold_multiplier = 100), props = props)
  
  sxs_l <- linear_smooth_interpolate(d1 = fourier_smooth(sx1), 
                                     d2 = fourier_smooth(sx2), props = props)
  sys_l <- linear_smooth_interpolate(d1 = fourier_smooth(sy1), 
                                     d2 = fourier_smooth(sy2), props = props)
  
  use_multicore <- T
  if(use_multicore){
    if(exists("cl")){
      stopCluster(cl)
      remove(cl)
    }
    if(!exists("cl")){
      cl <- makeCluster(12, outfile="")
      registerDoParallel(cl)
    }
    getDoParWorkers()
  }
  
  foreach(fi=1:length(plotting_indices), .packages = c("png")) %dopar% {
    print(i)
    i <- plotting_indices[fi]
    
    png(filename = paste0(frames_dir, paste0(rep(0, 5 - nchar(i)), collapse = ""), i, ".png"), 
        width = 800, height = 800, type = "cairo", pointsize = 35)
    
    par(xpd = NA, mar = c(4,4,3,2))
    
    ranges <- rbind(range(sxs[[i]]), range(sys[[i]]))
    
    plot(1, 1, xlim = ranges[1,], ylim = ranges[2,], type = "l", col = 0, 
         xaxt = "n", yaxt = "n", frame = F, xlab = "", ylab = "")
    
    lines(unlist(sxs[[i]]), unlist(sys[[i]]), type = "l", col = 1, lwd = 2)
    
    lines(unlist(sxs_l[[i]]), unlist(sys_l[[i]]), type = "l", col = 2, lwd = 4)
    
    dev.off()
    
  }
  
  
  #stich together animation
  file.remove(paste0(anim_dir, list.files(anim_dir, pattern = "*.mp4")))
  system(paste0("cd ", anim_dir, "; ffmpeg -r ", 
                fps / thin," -f image2 -s 1000x500 -i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p temp.mp4"))
  
  #reverse and append and loop
  system(paste0("cd ", anim_dir, "; ", 
                "ffmpeg -i temp.mp4 -vf reverse rev_temp.mp4; ",
                "touch input.txt;",
                "echo \"file temp.mp4\nfile rev_temp.mp4\" > input.txt;",
                "ffmpeg -f concat -i input.txt -codec copy 1t_temp.mp4; ",
                "ffmpeg -stream_loop 3 -i 1t_temp.mp4 -c copy final.mp4"))
  
  #read in data
  ts1 <- read.csv("~/data/nikvetr_sig_ts.csv")
  ts2 <- read.csv("~/data/nl_sig_ts.csv")
  
  #get block inds
  ts1$i <- c(1, cumsum(diff(ts1$t) > 50) + 1)
  ts2$i <- c(1, cumsum(diff(ts2$t) > 50) + 1)
  
  #compress away redundancy
  ts1 <- compress_redundancy_xy(ts1)
  ts2 <- compress_redundancy_xy(ts2)
  
  #match blocks
  nb1 <- table(ts1$i)
  nb2 <- table(ts2$i)
  
  if(length(nb1) > length(nb2)){
    b2u <- order(nb1, decreasing = T)[1:length(nb2)] #blocks to use
    b2nu <- setdiff(as.numeric(unique(names(nb1))), b2u) #blocks to not use
    matcher <- data.frame(b1 = sort(b2u), b2 = 1:length(nb2)) #match blocks between objs
    b2ui <- ts1[ts1$i %in% b2u,] #find locs of blocks to use in appr obj
    
    intersects <- lapply(b2nu, function(b2n){
      
      #get current block not in use
      b2ni <- ts1[ts1$i == b2n, c("x","y")] 
      
      #evaluate overlap of bounding boxes and then check for line intersects
      block_hits <- lapply(2:nrow(b2ni), function(ri){
        curr_seg_b2nu <- b2ni[ri:(ri-1),] #curr seg
        
        overlapping_rects <- sapply(2:nrow(b2ui), function(ri2){  
          check_overlap(rect1 = curr_seg_b2nu, 
                        rect2 = b2ui[ri2:(ri2-1),c("x", "y")])
        })
        hits <- which(overlapping_rects)
        
        if(length(hits) > 0){
          true_intersects <- sapply(hits, function(hi){
            intersect_func(line1 = curr_seg_b2nu, 
                           line2 = b2ui[hi:(hi+1),c("x", "y")])
          })
          if(any(true_intersects)){
            return(hits[true_intersects][1]) #returns first intersection index in b2ui
          } else {
            return(F)
          }
        } else {
          return(F)
        }
      })
      
      if(any(unlist(block_hits))){
        matched_intersects <- do.call(rbind, lapply(seq_along(block_hits), function(i) 
          if(block_hits[i] != F){return(c(b2ni = i, b2ui = block_hits[i]))}else{return(integer(0))}))
      }
      
      matched_intersects
      
    })
    
    
  } else {
    b2u <- order(nb2, decreasing = T)[1:length(nb1)] #blocks to use
    b2nu <- setdiff(as.numeric(unique(names(nb2))), b2u)
    matcher <- data.frame(b1 = 1:length(nb1), b2 = sort(b2u))
  }
  
  #transform to equal lengths
  orig_res <- c(ts2 = nrow(ts2), ts1 = nrow(ts1))
  orig_res_rat <- orig_res / max(orig_res)
  max_len <- max(c(nrow(ts2), nrow(ts1)))
  ts2.2 <- interpolate_to_length_kd(x = ts2[,c("x", "y")], n = max_len, 
                                    segs = ts2$i, equal_increments = T)
  ts1.2 <- interpolate_to_length_kd(x = ts1[,c("x", "y")], n = max_len, 
                                    segs = ts1$i, equal_increments = T)
  plot(ts2.2$x, -ts2.2$y, col = ts2.2$i)
  plot(ts1.2$x, -ts1.2$y, col = ts1.2$i)
  
  #center text on origin
  ts1w <- diff(range(ts1.2$x))
  ts1h <- diff(range(ts1.2$y))
  ts1.2$x <- ts1.2$x - min(ts1.2$x) - ts1w / 2
  ts1.2$y <- ts1.2$y - min(ts1.2$y) - ts1h / 2
  
  ts2w <- diff(range(ts2.2$x))
  ts2h <- diff(range(ts2.2$y))
  ts2.2$x <- ts2.2$x - min(ts2.2$x) - ts2w / 2
  ts2.2$y <- ts2.2$y - min(ts2.2$y) - ts2h / 2
  
  #get equal length text lines
  sig_lengths <- c(ts1 = sum(line_lengths(ts1.2[,c("x", "y")], ts1.2$i)), 
                   ts2 = sum(line_lengths(ts2.2[,c("x", "y")], ts2.2$i)))
  sig_props <- sig_lengths / max(sig_lengths)
  ts1.2[,c("x", "y")] <- ts1.2[,c("x", "y")] * sig_props["ts1"]
  ts2.2[,c("x", "y")] <- ts2.2[,c("x", "y")] * sig_props["ts2"]
  
  #replace w/ new vars
  ts1 <- ts1.2
  ts2 <- ts2.2
  
  #its1ert vertical axis
  ts1$y <- -ts1$y
  ts2$y <- -ts2$y  
}
