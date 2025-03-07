
log_posterior <- function(theta, dat) {
  
  # Component densities
  lp_banana <- dnorm(theta[,1], dat$banana_mu[1], dat$banana_sd[1], log=TRUE) +
    dnorm(theta[,2], dat$banana_mu[2] + theta[,1]^2, dat$banana_sd[2], log=TRUE)
  
  lp_blob1 <- dnorm(theta[,1], dat$blob_1_mu[1], dat$blob_1_sd[1], log=TRUE) +
    dnorm(theta[,2], dat$blob_1_mu[2], dat$blob_1_sd[2], log=TRUE)
  
  lp_blob2 <- dnorm(theta[,1], dat$blob_2_mu[1], dat$blob_2_sd[1], log=TRUE) +
    dnorm(theta[,2], dat$blob_2_mu[2], dat$blob_2_sd[2], log=TRUE)
  
  # logsumexp over these
  matrixStats::rowLogSumExps(cbind(log(dat$probs[1]) + lp_banana, 
                                   log(dat$probs[2]) + lp_blob1, 
                                   log(dat$probs[3]) + lp_blob2))
}

birth_component <- function(old_fit, data, new_weight = 0.05, up = T) {
  
  old_lambda  <- old_fit$lambda
  old_mu      <- old_fit$mu
  old_sigma   <- old_fit$sigma
  k           <- length(old_lambda)
  
  if(!up){
    best_components <- order(old_fit$lambda, decreasing = T)[1:2]
    new_lambda_vec <- old_lambda[best_components] / sum(old_lambda[best_components])
    new_mu_list    <- old_mu[best_components]
    new_sigma_list <- old_sigma[best_components]
    
    out <- list(lambda = new_lambda_vec,
                mu     = new_mu_list,
                sigma  = new_sigma_list)
    
  } else {
    
    # 1) Compute mixture PDF at each point under old_fit:
    pdf_vals <- numeric(nrow(data))
    for (j in seq_len(k)) {
      pdf_vals <- pdf_vals + old_lambda[j] *
        dmvnorm(data, mu = old_mu[[j]], sigma = old_sigma[[j]])
    }
    
    # 2) "Inverse-density" weighting => highlight underfit regions
    inv_weights <- 1 / (pdf_vals + 1e-12)
    inv_weights <- inv_weights / sum(inv_weights)
    
    # 3) Weighted average of data => new mean
    new_mu <- colSums(data * inv_weights)
    
    # 4) Simple covariance guess for the new component (identity)
    new_sigma <- diag(1, 2)
    
    # 5) Adjust old weights to sum to (1 - new_weight)
    lam_sum <- sum(old_lambda)
    old_lambda_scaled <- old_lambda * ((1 - new_weight) / lam_sum)
    new_lambda_vec <- c(old_lambda_scaled, new_weight)
    
    # Combine old means/sigmas with the new component
    new_mu_list    <- c(old_mu,    list(as.numeric(new_mu)))
    print(names(new_mu_list))
    new_sigma_list <- c(old_sigma, list(new_sigma))
    
    out <- list(lambda = new_lambda_vec,
                mu     = new_mu_list,
                sigma  = new_sigma_list)
    
  }
  
  return(out)
}


sample_from_mix <- function(n, gm_fit) {
  # gm_fit: result of mvnormalmixEM
  # randomly assign each new draw to a component j with prob=gm_fit$lambda[j]
  
  k  <- length(gm_fit$lambda)
  z  <- sample.int(k, size = n, replace = TRUE, prob = gm_fit$lambda)
  
  out <- matrix(NA_real_, n, 2)
  idx_current <- 1L
  for (j in seq_len(k)) {
    nj <- sum(z == j)
    if (nj > 0) {
      out[z == j, ] <- mvrnorm(nj, gm_fit$mu[[j]], gm_fit$sigma[[j]])
    }
  }
  colnames(out) <- c("x","y")
  out
}

plot_density_grid <- function(xy_sub, cell_dim, cols, xlim = NULL, ylim = NULL){
  plot.new()
  if(is.null(xlim)){
    xlim <- range(xy_sub$x)
  }
  if(is.null(ylim)){
    ylim <- range(xy_sub$y)
  }
  plot.window(xlim = xlim, ylim = ylim)
  axis(1); axis(2)
  mtext("x", 1, line = 2.5, cex = 2, font = 2)
  mtext("y", 2, line = 2.5, cex = 2, font = 2, )
  for(i in 1:length(cols)){
    rect(xleft = xy_sub$x[i] - cell_dim$w/2,
         ybottom = xy_sub$y[i] - cell_dim$h/2,
         xright = xy_sub$x[i] + cell_dim$w/2,
         ytop = xy_sub$y[i] + cell_dim$h/2,
         col = cols[i], border = NA
    )
  }  
}

#fit mixture model
compute_bic <- function(mixfit, data) {
  N    <- nrow(data)
  k    <- length(mixfit$lambda)   # number of components
  logL <- mixfit$loglik          # final log-likelihood from EM
  p    <- (k - 1) + 2*k + 3*k     # = 6k - 1
  bic  <- -2*logL + p*log(N)
  bic
}

log_mix_density <- function(theta, mixfit) {
  pi_vec    <- mixfit$lambda
  mu_list   <- mixfit$mu
  sigma_list<- mixfit$sigma
  k <- length(pi_vec)
  
  # We'll compute for each row i: sum_j pi_j * exp( logN(theta_i|mu_j,Sigma_j) )
  # Then take log of that sum.
  n <- nrow(theta)
  comp_vals <- matrix(NA, n, k)
  for (j in seq_len(k)) {
    comp_vals[,j] <- pi_vec[j] * dmvnorm(theta, mu=mu_list[[j]], sigma=sigma_list[[j]])
  }
  mix_dens <- rowSums(comp_vals)
  log(mix_dens)
}

density_grid <- function(n_grid = 100, xr = c(-5,5), yr = c(-10,10),
                         xy_grid = NULL, cell_dim = NULL, z_raw = NULL,
                         log_posterior = NULL, use_mass = F, theta = NULL, 
                         dat = NULL) {
  
  if(is.null(xy_grid)){
    xy_grid <- expand.grid(x = seq(xr[1],xr[2],length.out=n_grid+1), 
                           y = seq(yr[1],yr[2],length.out=n_grid+1))
  }
  if(is.null(cell_dim)){
    cell_dim <- list(w = diff(xr) / (n_grid), 
                     h = diff(yr) / (n_grid))  
  }
  if(is.null(z_raw)){
    z_raw <- z <- log_posterior(xy_grid, dat)  
  } else {
    z <- z_raw
  }
  
  #normalize z_raw to total cell count (to ensure total prob = 1)
  z_raw <- z_raw - log(sum(exp(z_raw)))
  if(!use_mass){
    #transform to an average density otherwise
    z_raw <- z_raw - log(prod(unlist(cell_dim)))
  }
  
  
  #transform these values to a color-space
  
  #for easier visual distinction, set cell values below the mean to 0
  z <- z - log(mean(exp(z)))
  z[z<0] <- -Inf
  
  #transform to colors
  z <- exp(z)
  z <- z / max(z) * 100
  xy_sub <- xy_grid[z>0,]
  z_sub <- z[z>0]
  zraw_sub <- z_raw[z>0]
  colpal <- colorRampPalette(c("white", "black"))(100)
  sqrtz <- sqrt(z_sub)
  sqrtz <- sqrtz - min(sqrtz)
  sqrtz <- sqrtz / max(sqrtz) * 99 + 1
  zcoli <- round(sqrtz)
  cols <- colpal[zcoli]
  return(list(
    xy_sub = xy_sub,
    cell_dim = cell_dim,
    cols = cols,
    colpal = colpal,
    zraw_sub = zraw_sub,
    sqrtz = sqrtz,
    xy_grid = xy_grid
  ))
}

compute_log_hist <- function(theta_final, n_grid, xr, yr, use_mass = F) {
  # 1) Some setup
  n <- nrow(theta_final)
  counts <- numeric(n_grid^2)  # 1D array of counts for each cell
  
  # 2) Grid spacing for x, y
  #    We'll treat [xr[1], xr[2]) as subdivided into n_grid uniform bins.
  #    So each bin is step_x wide in x, step_y wide in y.
  step_x <- (xr[2] - xr[1]) / n_grid
  step_y <- (yr[2] - yr[1]) / n_grid
  
  # 3) Bin each point in O(n) time
  x_data <- theta_final[,1]
  y_data <- theta_final[,2]
  
  # for i=1..n
  for (i in seq_len(n)) {
    xi <- x_data[i]
    yi <- y_data[i]
    
    # 3a) Compute bin indices along x, y
    #     i_x in [0, n_grid-1] if within range
    i_x <- floor((xi - xr[1]) / step_x)
    i_y <- floor((yi - yr[1]) / step_y)
    
    # 3b) Check boundary
    if (i_x < 0 || i_x >= n_grid) next
    if (i_y < 0 || i_y >= n_grid) next
    
    # 3c) Flatten 2D -> 1D index
    # by convention, let bin_id = i_y*n_grid + i_x + 1
    bin_id <- i_y * n_grid + i_x + 1  # +1 because R is 1-based
    counts[bin_id] <- counts[bin_id] + 1
  }
  
  # 4) Convert to log-frequency
  total_count <- sum(counts)
  if (total_count == 0) {
    # if all points were out of [xr, yr], everything is -Inf
    return(rep(-Inf, n_grid^2))
  }
  freq <- counts / total_count
  
  if(use_mass){
    log_freq <- log(freq)  # zero-count cells => -Inf
    return(log_freq)  
  } else {
    freq / (step_x * step_x)
  }
  
}

prettyScientificExpr <- function(numbers, digits = 2, pretty = T, scientific = T, exp_scientific = 4) {
  if(length(numbers) > 1){
    
    if(pretty) numbers <- pretty(numbers, length(numbers))
    
    if(scientific){
      out <- lapply(numbers, prettyScientificExpr, digits = digits)
      return(list(values = sapply(out, function(x) x[["values"]]), 
                  labels = lapply(out, function(x) x[["labels"]])))
      return()
    } else {
      return(list(values = numbers, 
                  labels = as.list(numbers)))
    }
    
  } else {
    
    split_number <- strsplit(formatC(numbers, format = "e", digits = digits), "e")[[1]]
    coefficient <- split_number[1]
    exponent <- as.numeric(split_number[2])
    
    
    
    if (exponent == 0) {
      result <- parse(text = coefficient)
    } else if(exponent >= exp_scientific || exponent < 0){
      result <- parse(text = paste0("list(", coefficient, " %*% 10^", exponent, ")"))
    } else {
      result <- parse(text = as.character(numbers))
    }
    
    return(list(values = numbers, 
                labels = result))
    
  }
  
}

add_continuous_legend <- function(colors, labels, positions, n_rect = 100, x = NA, y = NA, h = NA, w = NA, 
                                  vertical = TRUE, xpd = NA, n_labels = 5, 
                                  left_below = TRUE, log_scale = FALSE,
                                  n_digits = 3, pretty = T, scientific = T, main = NA) {
  
  usr <- par("usr")
  xyr <- diff(range(usr[1:2])) / par("pin")[1] / diff(range(usr[3:4])) * par("pin")[2]
  # Calculate plot dimensions if not provided
  if (is.na(h) || is.na(w)) {
    plot_width <- usr[2] - usr[1]
    plot_height <- usr[4] - usr[3]
  }
  
  # Adjust size for orientation
  if (vertical) {
    h <- ifelse(is.na(h), plot_height * 0.75, h)
    w <- ifelse(is.na(w), plot_width * 0.05, w)
  } else {
    h <- ifelse(is.na(h), plot_height * 0.05, h)
    w <- ifelse(is.na(w), plot_width * 0.75, w)
  }
  
  # Calculate legend position if not provided
  if (is.na(x) || is.na(y)) {
    if (vertical) {
      x <- usr[2] - w * 1.5
      y <- usr[4] - w * 0.5 / xyr
    } else {
      y <- usr[4] - h * 1.5
      x <- usr[2] - w - h * 1.5 * xyr
    }
  }
  
  # Adjust the color spectrum
  if (length(colors) != n_rect) {
    colors <- colorRampPalette(colors)(n_rect)
  }
  
  # Add rectangles
  rect_width <- w / n_rect
  rect_height <- h / n_rect
  
  for (i in 1:n_rect) {
    if (vertical) {
      rect_x_left <- x
      rect_x_right <- x + w
      rect_y_bottom <- y - (i-1) * rect_height
      rect_y_top <- y - i * rect_height
    } else {
      rect_x_left <- x + (i-1) * rect_width
      rect_x_right <- x + i * rect_width
      rect_y_bottom <- y
      rect_y_top <- y - h
    }
    
    rect(rect_x_left, rect_y_bottom, rect_x_right, rect_y_top, col = colors[i], border = NA, xpd = xpd)
  }
  rect(x, y-h, x+w, y, xpd = xpd)
  
  # Handling labels
  if (!all(is.na(labels)) && !all(is.na(positions))) {
    if (!all(is.na(n_labels)) && n_labels > 1) {
      
      if (log_scale) {
        min_val <- log10(min(labels))
        max_val <- log10(max(labels))
      } else {
        min_val <- min(labels)
        max_val <- max(labels)
      }
      max_pos <- positions[which.max(labels)]
      min_pos <- positions[which.min(labels)]
      diff_minmax <- max_val - min_val
      diff_scale <- max_pos - min_pos
      slope <- diff_minmax / diff_scale
      
      label_positions <- seq(0, 1, length.out = n_labels)
      label_relative_positions <- label_positions - min_pos
      extrapolated_values <- min_val + label_relative_positions * slope
      if (log_scale) {
        extrapolated_values <- 10^(extrapolated_values)
      }
      
      #now make pretty if we want that
      if(log_scale){
        order_mag <- (10^floor(log10(extrapolated_values)))
        lab_logmod10 <- extrapolated_values / order_mag
        if(pretty){
          prettylab_logmod10 <- ceiling(lab_logmod10 - 1E-6)
          prettylab_logmod10[length(prettylab_logmod10)] <- floor(lab_logmod10[length(lab_logmod10)] + 1E-6)
          labels_values <- prettyScientificExpr(prettylab_logmod10 * order_mag, n_digits, 
                                                pretty = F, scientific = scientific)
        } else {
          labels_values <- prettyScientificExpr(lab_logmod10 * order_mag, n_digits, 
                                                pretty = F, scientific = scientific)
        }
        
      } else {
        labels_values <- prettyScientificExpr(extrapolated_values, n_digits, 
                                              pretty = pretty, scientific = scientific)  
      }
      
      
      label_positions <- sapply(labels_values$values, function(xi){
        if(log_scale){
          (log10(xi) - min_val) / slope + min_pos
        } else {
          (xi - min_val) / slope + min_pos
        }
      })
      
      #check that none of the labels_vals are outside the range of the scale
      in_range <- label_positions >= 0 & label_positions <= 1
      labels_values$values <- labels_values$values[in_range]
      labels_values$labels <- labels_values$labels[in_range]
      label_positions <- label_positions[in_range]
      label_vals <- labels_values[["labels"]]
      
      # #if we got rid of boundary labels, add them back in?
      # interlab_dist <- 1 / n_labels
      # boundary_lab_relative_dists <- c(label_positions[1], 1 - tail(label_positions, 1)) / interlab_dist
      # if(any(boundary_lab_relative_dists > 0.5)){
      #   labpos_to_add <- c(0,1)[boundary_lab_relative_dists > 0.5]
      #   label_positions <- c(label_positions, labpos_to_add)
      #   labels_values
      # }
      
      #plot the labels
      lb <- left_below
      for (i in seq_along(label_vals)) {
        label_pos <- label_positions[i]
        label_val <- label_vals[[i]]
        if (vertical) {
          text_x <- x - (strwidth(label_val) / 2 + w / 2) * ifelse(lb, 1, -1) + ifelse(lb, 0, w)
          text_y <- y - h + h * label_pos
          segments(x0 = x + w * ifelse(lb, 0, 1), 
                   x1 = x + ifelse(lb, -w / 4, 5 * w / 4), 
                   y0 = text_y, 
                   y1 = text_y, 
                   xpd = xpd)
        } else {
          text_x <- x + label_pos * w
          text_y <- y - (h * ifelse(lb, 1, 0) + strheight(lb) / 2 + h / 2) * ifelse(lb, 1, -1)
          segments(x0 = text_x, 
                   x1 = text_x, 
                   y0 = y - h * ifelse(lb, 1, 0), 
                   y1 = y + ifelse(lb, - 5 * h / 4, h/ 4),
                   xpd = xpd)
        }
        text(text_x, text_y, label_val, xpd = xpd)
      }
    }
  }
  
  #write in the title
  if(!is.na(main)){
    if(left_below){
      text(x = x + w / 2, y = y + strheight(main) / 4 * 3, labels = main, xpd = xpd, font = 2)    
    } else {
      text(x = x + w, y = y - h / 2, labels = main, xpd = xpd, font = 2, pos = 4)  
    }
    
  }
  
}