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

bezier_curve <- function(x0, y0, x1, y1, p = 0.5, n = 128, k = 1, col = 1, ...) {
  
  #handle vector inputs
  arg_vals <- as.list(environment())
  nargs <- length(arg_vals)
  arg_lens <- sapply(arg_vals, length)
  if(any(arg_lens > 1)){
    for(i in 1:max(arg_lens)){
      arg_inds <- (i - 1) %% arg_lens + 1
      arg_vals_i <- lapply(setNames(1:nargs, names(arg_vals)), function(j) arg_vals[[j]][arg_inds[j]])
      do.call(bezier_curve, arg_vals_i)
    }
  } else {
    # midpoint for curve
    xm <- x0 + p * (x1 - x0)
    ym <- y0 + p * (y1 - y0)
    
    # bezier control points
    control_x1 <- x0 + p * (x1 - x0) * k
    control_y1 <- y0
    control_x2 <- x1 - (1 - p) * (x1 - x0) * k
    control_y2 <- y1
    
    # calculate bezier curve
    t <- seq(0, 1, length.out = n)
    curve_x <- (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * control_x1 + 3 * (1 - t) * t^2 * control_x2 + t^3 * x1
    curve_y <- (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * control_y1 + 3 * (1 - t) * t^2 * control_y2 + t^3 * y1
    
    # plot
    lines(curve_x, curve_y, col = col, ...)
  }
}

my_heatmap <- function(mat_cols, mat_sizes = NULL, dim_names = NULL,
                       plot_diagonal_labels = T, middle_line_index_range = c(47,53),
                       cell_cols = c("blue", "white", "red"),
                       space_scale_ylabs = 1.4, reorder_mat = T,
                       symmetric_cols = T, diag_matters = F){
  
  if(is.null(mat_sizes)){
    mat_sizes <- mat_cols^0
  }
  
  if(reorder_mat){
    mat_order <- order(cmdscale(1.1 * max(abs(mat_cols)) - mat_cols, k = 1))
    mat_cols <- mat_cols[mat_order, mat_order]
    mat_sizes <- mat_sizes[mat_order, mat_order]
  }
  
  if(is.null(dim_names)){
    dim_names <- rownames(mat_cols)
  }
  
  #get plot dimensions
  nr <- nrow(mat_cols)
  nc <- ncol(mat_cols)
  
  # Generate indices for rows and columns
  rows <- matrix(rep(1:nr, each = nc), nr, nc)
  cols <- matrix(rep(1:nc, times = nr), nr, nc)
  
  # Create the plot area
  par(mar = c(2,18,8,2))
  plot(rows, cols, type = "n", xlab = "", ylab = "", 
       xlim = c(0.5, nr + 0.5), ylim = c(0.5, nc + 0.5), xaxt = "n", yaxt = "n", frame = F)
  
  # Adding a color palette (you can customize this)
  colors <- colorRampPalette(cell_cols)(100)
  
  # Map mat_cols values to color indices
  if(symmetric_cols){
    
    max_abs_val <- ifelse(diag_matters, max(abs(mat_cols)), max(abs(mat_cols[upper.tri(mat_cols)])))  # Maximum absolute value for scaling
    breaks <- seq(-max_abs_val, max_abs_val, length.out = 101) * (1+1E-6)
  } else {
    range_mat_cols <- c(ifelse(diag_matters, min(mat_cols), min(mat_cols[upper.tri(mat_cols)])),
                        ifelse(diag_matters, max(mat_cols), max(mat_cols[upper.tri(mat_cols)])))
    breaks <- seq(range_mat_cols[1], range_mat_cols[2], length.out = 101) * (1+1E-6)
  }
  color_indices <- matrix(cut(mat_cols, breaks = breaks, labels = FALSE), nr, nc)
  
  # Normalize mat_sizes for size mapping (scaling from range [0.5, 1])
  size <- abs(mat_sizes - 0.5) * 2
  
  # Plot the heatmap-like rectangles
  for (i in 1:nr^2) {
    x_center <- cols[i]
    y_center <- rows[i]
    
    #skip plotting diagonal if it does not matter (eg, in correlation matrices)
    if(!diag_matters){
      if(cols[i] == rows[i]){
        next()
      }
    }
    
    half_size <- size[i] / 2  # Half of the side length to adjust for centering
    
    # Draw rectangles
    rect(x_center - half_size, y_center - half_size, 
         x_center + half_size, y_center + half_size, 
         col = colors[color_indices[i]], border = NA)
  }
  
  #label heatmap
  
  #left labels (for names)
  xlocs_ll <- par("usr")[1] - diff(par("usr")[1:2])/10
  ylocs_ll <- (1:nr - nr/2) * space_scale_ylabs + nr / 2 * space_scale_ylabs - nr/10
  text(x = xlocs_ll, y = ylocs_ll, labels = paste0(nr:1, ". ", dim_names), 
       xpd = NA, pos = 2, col = 1, cex = 0.8)
  segments(x0 = -0.5, x1 = xlocs_ll - 0.5, y0 = 1:nr, y1 = ylocs_ll, xpd = NA)
  
  #top labels
  text(1:nr, nr + rep(0:1, ceiling(nr/2))[1:nr], 
       1:nr, pos = 3, cex = 0.5, xpd = NA)
  segments(x0 = 1:(nr/2) * 2, 
           x1 = 1:(nr/2) * 2, 
           y0 = nr+0.75,
           y1 = nr+1.75,
           xpd = NA)
  
  #bottom labels
  text(1:nr, 1 - rep(0:1, ceiling(nr/2))[1:nr], 
       1:nr, pos = 1, cex = 0.5, xpd = NA)
  segments(x0 = 1:(nr/2) * 2, 
           x1 = 1:(nr/2) * 2, 
           y0 = 0.25,
           y1 = -0.75,
           xpd = NA)
  
  #diagonal labels?
  if(plot_diagonal_labels){
    #get the indices of the rows that are blank enough and contain the diagonal
    coords_middle_lines <- do.call(rbind, lapply(1:nr, function(i){
      cis <- color_indices[,i]
      rlecont <- rle(cis > middle_line_index_range[1] & cis < middle_line_index_range[2])
      rlec <- data.frame(val = rlecont$values, len = rlecont$lengths)
      rlec$i0 <- cumsum(c(1, rlecont$lengths[-nrow(rlec)]))
      rlec$i1 <- cumsum(rlecont$lengths)
      rlec <- rlec[rlec$val,]
      rlec_mid <- rlec[rlec$i0 <= i & rlec$i1 >= i,]
      if(nrow(rlec_mid) > 0){
        out <- data.frame(x0 = i, x1 = i, y0 = rlec_mid$i0, y1 = rlec_mid$i1)  
      } else {
        out <- NULL
      }
      return(out)
    }))
    
    #draw guiding lines
    if(!is.null(coords_middle_lines)){
      segments(x0 = coords_middle_lines$x0, 
               y0 = coords_middle_lines$y0, 
               x1 = coords_middle_lines$x1, 
               y1 = coords_middle_lines$y1,
               adjustcolor(1, 0.1))
      
      #plot the labels 
      rect(xleft = 1:nr - 0.5, ybottom = 1:nr - 0.5, 
           xright = 1:nr + 0.5, ytop = 1:nr + 0.5, 
           col = colors[round(length(colors)/2)], border = NA)
      text(1:nr, 1:nr, 1:nr, cex = 0.5, xpd = NA, col = adjustcolor(1, 0.5))
    }
    
  }
  
  #top legend
  add_continuous_legend(colors = colors, labels = breaks[-1] + diff(breaks)/2, 
                        x = 1, y = par("usr")[4] + diff(par("usr")[3:4])/10, 
                        vertical = F, xpd = NA, positions = 1:100/100, left_below = F,
                        scientific = F, main = pn, w = diff(par("usr")[1:2])/2)
  
  
}

my_ridgeline <- function(
    d,
    ndens = 256,             # Number of density points
    kern_adj = 2,            # Kernel adjustment for density smoothing
    scale_vertical = 16,     # Vertical scaling factor for KDEs
    dens_thres = 1 / 256 / 4, # Density threshold for polygon display
    mass_thresh = 0.999,     # Threshold for cumulative density mass
    hard_bounds_range = NULL, # Hard x-axis bounds for density plots
    colors = c("brown2", "darkviolet", "deepskyblue4"), # Color palette
    min_sample_size = 20,
    decreasing = F,
    xlab = NULL,
    space_scale_ylabs = 1.4,
    ylab_disp = 0.1, 
    alternate_border_col = NULL,
    ylabs = T,
    panel_coords = NULL,
    order_inputs = NULL, 
    resc_to_maxh = F
) {
  
  # Filter out bins with low sample size
  d <- d[sapply(d, length) >= min_sample_size]
  n <- length(d) 
  
  # Get the sorted group levels 
  if(is.null(order_inputs)){
    group_order <- order(sapply(d, mean), decreasing = decreasing)
    d <- d[group_order]
  } else {
    d <- d[order_inputs]
  }
  
  # Calculate densities for each group
  m_dens <- lapply(d, function(x) {
    if (length(x) < min_sample_size) return(NULL)
    # Adjust for smoothing
    # x <- x + runif(length(x), 0, 1) - runif(length(x), 0, 1)
    out <- density(x, n = ndens, adjust = kern_adj)
    return(data.frame(x = out$x, y = out$y))
  })
  
  # Determine plotting range and color gradient
  x_range <- range(sapply(d, range), na.rm = TRUE)
  if (is.null(hard_bounds_range)) {
    hard_bounds_range <- x_range
  }
  
  #get color palette
  color_palette <- colorRampPalette(colors)(100)
  if(!is.null(alternate_border_col)){
    border_cols <- rep(c("black", alternate_border_col), ceiling(n/2))[1:n]
  } else {
    border_cols <- rep("black", n)
  }
  
  #convert to subpanel coordinates
  pusr <- c(hard_bounds_range, par("usr")[3:4])
  if(is.null(panel_coords)){
    pc <- pusr
  } else {
    pc <- panel_coords
  }
  dx <- function(x){
    wusr <- diff(pusr[1:2])
    pw <- diff(pc[1:2])
    ((x - pusr[1]) / wusr * pw) + pc[1]
  }
  dy <- function(y){
    husr <- diff(pusr[3:4])
    ph <- diff(pc[3:4])
    ((y - pusr[3]) / husr * ph) + pc[3]
  }
  
  # Draw KDE polygons
  lbs <- rbs <- yb <- numeric(n)
  for (i in n:1) {
    
    ad <- m_dens[[i]]
    if (is.null(ad)) next
    
    # Determine the density subset based on thresholds
    ad_cumul_mass <- cumsum((ad$y[-1] + ad$y[-nrow(ad)]) / 2 * diff(ad$x))
    ad_subset <- ad_cumul_mass > (1 - mass_thresh) / 2 & 
      ad_cumul_mass < mass_thresh + (1 - mass_thresh) / 2
    ad_subset <- c(ad_subset[1], ad_subset)
    ad_subset <- ad_subset & ad$x > hard_bounds_range[1] & ad$x < hard_bounds_range[2]
    
    ad_subset_dens_thresh <- rep(FALSE, nrow(ad))
    ad_subset_dens_thresh[min(which(ad$y > dens_thres)):max(which(ad$y > dens_thres))] <- TRUE
    ad_subset <- ad_subset & ad_subset_dens_thresh
    
    ad <- ad[ad_subset, ]
    scale_vertical_adj <- scale_vertical / sum(dy(ad$y[-1] + diff(ad$y)/2) * diff(dx(ad$x)))
    if(resc_to_maxh){
      max(ad$y)
      y_coords <- dy(i + c(ad$y, rep(0, nrow(ad))) / max(ad$y) * scale_vertical_adj)
    } else {
      y_coords <- dy(i + c(ad$y, rep(0, nrow(ad))) * scale_vertical_adj)
    }
    
    print(max(y_coords))
    
    polygon(x = dx(c(ad$x, rev(ad$x))),
            y = y_coords, 
            border = border_cols[i], col = color_palette[i])
    
    lbs[i] <- ad$x[1]
    rbs[i] <- ad$x[length(ad$x)]
    yb[i] <- i
  }

  # Add axis and labels
  hax_xlocs <- pretty(c(lbs, rbs))
  hax_xlocs <- hax_xlocs[hax_xlocs > min(lbs) & hax_xlocs < max(rbs)]
  hax_ylocs <- 0
  lines(dx(range(hax_xlocs)), rep(dy(hax_ylocs), 2))
  husr <- diff(pusr[3:4])
  tickh <- husr/50
  segments(x0 = dx(hax_xlocs), x1 = dx(hax_xlocs), y0 = dy(0), y1 = dy(-tickh))
  text(x = dx(hax_xlocs), y = dy(hax_ylocs - tickh), labels = hax_xlocs, pos = 1)
  text(x = dx(mean(range(hax_xlocs))), 
       y = dy(hax_ylocs - tickh * 2 - strheight("0", "user") * 1), 
       labels = xlab, pos = 1)
  
  #left labels (for names)
  if(ylabs){
    xlocs_ll <- dx(min(lbs) - diff(range(c(lbs, rbs))) / 5
    )
    ylocs_ll <- dy((1:n - n/2) * space_scale_ylabs + n / 2 * space_scale_ylabs - n/10 - 
      mean(ylocs_ll) + n/2 + ylab_disp * diff(par("usr")[3:4]))
    text(x = xlocs_ll, y = ylocs_ll, labels = paste0(n:1, ". ", names(d)), 
         xpd = NA, pos = 2, col = border_cols, cex = 0.8)
    
    #connecting segments -- first drawing the line over
    segments(x0 = dx(min(lbs)), x1 = dx(lbs), 
             y0 = dy(1:nr), y1 = dy(1:nr), xpd = NA, col = border_cols)
    #then using a bezier curve to smoothly connect to the label
    bezier_curve(x0 = dx(min(lbs)), x1 = xlocs_ll, 
                 y0 = dy(1:nr), y1 = ylocs_ll, k = 1, col = border_cols)
  }
  
  return(data.frame(i = 1:50, name = names(d), lbs = dx(lbs), rbs = dx(rbs), 
                    yb = dy(yb), bc = border_cols))
  
}

connect_ridges <- function(bounds, lab_cex = 0.7){
  
  n_plots <- length(bounds)
  name_inds <- setNames(bounds[[1]]$i, bounds[[1]]$name)
  n_dens <- length(name_inds)
  
  for(i in 2:n_plots){
    lp <- bounds[[i-1]]
    rp <- bounds[[i]]
    
    #connecting segments -- first drawing rbs and lbs over
    segments(x0 = lp$rbs, x1 = max(lp$rbs), 
             y0 = lp$yb, y1 = lp$yb, xpd = NA, col = lp$bc)
    segments(x0 = rp$lbs, x1 = min(rp$lbs), 
             y0 = rp$yb, y1 = rp$yb, xpd = NA, col = rp$bc)
    
    
    #then using a bezier curve to smoothly connect to the label
    mp <- merge(lp, rp, by = "name")
    bezier_curve(x0 = max(mp$rbs.x), x1 = min(mp$lbs.y), 
                 y0 = mp$yb.x, y1 = mp$yb.y, k = 1, col = mp$bc.y)
    
    #label the inds of the points on the right set of plots
    lab_x <- (min(rp$lb) + rp$lb)/2
    lab_h <- strheight(s = name_inds[rp$name], units = "user", cex = lab_cex)
    lab_w <- strwidth(s = name_inds[rp$name], units = "user", cex = lab_cex) * 1.15
    rect(xleft = lab_x - lab_w/2, 
         xright = lab_x + lab_w/2, 
         ybottom = rp$yb - lab_h/2, 
         ytop = rp$yb + lab_h/2,
         col = "white", border = NA)
    text(x = lab_x, y = rp$yb, labels = name_inds[rp$name], cex = lab_cex)
    
    #finish off the segments at the right of the figure
    if(i == n_plots){
      segments(x0 = rp$rbs, x1 = max(rp$rbs), 
               y0 = rp$yb, y1 = rp$yb, xpd = NA, col = rp$bc)
      lab_x <- max(rp$rb)
      lab_h <- strheight(s = name_inds[rp$name], units = "user", cex = lab_cex)
      lab_w <- strwidth(s = name_inds[rp$name], units = "user", cex = lab_cex) * 1.15
      rect(xleft = lab_x - lab_w/2, 
           xright = lab_x + lab_w/2, 
           ybottom = rp$yb - lab_h/2, 
           ytop = rp$yb + lab_h/2,
           col = "white", border = NA)
      text(x = lab_x, y = rp$yb, labels = name_inds[rp$name], cex = lab_cex)
    }
    
    #add in left hand indices for good measure
    if(i == 2){
      lab_x <- (lp$lb + min(lp$lb))/2
      lab_h <- strheight(s = name_inds[rp$name], units = "user", cex = lab_cex)
      lab_w <- strwidth(s = name_inds[rp$name], units = "user", cex = lab_cex) * 1.15
      rect(xleft = lab_x - lab_w/2, 
           xright = lab_x + lab_w/2, 
           ybottom = rp$yb - lab_h/2, 
           ytop = rp$yb + lab_h/2,
           col = "white", border = NA)
      text(x = lab_x, y = lp$yb, labels = name_inds[lp$name], cex = lab_cex)
    }
    
    
  }
  
}

ifelse2 <- function(test, yes, no) if(test){return(yes)}else{return(no)}
