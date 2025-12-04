
#helpers 
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
                                  left_below = TRUE, log_scale = FALSE, cex.main = 1,
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
    if(left_below || vertical){
      text(x = x + w / 2, y = y + strheight(main) / 4 * 3, labels = main, xpd = xpd, font = 2, cex = cex.main)    
    } else {
      text(x = x + w, y = y - h / 2, labels = main, xpd = xpd, font = 2, pos = 4, cex = cex.main)  
    }
    
  }
  
}


xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

indvec <- function(xind, nxind, ind, n){
  out <- rep(x = nxind, n)
  out[ind] <- xind
  return(out)
}

optimize_ticks <- function(x, y,
                           target.ticks = 5,
                           range_tt = NULL,
                           ppin  = par("pin"),
                           pusr  = par("usr"),
                           weight_pow_tt = 0,
                           xyr   = NULL){
  if(is.null(xyr)){
    xyr <- xyrat()
  }
  if(is.null(range_tt)){
    range_tt <- (ceiling(target.ticks/2)):(floor(target.ticks*2))  
  }
  poss_ticks <- lapply(range_tt, function(tt_opt){
    match_ticks(x, y, tt_opt, 
                ppin, pusr, xyr)
  })
  vrs <- sapply(poss_ticks, function(posst) posst$visual_ratio)
  lvrs <- abs(log(vrs)) * (abs(log(range_tt) - log(target.ticks)) + 1)^weight_pow_tt
  best_tick_opts <- which((lvrs - min(lvrs)) < 1E-3)
  best_ticks <- best_tick_opts[which.min(abs(range_tt[best_tick_opts] - target.ticks))]
  return(poss_ticks[[best_ticks]])
}

match_ticks <- function(x, y,
                        target.ticks = 5,
                        ppin  = par("pin"),
                        pusr  = par("usr"),
                        xyr   = NULL,
                        high.u.bias = 2,
                        u5.bias     = 1){
  
  # number of inches available in x and y
  pinx <- ppin[1];  piny <- ppin[2]
  # user ranges
  ux   <- diff(pusr[1:2]);  uy   <- diff(pusr[3:4])
  
  # which axis is physically smaller?
  cramped <- if (pinx < piny) 1 else 2          # 1 = x, 2 = y
  rng      <- list(range(x, na.rm = TRUE),
                   range(y, na.rm = TRUE))
  
  ## 1. get "pretty" ticks for the cramped axis
  cramped_ticks <- pretty(rng[[cramped]], n = target.ticks,
                          high.u.bias = high.u.bias,
                          u5.bias     = u5.bias)
  
  ## 2. convert their spacing from user units to inches
  du_per_inch   <- c(ux / pinx, uy / piny)      # data‑units per inch
  inch_spacing  <- median(diff(cramped_ticks)) / du_per_inch[cramped]
  
  ## 3. desired spacing (in user units) on the other axis
  other         <- 3 - cramped                  # 2 if 1, 1 if 2
  desired_du    <- inch_spacing * du_per_inch[other]
  
  ## 4. how many intervals does that imply?
  n_other <- max(1, round(diff(rng[[other]]) / desired_du))
  
  ## 5. final "pretty" ticks for both axes
  other_ticks <- pretty(rng[[other]], n = n_other,
                        high.u.bias = high.u.bias,
                        u5.bias     = u5.bias)
  
  xticks <- if (cramped == 1) cramped_ticks else other_ticks
  yticks <- if (cramped == 2) cramped_ticks else other_ticks
  
  ## 6. visual similarity measure (ideal = 1)
  phys_spacing <- c(
    median(diff(xticks)) / du_per_inch[1],
    median(diff(yticks)) / du_per_inch[2]
  )
  visual_ratio <- phys_spacing[1] / phys_spacing[2]
  
  list(xticks = xticks,
       yticks = yticks,
       visual_ratio = visual_ratio)
}

choose_scatter_alpha <- function(x, y,
                                 cex.pts      = par("cex"),
                                 high_q      = 0.9,
                                 A_high      = 0.85,
                                 alpha_min   = 0.02,
                                 alpha_max   = 0.5) {
  # guard: no data
  if (length(x) == 0L || length(y) == 0L) {
    return(A_high)
  }
  
  # estimate point footprint in user units
  w_usr <- strwidth("o", cex = cex.pts)
  h_usr <- strheight("o", cex = cex.pts)
  if (!is.finite(w_usr) || w_usr <= 0 || !is.finite(h_usr) || h_usr <= 0) {
    return(A_high)
  }
  
  xr <- range(x, finite = TRUE)
  yr <- range(y, finite = TRUE)
  if (!all(is.finite(xr)) || !all(is.finite(yr))) {
    return(A_high)
  }
  
  # choose step sizes ~ point footprint, but ensure at least a few bins
  if (diff(xr) <= 0 || diff(yr) <= 0) {
    return(A_high)
  }
  
  nx_approx <- diff(xr) / w_usr
  ny_approx <- diff(yr) / h_usr
  
  if (nx_approx < 5) {
    step_x <- diff(xr) / 5
  } else {
    step_x <- w_usr
  }
  if (ny_approx < 5) {
    step_y <- diff(yr) / 5
  } else {
    step_y <- h_usr
  }
  
  bx <- seq(xr[1], xr[2] + 1e-9, by = step_x)
  by <- seq(yr[1], yr[2] + 1e-9, by = step_y)
  if (length(bx) < 2L || length(by) < 2L) {
    return(A_high)
  }
  
  xbin <- cut(x, breaks = bx, include.lowest = TRUE, right = FALSE)
  ybin <- cut(y, breaks = by, include.lowest = TRUE, right = FALSE)
  
  counts <- table(xbin, ybin)
  n_vec  <- as.numeric(counts)
  pos_idx <- n_vec > 0
  
  if (!any(pos_idx)) {
    return(A_high)
  }
  
  n_pos <- n_vec[pos_idx]
  n_high <- as.numeric(stats::quantile(n_pos, probs = high_q, type = 1))
  n_high <- max(1L, round(n_high))
  
  if (!is.finite(n_high) || n_high <= 0) {
    return(A_high)
  }
  
  # enforce sane target and compute alpha for high-density anchor
  if (A_high <= 0) {
    return(alpha_min)
  }
  if (A_high >= 1) {
    A_high <- 0.99
  }
  
  a <- 1 - (1 - A_high)^(1 / n_high)
  
  if (!is.finite(a) || a <= 0) {
    a <- alpha_min
  }
  
  a <- max(alpha_min, min(alpha_max, a))
  a
}


scatter_hist <- function(x, y,
                         # embellishment
                         highlight_pt      = NULL,
                         col.hpt           = "red",
                         plot_pt = T,
                         shade_hist_tails  = FALSE,
                         shade_right_tail = T,
                         label_tail_proportions = T,
                         # appearance / layout (your original defaults)
                         col.pts = adjustcolor(1, 0.5),
                         cex.pts = 1,
                         target.ticks   = 10,
                         weight_pow_tt = 0,
                         pad.prop       = 0.04,
                         gap.prop       = 0.06,
                         hist.prop      = 0.25,
                         pch.pt         = 19,
                         equal_axes     = FALSE,
                         asp            = NULL,
                         # histogram + grid tweaks
                         equal_dim_hist_breaks = TRUE,
                         equal_mass_hists      = TRUE,
                         minor_tick_ratio      = 0.5,
                         draw_tick_grid        = TRUE,
                         xlab = expression(X[1]),
                         ylab = expression(X[2]),
                         labgap_scales = c(x = 1, y = 1),
                         xlim = NULL,
                         ylim = NULL,
                         xbreaks = NULL,
                         ybreaks = NULL,
                         add_one_to_one_line = F,
                         cex.lab = 1.5,
                         mar = c(4,4,1,1),
                         plot_as_heatmap = F,
                         kde_heatmap = F,
                         heatmap_max_dens_col = "steelblue",
                         ...) {
  
  
  #initiallize plotting window
  stopifnot(length(x) == length(y))
  op <- par(c("mar", "xpd"))
  on.exit(par(op), add = TRUE)
  
  par(mar = mar, xpd = NA)
  plot.new()
  ppin <- par("pin")
  propin.xy <- ppin[1] / ppin[2]
  
  ## ranges & layout maths
  xr <- range(x)
  yr <- range(y)
  
  dxr <- diff(xr)
  dyr <- diff(yr)
  
  xm <- mean(xr)
  ym <- mean(yr)
  
  pad.x <- pad.prop * dxr
  pad.y <- pad.prop * dyr * propin.xy
  
  x0 <- min(x) - pad.x
  x1 <- max(x) + pad.x
  y0 <- min(y) - pad.y
  y1 <- max(y) + pad.y
  
  if (equal_axes) {
    x0 <- y0 <- min(c(x0, y0)) 
    x1 <- y1 <- max(c(x1, y1))
  }
  
  hist.sides <- hist.prop * c(x = dxr + pad.x * 2,
                              y = (dyr + pad.y * 2) * propin.xy)
  
  if(!is.null(xlim)){
    x0 <- xlim[1]
    x1 <- xlim[2]
    dxr <- diff(xlim)
    pad.x <- pad.prop * dxr
    hist.sides["x"] <- hist.prop * c(x = dxr + pad.x * 2)
  }
  
  if(!is.null(ylim)){
    y0 <- ylim[1]
    y1 <- ylim[2]
    dyr <- diff(ylim)
    pad.y <- pad.prop * dyr * propin.xy
    hist.sides["y"] <- hist.prop * c(y = dyr + pad.y * 2)
  }
  
  gaps <- gap.prop * c(x = dxr + pad.x * 2,
                       y = (dyr + pad.y * 2) * propin.xy)
  
  #initialize plot
  xlim <- c(x0, x1 + gaps["x"] + hist.sides["x"])
  ylim <- c(y0, y1 + gaps["y"] + hist.sides["y"])
  plot.window(xlim = xlim,
              ylim = ylim,
              asp = asp)
  pusr <- par("usr"); ppin <- par("pin"); xyr <- xyrat()
  
  ## choose “pretty” breaks shared by scatter & hists
  if (equal_dim_hist_breaks) {
    ok <- optimize_ticks(x = x, y = y, 
                         target.ticks = target.ticks,
                         weight_pow_tt = weight_pow_tt, range_tt = NULL, 
                         ppin = ppin, pusr = pusr, xyr = xyr)
    xticks_orig <- ok$xticks; yticks_orig <- ok$yticks
  } else {
    xticks_orig <- pretty(x, n = target.ticks)
    yticks_orig <- pretty(y, n = target.ticks)
  }
  xticks <- xticks_orig[xticks_orig >= x0 & xticks_orig <= x1]
  yticks <- yticks_orig[yticks_orig >= y0 & yticks_orig <= y1]
  
  # or substitute in use-specified breaks
  if(!is.null(xbreaks)){
    xticks <- xticks_orig <- xbreaks;
  }
  if(!is.null(ybreaks)){
    yticks <- yticks_orig <- ybreaks;
  }
  
  ## histogram objects
  hx <- hist(x, breaks = xticks_orig, plot = FALSE)
  hx$mass <- hx$counts / sum(hx$counts)
  hx$area <- sum(hx$mass * diff(hx$breaks))
  
  hy <- hist(y, breaks = yticks_orig, plot = FALSE)
  hy$mass <- hy$counts / sum(hy$counts)
  hy$area <- sum(hy$mass * diff(hy$breaks))
  
  ## highlight bookkeeping
  if (!is.null(highlight_pt)) {
    stopifnot(length(highlight_pt) == 2, is.numeric(highlight_pt))
    redx <- highlight_pt[1]; redy <- highlight_pt[2]
  }
  
  #draw histograms
  
  #ensure bars fit in the alotted space and have the same area overall
  if(equal_mass_hists){
    hyhx_area_rat <- hy$area / hx$area
    hx_heights <- (hx$counts / sum(hx$counts)) * hyhx_area_rat
    hy_widths <- (hy$counts / sum(hy$counts))
    hx_max_hrats <- max(hx_heights / hist.sides["y"])
    hy_max_wrats <- max(hy_widths / hist.sides["x"])
    scale_mass <- max(hx_max_hrats, hy_max_wrats)
    hx_heights <- hx_heights / scale_mass
    hy_widths <- hy_widths / scale_mass
  } else {
    hx_heights <- hx$counts / sum(hx$counts)
    hy_widths <- hy$counts / sum(hy$counts)
    hx_max_hrats <- max(hx_heights / hist.sides["y"])
    hy_max_wrats <- max(hy_widths / hist.sides["x"])
    hx_heights <- hx_heights / hx_max_hrats
    hy_widths <- hy_widths / hy_max_wrats
  }
  
  y.base <- y1 + gaps["y"]
  for (i in seq_along(hx$counts)) {
    h <- hx_heights[i]
    if (h == 0) next
    col.i <- "grey80"; bord <- "grey40"
    if (!is.null(highlight_pt) && shade_hist_tails && 
        ifelse(shade_right_tail, redx <= hx$breaks[i], redx > hx$breaks[i])) {
      col.i <- adjustcolor(col.hpt, 0.5); bord <- adjustcolor(col.hpt, 0.8)
    }
    rect(hx$breaks[i], y.base, hx$breaks[i + 1], y.base + h,
         col = col.i, border = bord)
  }
  
  x.base <- x1 + gaps["x"]
  for (i in seq_along(hy$counts)) {
    w <- hy_widths[i]
    if (w == 0) next
    col.i <- "grey80"; bord <- "grey40"
    if (!is.null(highlight_pt) && shade_hist_tails && 
        ifelse(shade_right_tail, redy <= hy$breaks[i], redy > hy$breaks[i])) {
      col.i <- adjustcolor(col.hpt, 0.5); bord <- adjustcolor(col.hpt, 0.8)
    }
    rect(x.base, hy$breaks[i], x.base + w, hy$breaks[i + 1],
         col = col.i, border = bord)
  }
  
  if (label_tail_proportions && !is.null(highlight_pt)) {
    
    # left-tail probabilities
    prop.x.num <- mean(x < highlight_pt[1])
    prop.y.num <- mean(y < highlight_pt[2])
    if(prop.x.num < 0.001){
      prop.x <- prettyScientificExpr(prop.x.num)$labels  
    } else {
      prop.x <- paste0(round(prop.x.num * 100, 1), "%")
    }
    if(prop.y.num < 0.01){
      prop.y <- prettyScientificExpr(prop.y.num)$labels
    } else {
      prop.y <- paste0(round(prop.y.num * 100, 1), "%")
    }
    
    # midpoints of tail regions (left tails)
    tail.x.lo <- min(hx$breaks)
    tail.y.lo <- min(hy$breaks)
    midp.x <- mean(c(highlight_pt[1], tail.x.lo))
    midp.y <- mean(c(highlight_pt[2], tail.y.lo))
    
    # string extents in user units
    # for X label (horizontal text)
    w.x <- strwidth(prop.x)
    h.x <- strheight(prop.x)
    
    # for Y label, we will rotate by 90 degrees,
    # so its vertical extent comes from the width, scaled by xyr
    w.y <- strwidth(prop.y)
    h.y <- strheight(prop.y)
    len.y <- w.y * xyr        # vertical length of rotated Y-label
    
    x.left  <- midp.x - w.x / 2
    x.right <- midp.x + w.x / 2
    
    idx.x <- which(
      hx$breaks[-length(hx$breaks)] < x.right &
        hx$breaks[-1]               > x.left
    )
    
    if (length(idx.x) == 0L) {
      top.x <- y.base
    } else {
      top.x <- y.base + max(hx_heights[idx.x])
    }
    
    # a small offset so the label clears the highest overlapped bar
    off.x <- h.x * 0.4
    
    text(x = midp.x,
         y = top.x + off.x,
         labels = prop.x,
         adj = c(0.5, 0),
         col = col.hpt)
    
    y.low  <- midp.y - len.y / 2
    y.high <- midp.y + len.y / 2
    
    idx.y <- which(
      hy$breaks[-length(hy$breaks)] < y.high &
        hy$breaks[-1]               > y.low
    )
    
    if (length(idx.y) == 0L) {
      right.y <- x.base
    } else {
      right.y <- x.base + max(hy_widths[idx.y])
    }
    
    # offset so the (rotated) label clears the widest overlapped bar
    off.y <- h.y * 0.4
    
    text(x = right.y + off.y,
         y = midp.y,
         labels = prop.y,
         srt = 270,
         adj = c(0.5, 0),
         col = col.hpt)
  }
  
  y_lb <- min(hy$breaks)
  x_lb <- min(hx$breaks)
  segments(x.base, y_lb, x.base, y.base, lwd = 1.5)
  segments(x_lb, y.base, x.base, y.base, lwd = 1.5)
  
  #### axes 
  prop_tick_len <- 0.05
  max.tl <- prop_tick_len * c(diff(pusr[3:4]), diff(pusr[1:2]))[which.max(ppin)]
  tick.len <- c(x = max.tl, y = max.tl * xyr)
  lab.step.x <- max(1, ceiling(length(xticks) / target.ticks * 2))
  x.lab.id <- seq(1, length(xticks), by = lab.step.x)
  xtl <- indvec(tick.len["x"], tick.len["x"] * minor_tick_ratio,
                x.lab.id, length(xticks))
  segments(xticks, y0, xticks, y0 - xtl)
  text(xticks[x.lab.id], y0 - tick.len["x"], labels = xticks[x.lab.id], pos = 1)
  text(mean(c(x0, x1)), y0 - tick.len["x"] - max(strheight(xticks[x.lab.id])) * 2 - gaps["y"] * labgap_scales["y"],
       xlab, font = 1, cex = cex.lab)
  
  if(equal_axes){
    y.lab.id <- x.lab.id  
  } else {
    lab.step.y <- max(1, ceiling(length(yticks) / target.ticks * 2))
    lab.step.y <- lab.step.x
    y.lab.id <- seq(1, length(yticks), by = lab.step.y)
  }
  
  ytl <- indvec(tick.len["y"], tick.len["y"] * minor_tick_ratio,
                y.lab.id, length(yticks))
  segments(x0, yticks, x0 - ytl, yticks)
  text(x0 - tick.len["y"], yticks[y.lab.id], yticks[y.lab.id], adj = 1, pos = 2)
  text(x0 - tick.len["y"] - max(strwidth(yticks[y.lab.id])) * 1.5 - gaps["x"] * labgap_scales["x"],
       mean(c(y0, y1)), ylab, font = 1, cex = cex.lab, srt = 90)
  
  #### scatter, grid and highlight 
  if (plot_as_heatmap) {
    # bin x,y on the same grid used for the marginal histograms
    xbin <- cut(x, breaks = hx$breaks, include.lowest = TRUE, right = FALSE)
    ybin <- cut(y, breaks = hy$breaks, include.lowest = TRUE, right = FALSE)
    
    zmat <- xtabs(~ xbin + ybin, drop.unused.levels = FALSE)
    
    # mass per cell and log10 mass
    if (kde_heatmap) {
      nx <- nrow(zmat)
      ny <- ncol(zmat)
      
      dens <- MASS::kde2d(
        x, y, h = c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y)) * 0.1,
        n    = c(nx, ny),
        lims = c(min(hx$breaks), max(hx$breaks),
                 min(hy$breaks), max(hy$breaks))
      )
      
      mass_mat <- dens$z
      mass_vec <- as.numeric(mass_mat)
      mass_vec <- mass_vec / sum(mass_vec)
    } else {
      mass_vec <- as.numeric(zmat) / length(x)
    }
    pos_idx <- mass_vec > 0
    
    if (any(pos_idx)) {
      logmass_vec <- rep(NA_real_, length(mass_vec))
      logmass_vec[pos_idx] <- log10(mass_vec[pos_idx])
      
      logmin <- min(logmass_vec[pos_idx])
      logmax <- max(logmass_vec[pos_idx])
      
      ncol_hm <- 100
      hm_cols <- colorRampPalette(c("white", heatmap_max_dens_col))(ncol_hm)
      
      col_index_all <- rep(NA_integer_, length(mass_vec))
      if (logmax > logmin) {
        log_breaks <- seq(logmin, logmax, length.out = ncol_hm + 1)
        log_cut <- cut(logmass_vec[pos_idx], breaks = log_breaks, include.lowest = TRUE)
        col_index_all[pos_idx] <- as.integer(log_cut)
      } else {
        # all positive masses identical: just use the max-intensity color
        col_index_all[pos_idx] <- ncol_hm
      }
      
      # draw heatmap
      for (i in seq_len(nrow(zmat))) {
        for (j in seq_len(ncol(zmat))) {
          idx <- (j - 1L) * nrow(zmat) + i
          if (is.na(col_index_all[idx])) next
          rect(hx$breaks[i], hy$breaks[j],
               hx$breaks[i + 1L], hy$breaks[j + 1L],
               col = hm_cols[col_index_all[idx]], border = NA)
        }
      }
      
      # legend in log10(cell mass), colors aligned with heatmap,
      # placed just to the right of the furthest histogram bar,
      # with ticks/labels on the right side
      if (logmax > logmin) {
        usr_leg <- par("usr")
        pin_leg <- par("pin")
        xyr_leg <- diff(range(usr_leg[1:2])) / pin_leg[1] /
          diff(range(usr_leg[3:4])) * pin_leg[2]
        
        plot_width  <- usr_leg[2] - usr_leg[1]
        plot_height <- usr_leg[4] - usr_leg[3]
        
        leg_h <- plot_height * 0.5
        leg_w <- plot_width  * 0.05
        
        # initial vertical placement (top of legend)
        legend_y <- usr_leg[4] - leg_w * 0.5 / xyr_leg
        leg_y_bottom <- legend_y - leg_h
        leg_y_top    <- legend_y
        
        # which vertical histogram bars overlap legend vertically?
        y_low  <- hy$breaks[-length(hy$breaks)]
        y_high <- hy$breaks[-1]
        idx_overlap <- which(
          y_low  < leg_y_top &
            y_high > leg_y_bottom
        )
        
        if (length(idx_overlap) == 0L) {
          # if nothing overlaps, fall back to the base
          hist_far_right <- x.base
        } else {
          hist_far_right <- x.base + max(hy_widths[idx_overlap])
        }
        
        legend_gap <- tick.len["y"]
        legend_x   <- hist_far_right + legend_gap
        
        add_continuous_legend(
          colors    = rev(hm_cols),                 # reverse so high values are dark at the top
          labels    = c(logmin, logmax),            # log10 mass range
          positions = c(0, 1),
          x         = legend_x, y = 1.4, 
          h = leg_h, w = leg_w,
          left_below = F,                       # ticks and labels on the right
          main      = expression(log[10]("mass"))
        )
      }
    }
  } else {
    if(is.null(col.pts)){
      col.alpha <- choose_scatter_alpha(x, y, cex.pts = cex.pts)
      col.pts <- adjustcolor(1, col.alpha)
    }
    points(x, y, pch = pch.pt, col = col.pts)
  }
  rect(x0, y0, x1, y1)
  
  if (draw_tick_grid) {
    segments(xticks, y0, xticks, y1, lty = 3, col = adjustcolor(1, 0.25))
    segments(x0, yticks, x1, yticks, lty = 3, col = adjustcolor(1, 0.25))
  }
  
  if (!is.null(highlight_pt)) {
    if(plot_pt){
      points(redx, redy, pch = pch.pt, col = col.hpt)  
    }
    
    # guides up/right
    if(shade_hist_tails){
      segments(redx, redy, redx, y.base, col = col.hpt, lty = 3, lwd = 2)
      segments(redx, redy, x.base, redy, col = col.hpt, lty = 3, lwd = 2)
    }
  }
  
  if(add_one_to_one_line){
    xmin <- min(x0, x1); xmax <- max(x0, x1)
    ymin <- min(y0, y1); ymax <- max(y0, y1)
    pts_121 <- cbind(x = c(xmin, xmax, ymin, ymax),
                     y = c(xmin, xmax, ymin, ymax))[ 
                       c(ymin <= xmin & xmin <= ymax,
                         ymin <= xmax & xmax <= ymax,
                         xmin <= ymin & ymin <= xmax,
                         xmin <= ymax & ymax <= xmax), , drop = FALSE]
    segments(pts_121[1,1], pts_121[1,2], pts_121[2,1], pts_121[2,2],
             lty = 2, col = adjustcolor(2, 0.85), lwd = 2)
  }
}

#### simulate data ####
run_example <- F
if(run_example){
  n            <- 1000      # sample size
  r            <- -0.8       # correlation  (‑1 < r < 1)
  x <- rnorm(n)
  y <- r * x + sqrt(1 - r^2) * rnorm(n)      # correlated partner
  y <- y * 3
  scatter_hist(x, y, equal_mass_hists = T, 
               weight_pow_tt = 1, target.ticks = 15, 
               col.pts = adjustcolor(1,0.4), highlight_pt = c(1, 1.5), 
               shade_hist_tails = T)
}
