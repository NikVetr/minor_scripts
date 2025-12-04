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
                                  n_digits = 3, pretty = T, scientific = T, main = NULL,
                                  return_BB = F) {
  
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
      
      #hack to fix misaligned positions, need to fix properly later
      label_positions <- label_positions - min(label_positions)
      
      #check that none of the labels_vals are outside the range of the scale
      eps <- diff(label_positions)[1]/100
      in_range <- label_positions >= (0-eps) & label_positions <= (1+eps)
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
  if(!is.null(main)){
    if(left_below){
      text(x = x + w / 2, y = y + strheight(main) / 4 * 3, labels = main, xpd = xpd, font = 2, cex = cex.main)    
    } else {
      text(x = x + w, y = y - h / 2, labels = main, xpd = xpd, font = 2, pos = 4, cex = cex.main)  
    } 
  }
  
  #return bounding box if requested
  if (return_BB) {
    # start with color bar extent
    bb_xmin <- x
    bb_xmax <- x + w
    bb_ymin <- y - h
    bb_ymax <- y
    
    # include tick labels if they were drawn
    if (exists("label_vals") && length(label_vals) > 0) {
      lb <- left_below
      
      if (vertical) {
        for (i in seq_along(label_vals)) {
          label_pos <- label_positions[i]
          label_val <- label_vals[[i]]
          
          text_x <- x - (strwidth(label_val) / 2 + w / 2) * ifelse(lb, 1, -1) + ifelse(lb, 0, w)
          text_y <- y - h + h * label_pos
          
          lw <- strwidth(label_val)
          lh <- strheight(label_val)
          
          bb_xmin <- min(bb_xmin, text_x - lw / 2)
          bb_xmax <- max(bb_xmax, text_x + lw / 2)
          bb_ymin <- min(bb_ymin, text_y - lh / 2)
          bb_ymax <- max(bb_ymax, text_y + lh / 2)
        }
      } else {
        for (i in seq_along(label_vals)) {
          label_pos <- label_positions[i]
          label_val <- label_vals[[i]]
          
          text_x <- x + label_pos * w
          text_y <- y - (h * ifelse(lb, 1, 0) + strheight(lb) / 2 + h / 2) * ifelse(lb, 1, -1)
          
          lw <- strwidth(label_val)
          lh <- strheight(label_val)
          
          bb_xmin <- min(bb_xmin, text_x - lw / 2)
          bb_xmax <- max(bb_xmax, text_x + lw / 2)
          bb_ymin <- min(bb_ymin, text_y - lh / 2)
          bb_ymax <- max(bb_ymax, text_y + lh / 2)
        }
      }
    }
    
    # include title if present
    if (!is.null(main)) {
      if (left_below) {
        title_x <- x + w / 2
        title_y <- y + strheight(main, cex = cex.main) * 3 / 4
        mw <- strwidth(main, cex = cex.main)
        mh <- strheight(main, cex = cex.main)
        
        bb_xmin <- min(bb_xmin, title_x - mw / 2)
        bb_xmax <- max(bb_xmax, title_x + mw / 2)
        bb_ymin <- min(bb_ymin, title_y - mh / 2)
        bb_ymax <- max(bb_ymax, title_y + mh / 2)
      } else {
        title_x <- x + w
        title_y <- y - h / 2
        mw <- strwidth(main, cex = cex.main)
        mh <- strheight(main, cex = cex.main)
        
        bb_xmin <- min(bb_xmin, title_x)
        bb_xmax <- max(bb_xmax, title_x + mw)
        bb_ymin <- min(bb_ymin, title_y - mh / 2)
        bb_ymax <- max(bb_ymax, title_y + mh / 2)
      }
    }
    
    return(c(xmin = bb_xmin, xmax = bb_xmax, ymin = bb_ymin, ymax = bb_ymax))
  }
  
}


bezier_curve <- function(x0, y0, x1, y1,
                         p = 0.5, n = 128, k = 1,
                         ends = c("flat", "steep")[1],
                         col = 1, return_points = FALSE, debug = FALSE,
                         xpd = NA, ...) {
  # choose behavior: "flat" = horizontal tangents at ends (original);
  #                  "steep" = vertical tangents at ends (the “other direction”)
  if(is.numeric(ends)){
    ends <- c("flat", "steep")[ends]
  }
  ends <- match.arg(ends)  
  
  
  # handle vector inputs (recurse over elementwise args)
  arg_vals <- as.list(environment())
  nargs <- length(arg_vals)
  arg_lens <- sapply(arg_vals, length)
  if (any(arg_lens > 1)) {
    for (i in 1:max(arg_lens)) {
      arg_inds <- (i - 1) %% arg_lens + 1
      arg_vals_i <- lapply(setNames(1:nargs, names(arg_vals)),
                           function(j) arg_vals[[j]][arg_inds[j]])
      do.call(bezier_curve, arg_vals_i)
    }
    return(invisible(NULL))
  }
  
  # control points
  if (ends == "flat") {
    # horizontal tangents at endpoints (original behavior)
    control_x1 <- x0 + p * (x1 - x0) * k
    control_y1 <- y0
    control_x2 <- x1 - (1 - p) * (x1 - x0) * k
    control_y2 <- y1
  } else {
    # vertical tangents at endpoints (steep at start/end)
    control_x1 <- x0
    control_y1 <- y0 + p * (y1 - y0) * k
    control_x2 <- x1
    control_y2 <- y1 - (1 - p) * (y1 - y0) * k
  }
  
  if (debug) {
    cat("# control points:\n")
    cat(sprintf("# P0=(%.4f, %.4f)  P1=(%.4f, %.4f)\n", x0, y0, control_x1, control_y1))
    cat(sprintf("# P2=(%.4f, %.4f)  P3=(%.4f, %.4f)\n", control_x2, control_y2, x1, y1))
    # endpoint derivatives of a cubic Bézier: 3*(P1-P0) at t=0; 3*(P3-P2) at t=1
    d0x <- 3 * (control_x1 - x0); d0y <- 3 * (control_y1 - y0)
    d1x <- 3 * (x1 - control_x2); d1y <- 3 * (y1 - control_y2)
    cat(sprintf("# dB/dt at t=0  ≈ (%.4f, %.4f)\n", d0x, d0y))
    cat(sprintf("# dB/dt at t=1  ≈ (%.4f, %.4f)\n", d1x, d1y))
  }
  
  # sample curve
  t <- seq(0, 1, length.out = n)
  curve_x <- (1 - t)^3 * x0 +
    3 * (1 - t)^2 * t * control_x1 +
    3 * (1 - t) * t^2 * control_x2 +
    t^3 * x1
  curve_y <- (1 - t)^3 * y0 +
    3 * (1 - t)^2 * t * control_y1 +
    3 * (1 - t) * t^2 * control_y2 +
    t^3 * y1
  
  # draw
  lines(curve_x, curve_y, col = col, xpd = xpd, ...)
  
  if (return_points) {
    return(invisible(list(x = curve_x, y = curve_y)))
  } else {
    return(invisible(NULL))
  }
}

darken_rgb <- function(cols, amount = 0.3) {
  # amount: 0 = no change, 1 = full black
  if (amount < 0 || amount > 1) {
    stop("amount must be between 0 and 1")
  }
  
  cols_rgb <- col2rgb(cols)  # 3 x n matrix (r, g, b)
  
  # scale toward black
  scale_factor <- 1 - amount
  darker_rgb <- cols_rgb * scale_factor
  
  # clamp just in case
  darker_rgb[darker_rgb < 0] <- 0
  darker_rgb[darker_rgb > 255] <- 255
  
  darker_cols <- rgb(
    red   = darker_rgb[1, ],
    green = darker_rgb[2, ],
    blue  = darker_rgb[3, ],
    maxColorValue = 255
  )
  
  return(darker_cols)
}

my_heatmap <- function(mat_cols, mat_size_rule = NULL, dim_names = NULL,
                       plot_diagonal_labels = F, middle_line_index_range = c(47,53),
                       cell_cols = c("blue", "white", "red"),
                       space_scale_ylabs = 1.4, reorder_mat = T, 
                       cex.leg = 1, cex_lab = 1, lwd_for_lines = 1, 
                       symmetric_cols = T, diag_matters = F, legend_title = "", 
                       plot_labels = T, plot_legend = T, col_scale = NULL,
                       leg_loc_scale = c(x=0,y=0),
                       plot_guiding_lines = F, tlbr_diag = T, use_bezier = F,
                       mds_method = c("cailliez", "isoMDS", "smacof", "linear", "OLO", "hclust")[3],
                       blockmodel = F, highlight_blocks = T, blocks_to_highlight = "all",
                       darken_lc_by = 0.1, asp = NA, space_for_sep = 1,
                       demarcate_maps = T, outer_border = F, 
                       outer_border_brackets = T, border.lwd = 2,
                       plot_numbers = T, plot_rects = T, number_col = NULL,
                       main = NULL, cex.main = 1.25){
  
  #quick helper
  ifelse2 <- function(test, yes, no) if(test){return(yes)}else{return(no)}
  
  #do we have one matrix, or a list of matrices?
  multimat <- "list" %in% class(mat_cols)
  if(multimat){
    mat_cols_list <- mat_cols
  } else {
    mat_cols_list <- list(mat_cols)
  }
  
  #name rows / columns if unnamed
  if(is.null(colnames(mat_cols_list[[1]])) & 
     is.null(rownames(mat_cols_list[[1]]))){
    mat_cols_list <- lapply(mat_cols_list, function(mci){
      mci_new <- mci
      rownames(mci_new) <- paste0("dim_", 1:nrow(mci))
      colnames(mci_new) <- paste0("dim_", 1:ncol(mci))
      return(mci_new)
    })
  } else {
    stop("both rows and columns need to be named")
  }
  nmats <- length(mat_cols_list)
  mat_names <- names(mat_cols_list)
  mat_cols <- mat_cols_list[[1]] #do all the ordering etc. off the first matrix
  
  if(is.null(mat_size_rule)){
    mat_sizes_list <- lapply(mat_cols_list, function(ms) ms^0)
  } else if(mat_size_rule == "abs"){
    mat_sizes_list <- lapply(mat_cols_list, function(ms) abs(ms)^0.5)
  }
  mat_sizes <- mat_sizes_list[[1]]
  
  if(reorder_mat){
    dists <- max(mat_cols) - mat_cols
    if(mds_method == "cailliez") {
      fit <- cmdscale(as.dist(dists), k = 1, add = T)
      if("list" %in% class(fit)) fit <- fit$points
      x1d <- setNames(as.numeric(fit), rownames(fit))
    } else if(mds_method == "isoMDS") {
      fit <- MASS::isoMDS(as.dist(dists), k = 1)$points
      x1d <- setNames(as.numeric(fit), rownames(fit))
    } else if(mds_method == "smacof") {
      fit <- smacof::smacofSym(as.dist(dists), ndim = 1, type = "ordinal")$conf
      x1d <- setNames(as.numeric(fit), rownames(fit))
    } else if(mds_method == "linear") {
      x1d <- setNames(1:nrow(dists)/nrow(dists), rownames(dists))
    } else if(mds_method == "OLO") {
      tree <- seriation::seriate(as.dist(dists), method = "OLO")[[1]]
      x1d <- setNames(1:nrow(dists)/nrow(dists), tree$labels)
    } else if(mds_method == "hclust") { #not an mds method lol
      tree <- hclust(as.dist(dists))
      tree <- dendsort::dendsort(as.dendrogram(tree))
      tree <- as.hclust(tree)
      x1d <- setNames(1:nrow(dists)/nrow(dists), tree$labels)
    }
    mat_order <- names(x1d)[order(x1d)]
    mat_cols <- mat_cols[mat_order, mat_order]
    mat_sizes <- mat_sizes[mat_order, mat_order]
  } else {
    x1d <- setNames(0:(nrow(mat_cols)-1)/(nrow(mat_cols)-1),
                    rownames(mat_cols))
  }
  
  #specify color vector for labels
  lab_cols <- rep(1, nrow(mat_cols))
  
  #figure out range of values to omit leading 0s or not for cell labels
  eps <- 1E-6
  vals_for_range <- as.numeric(mat_cols[!is.na(mat_cols)])
  omit_leading_zero <- length(vals_for_range) > 0 &&
    all(vals_for_range >= -(1+eps) & vals_for_range <= (1+eps))
  
  #run blockmodeling algorithm
  lab_cols <- rep(1, nrow(mat_cols))
  lab_fonts <- rep(1, nrow(mat_cols))
  line_lwds <- rep(lwd_for_lines, nrow(mat_cols))
  if(blockmodel) {
    adj_mat <- mat_cols - min(mat_cols)
    diag(adj_mat) <- 0
    bm <- blockmodels::BM_gaussian(
      membership_type = "SBM_sym",
      adj            = adj_mat,
      verbosity      = 0,
      plotting       = ""
    )
    bm$estimate()
    best_q <- which.max(bm$ICL)
    membs <- bm$memberships[[best_q]]$Z
    block_ids <- setNames(max.col(membs, ties.method = "first"), 
                          rownames(adj_mat))
    # block_sets <- split(names(block_ids), block_ids)
    
    #order entries within and between blocks to mds result
    block_ids <- block_ids[names(x1d)]
    block_meanlocs <- sapply(split(x1d, block_ids), mean)
    block_ids <- setNames(order(block_meanlocs)[block_ids],
                          names(block_ids)) #between blocks
    block_factor <- factor(block_ids, levels = sort(unique(block_ids))) #within blocks
    mat_order <- names(block_factor)[order(block_factor, x1d)]
    
    #and reorder to final configuration
    block_ids <- block_ids[mat_order]
    mat_cols  <- mat_cols[mat_order, mat_order, drop = FALSE]
    mat_sizes <- mat_sizes[mat_order, mat_order, drop = FALSE]
    
    #get block colors too
    block_cols <- rep(1, best_q)
    block_fonts <- rep(1, best_q)
    lwd_lines <- rep(lwd_for_lines, best_q)
    block_cols_transp <- adjustcolor(block_cols, 0)
    if(highlight_blocks){
      block_cols_hl <- ifelse2(best_q > 2 & best_q < 9,
                            RColorBrewer::brewer.pal(best_q, "Dark2"),
                            rainbow(best_q))
      block_cols_transp_hl <- adjustcolor(block_cols_hl, 0.25)  
      if(all(blocks_to_highlight == "all")){
        block_cols <- block_cols_hl
        block_cols_transp <- block_cols_transp_hl
      } else {
        blocks_to_highlight <- pmax(pmin(blocks_to_highlight, best_q), 1)
        block_cols[blocks_to_highlight] <- block_cols_hl[blocks_to_highlight]
        block_cols_transp[blocks_to_highlight] <- block_cols_transp_hl[blocks_to_highlight]
        block_fonts[blocks_to_highlight] <- 2
        lwd_lines[blocks_to_highlight] <- lwd_lines * 1.5
      }
    }
    block_cols_dark <- darken_rgb(block_cols, amount = darken_lc_by)
    lab_cols <- block_cols_dark[block_ids]
    lab_fonts <- block_fonts[block_ids]
    line_lwds <- lwd_lines[block_ids]
  }
  
  
  if(is.null(dim_names)){
    dim_names <- rownames(mat_cols)
  }
  
  #get plot dimensions
  nr <- nrow(mat_cols)
  nc <- ncol(mat_cols)
  
  # Generate indices for rows and columns
  rows <- row(mat_cols)
  cols <- col(mat_cols)
  
  if (tlbr_diag) {
    # flip rows so that row 1 is at the top
    rows <- nr + 1 - rows
  }
  
  # Create the plot area
  # par(mar = c(2,18,8,2))
  if(multimat){
    plot(rows, cols, type = "n", xlab = "", ylab = "", 
         xlim = c(0.5, (nr + 1.5) * nmats), 
         ylim = c(0.5, nc + 0.5), 
         xaxt = "n", yaxt = "n", frame = F, asp = asp)
  } else {
    plot(rows, cols, type = "n", xlab = "", ylab = "", 
         xlim = c(0.5, nr + 0.5), ylim = c(0.5, nc + 0.5), 
         xaxt = "n", yaxt = "n", frame = F, asp = asp)  
  }
  pusr <- par("usr")
  w <- diff(pusr[1:2])
  h <- diff(pusr[3:4])
  
  # Adding a color palette (you can customize this)
  colors <- colorRampPalette(cell_cols)(100)
  
  # Map mat_cols values to color indices
  if(symmetric_cols){
    if(is.null(col_scale)){
      max_abs_val <- ifelse(diag_matters, 
                            max(abs(unlist(mat_cols_list))), 
                            max(abs(unlist(lapply(mat_cols_list, function(mc) mc[upper.tri(mc)]))))
                            )  # Maximum absolute value for scaling
    } else {
      max_abs_val <- max(abs(col_scale))
    }
    breaks <- seq(-max_abs_val, max_abs_val, length.out = 101) * (1+1E-6)
  } else {
    if(is.null(col_scale)){
      range_mat_cols <- c(ifelse(diag_matters, min(mat_cols), min(mat_cols[upper.tri(mat_cols)])),
                          ifelse(diag_matters, max(mat_cols), max(mat_cols[upper.tri(mat_cols)])))  
    } else {
      range_mat_cols <- col_scale
    }
    breaks <- seq(range_mat_cols[1], range_mat_cols[2], length.out = 101) * (1+1E-6)
  }
  #guard against constant matrix
  if(length(unique(breaks)) == 1){
    breaks <- seq(max_abs_val-1, max_abs_val+1, length.out = 101) * (1+1E-6)
  }
  
  
  #iterate through all of our listed matrices
  for(mi in 1:nmats){
    
    #get core level variables
    disp_x_by <- (mi-1) * (nr + space_for_sep)
    mat_cols <- mat_cols_list[[mi]]
    mat_sizes <- mat_sizes_list[[mi]]
    
    #permute everything correctly
    mat_cols  <- mat_cols[mat_order, mat_order, drop = FALSE]
    mat_sizes <- mat_sizes[mat_order, mat_order, drop = FALSE]
    
    #get cell colors
    color_indices <- matrix(cut(mat_cols, breaks = breaks, labels = FALSE), nr, nc)
    
    #get number colors
    if(is.null(number_col)){
      if(plot_rects){
        num_col_mat <- matrix("white", nrow = nr, ncol = nc)
      } else {
        num_col_mat <- matrix(colors[color_indices], nrow = nr, ncol = nc)
      }
    } else {
      num_col_mat <- matrix(number_col, nrow = nr, ncol = nc)
    }
    
    # Normalize mat_sizes for size mapping (scaling from range [0.5, 1])
    if(max(mat_sizes) > 1){
      size <- mat_sizes / max(mat_sizes)  
    } else {
      size <- mat_sizes
    }
    
    
    #label heatmap
    
    #left labels (for names)
    if(plot_labels){
      xlocs_ll <- -1/2 - diff(pusr[1:2])/10/nmats
      ylocs_ll <- (1:nr - nr/2) * space_scale_ylabs + nr / 2 * space_scale_ylabs - nr/10
      
      if(mi == 1){
        text(x = xlocs_ll, y = ylocs_ll, labels = rev(paste0(1:nr, ". ", dim_names)), 
             xpd = NA, pos = 2, cex = cex_lab * 0.8, col = rev(lab_cols), font = rev(lab_fonts))
        if(use_bezier){
          bezier_curve(x0 = -1/2, x1 = xlocs_ll - 0.5 * nr / 30, y0 = 1:nr, 
                       y1 = ylocs_ll, xpd = NA, col = rev(lab_cols), lwd = rev(line_lwds))
        } else {
          segments(x0 = -1/2, x1 = xlocs_ll - 0.5 * nr / 30, y0 = 1:nr, y1 = ylocs_ll, 
                   xpd = NA, col = rev(lab_cols), lwd = rev(line_lwds))
        }  
      }
      
      #top labels
      disp_size <- nr / 30
      text(1:nr + disp_x_by, nr + rep(0:1, ceiling(nr/2))[1:nr] * disp_size + 0.1, 
           labels = ifelse2(tlbr_diag, 1:nr, nr:1), pos = 3, cex = cex_lab * 0.5, xpd = NA, 
           col = lab_cols, font = lab_fonts)
      segments(x0 = 1:(nr/2) * 2 + disp_x_by, 
               x1 = 1:(nr/2) * 2 + disp_x_by, 
               y0 = nr + 0.75,
               y1 = nr + disp_size + 0.6,
               xpd = NA, col = lab_cols[1:(nr/2) * 2], 
               lwd = line_lwds[1:(nr/2) * 2])
      
      #bottom labels
      text(1:nr + disp_x_by, 1 - rep(0:1, ceiling(nr/2))[1:nr] * disp_size - 0.25, 
           labels = ifelse2(tlbr_diag, 1:nr, nr:1), pos = 1, cex = cex_lab * 0.5, xpd = NA, 
           col = lab_cols, font = lab_fonts)
      segments(x0 = 1:(nr/2) * 2 + disp_x_by, 
               x1 = 1:(nr/2) * 2 + disp_x_by, 
               y0 = 0.25,
               y1 = -disp_size + 0.4,
               xpd = NA, col = lab_cols[1:(nr/2) * 2], 
               lwd = line_lwds[1:(nr/2) * 2])
      
    }
    
    
    #diagonal labels?
    if (plot_diagonal_labels) {
      
      # find vertical middle-band segments that contain the diagonal element in each column
      coords_middle_lines <- do.call(rbind, lapply(1:nr, function(i) {
        cis <- color_indices[, i]
        
        in_mid <- cis > middle_line_index_range[1] & cis < middle_line_index_range[2]
        
        rlecont <- rle(in_mid)
        rlec <- data.frame(val = rlecont$values, len = rlecont$lengths)
        
        rlec$i0 <- cumsum(c(1, rlec$len[-nrow(rlec)]))
        rlec$i1 <- cumsum(rlec$len)
        
        rlec <- rlec[rlec$val, ]
        
        # keep only segments that contain the diagonal row i
        rlec_mid <- rlec[rlec$i0 <= i & rlec$i1 >= i, ]
        
        if (nrow(rlec_mid) > 0) {
          data.frame(x0 = i, x1 = i,
                     y0 = rlec_mid$i0, y1 = rlec_mid$i1)
        } else {
          NULL
        }
      }))
      
      #full guiding grid (under everything)
      if (plot_guiding_lines) {
        #horiz guides
        segments(x0 = ifelse(mi==1, -1/2, 0) + disp_x_by, 
                 x1 = nr + disp_x_by + 1/2, 
                 y0 = 1:nr, 
                 y1 = 1:nr,
                 col = adjustcolor(1, 0.5), lty = 3)
        #vert guides
        segments(x0 = 1:nr + disp_x_by, 
                 x1 = 1:nr + disp_x_by, 
                 y0 = 0.25, 
                 y1 = nr + 0.75,
                 col = adjustcolor(1, 0.5), lty = 3)
      }
      
      #extra vertical guides along "blank" parts that cross the diagonal
      if (!is.null(coords_middle_lines)) {
        
        # map matrix row indices to plotting y-coordinates if tlbr_diag flips rows
        if (tlbr_diag) {
          y0_plot <- nr + 1 - coords_middle_lines$y1
          y1_plot <- nr + 1 - coords_middle_lines$y0
        } else {
          y0_plot <- coords_middle_lines$y0
          y1_plot <- coords_middle_lines$y1
        }
        
        segments(x0 = coords_middle_lines$x0 + disp_x_by,
                 y0 = y0_plot,
                 x1 = coords_middle_lines$x1 + disp_x_by,
                 y1 = y1_plot,
                 col = adjustcolor(1, 0.1))
      }
    }
    
    #actually plot the diagonal labels
    if (plot_diagonal_labels && !is.null(color_indices)) {
      
      # compute diagonal positions in plotting coordinates
      diag_x <- 1:nr
      
      if (tlbr_diag) {
        diag_y <- nr + 1 - (1:nr)
        diag_lab <- 1:nr
      } else {
        diag_y <- 1:nr
        diag_lab <- nr:1
      }
      
      # neutral background on the diagonal
      rect(xleft = diag_x - 0.5 + disp_x_by, xright = diag_x + 0.5 + disp_x_by,
           ybottom = diag_y - 0.5, ytop = diag_y + 0.5,
           col = colors[round(length(colors) / 2)], border = NA)
      
      # diagonal index labels
      text(diag_x + disp_x_by, diag_y, diag_lab,
           cex = cex_lab * 0.5, xpd = NA,
           col = ifelse2(blockmodel, lab_cols, 
                         adjustcolor(1, 0.5)),
           font = ifelse2(blockmodel, lab_fonts, 
                         1),
      )
    }
    
    
    #plot block highlights
    if(blockmodel & mi == 1){
      for(bi in 1:best_q){
        block_inds <- range(which(block_ids == bi))
        
        if(tlbr_diag) {
          ybottom <- nr - block_inds[1] + 1 + 1/2
          ytop    <- nr - block_inds[2] + 1 - 1/2
        } else {
          ybottom <- block_inds[1] - 1/2
          ytop    <- block_inds[2] + 1/2
        }
        
        rect(xleft = block_inds[1] - 1/2,
             xright = block_inds[2] + 1/2,
             ybottom = ybottom,
             ytop = ytop,
             col = block_cols_transp[bi],
             border = NA
        )
      }
    }
    
    # Plot the heatmap-like rectangles
    for (i in 1:nr^2) {
      x_center <- cols[i]
      y_center <- rows[i]
      half_size <- size[i] / 2  # Half of the side length to adjust for centering
      
      #skip plotting diagonal if it does not matter (eg, in correlation matrices)
      if(tlbr_diag){
        on_diag <- cols[i] == (nr - rows[i] + 1)
      } else {
        on_diag <- cols[i] == rows[i]
      }
      
      if(!diag_matters){
        if(on_diag){
          next()
        }
      }
      
      # Draw rectangles
      
      # Draw rectangles
      if(plot_rects){
        rect(x_center - half_size + disp_x_by, y_center - half_size, 
             x_center + half_size + disp_x_by, y_center + half_size, 
             col = colors[color_indices[i]], border = NA, xpd = NA)
      }
      
      # Draw numbers
      if(plot_numbers && !is.na(mat_cols[i])){
        raw_lab <- sprintf("%.2f", mat_cols[i])
        if(omit_leading_zero){
          lab <- sub("^(-?)0\\.", "\\1.", raw_lab)
        } else {
          lab <- raw_lab
        }
        
        # choose cex so text fits within 90% of the cell in both width and height
        base_cex <- cex_lab * 0.6
        cex_h <- (0.8 * size[i]) / strheight(lab, cex = 1)
        cex_w <- (0.8 * size[i]) / strwidth(lab, cex = 1)
        cex_use <- min(base_cex, cex_h, cex_w)
        
        if(!plot_rects){
          rect(x_center - half_size + disp_x_by, y_center - half_size, 
               x_center + half_size + disp_x_by, y_center + half_size, 
               col = "white", border = NA, xpd = NA)
        }
        
        if(cex_use > 0){
          text(x_center + disp_x_by, y_center,
               labels = lab,
               cex = cex_use,
               col = num_col_mat[i],
               xpd = NA,
               adj = c(0.5, 0.5))
        }
      }
      
    }
    
  }
  
  #top legend
  if(plot_legend){
    leg_y <- pusr[4] + h/(10+ifelse(plot_labels, 0, 10)) + leg_loc_scale[2] * h
    legend_bb <- add_continuous_legend(colors = colors, 
                          labels = breaks[-1] + diff(breaks)/2, 
                          x = 1 + leg_loc_scale[1] * w, 
                          y = leg_y, 
                          vertical = F, xpd = NA, positions = 0:100/100, left_below = F,
                          scientific = F, main = legend_title, 
                          w = diff(pusr[1:2])/2/nmats, cex.main = cex.leg, 
                          return_BB = T)
  }
  
  #plot title above legend
  if(!is.null(main)){
    main.h <- strheight(main, cex = cex.main)
    text(x = nc / 2, 
         y = ifelse(plot_legend, legend_bb[4], nr) + main.h / 2, 
         labels = main, xpd = NA, font = 2, pos = 3, cex = cex.main)  
  }
  
  if(demarcate_maps & nmats > 1){
    vdemarc_xlocs <- ((1:nmats-1) * (nr + space_for_sep)) - space_for_sep / 2 + 1/2
    segments(x0 = vdemarc_xlocs[-1], x1 = vdemarc_xlocs[-1],
             y0 = pusr[3] - h / 20,
             y1 = pusr[4]  + h / 20, lwd = lwd_for_lines * 1.5, xpd = NA)
    text(x = vdemarc_xlocs + (nr+1)/2, y = pusr[3] - h / 10, labels = mat_names, xpd = NA)
  }
  
  #draw a bounding box around the entire heatmap
  if (outer_border) {
    if (outer_border_brackets) {
      
      vert_off <- nr / 50
      lat_off  <- nc / 50
      lat_in   <- nc / 20
      
      # LEFT BRACKET
      segments(0.5 - lat_off, 
               0.5 - vert_off,
               0.5 - lat_off,
               nr + 0.5 + vert_off,
               xpd = NA, lwd = border.lwd)
      
      segments(0.5 - lat_off,
               nr + 0.5 + vert_off,
               0.5 + lat_in,
               nr + 0.5 + vert_off,
               xpd = NA, lwd = border.lwd)
      
      segments(0.5 - lat_off,
               0.5 - vert_off,
               0.5 + lat_in,
               0.5 - vert_off,
               xpd = NA, lwd = border.lwd)
      
      
      # RIGHT BRACKET
      segments(nc + 0.5 + lat_off,
               0.5 - vert_off,
               nc + 0.5 + lat_off,
               nr + 0.5 + vert_off,
               xpd = NA, lwd = border.lwd)
      
      segments(nc + 0.5 - lat_in,
               nr + 0.5 + vert_off,
               nc + 0.5 + lat_off,
               nr + 0.5 + vert_off,
               xpd = NA, lwd = border.lwd)
      
      segments(nc + 0.5 - lat_in,
               0.5 - vert_off,
               nc + 0.5 + lat_off,
               0.5 - vert_off,
               xpd = NA, lwd = border.lwd)
      
    } else {
      rect(0.5, 0.5, nc + 0.5, nr + 0.5, lwd = border.lwd)
    }
  }
  
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
    extra_xlab = NULL,
    xlab_col = "black",
    extra_xlab_col = "black",
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
       labels = xlab, pos = 1, col = xlab_col)
  text(x = dx(mean(range(hax_xlocs))), 
       y = dy(hax_ylocs - tickh * 2 - strheight("0", "user") * 1), 
       labels = extra_xlab, pos = 1, col = extra_xlab_col)
  
  #left labels (for names)
  if(ylabs){
    xlocs_ll <- dx(min(lbs) - diff(range(c(lbs, rbs))) / 5
    )
    ylocs_ll <- (1:n - n/2) * space_scale_ylabs + n / 2 * space_scale_ylabs - n/10 
    ylocs_ll <- dy(ylocs_ll - mean(ylocs_ll) + n/2 + ylab_disp * diff(par("usr")[3:4]))
    text(x = xlocs_ll, y = ylocs_ll, labels = paste0(n:1, ". ", names(d)), 
         xpd = NA, pos = 2, col = border_cols, cex = 0.8)
    
    #connecting segments -- first drawing the line over
    segments(x0 = dx(min(lbs)), x1 = dx(lbs), 
             y0 = dy(1:n), y1 = dy(1:n), xpd = NA, col = border_cols)
    #then using a bezier curve to smoothly connect to the label
    bezier_curve(x0 = dx(min(lbs)), x1 = xlocs_ll, 
                 y0 = dy(1:n), y1 = ylocs_ll, k = 1, col = border_cols)
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
