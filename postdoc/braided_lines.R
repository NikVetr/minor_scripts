#' @title Draw a Braided Line Segment
#'
#' @description This function draws a line segment between two points, (x0, y0)
#'   and (x1, y1), as a "braid" of multiple colored strands. It ensures the
#'   braid has a consistent *perceived width* regardless of the segment's
#'   direction or the plot's aspect ratio.
#'
#' @param x0 The x-coordinate of the starting point.
#' @param y0 The y-coordinate of the starting point.
#' @param x1 The x-coordinate of the ending point.
#' @param y1 The y-coordinate of the ending point.
#' @param cols A character vector of colors for the strands. The number of
#'   strands is `length(cols)`.
#' @param bwidth The total width of the braid in **points** (1/72.27 inches).
#'   This is similar to the `lwd` parameter in `par`. A value of 10-20 is
#'   good for visibility.
#' @param twists The total number of full twists (cycles) the braid should
#'   make along the length of the segment.
#' @param min_segments The minimum number of polygon segments used to
#'   approximate the curved braid. The function will *add* segments at
#'   all crossover points, so the true number of segments will be
#'   `min_segments` + `n_crossovers`.
#' @param ... Additional arguments passed to the `polygon()` function for
#'   drawing each segment (e.g., `border`, `lty`). `border = NA` is the
#'   default to ensure a clean occlusion effect.
#'
#' @details
#' This function is an advanced version that calculates the *exact*
#' crossover points for all strands. It then merges these "topological"
#' breaks with a "smoothness" grid (`min_segments`) to create the
#' minimum number of polygon segments required to be both smooth
#' and topologically correct.
#'
#' @return Returns `NULL` invisibly. This function is called for its
#'   side effect of drawing on the current plot device.
#'
#' @export
braid_segment <- function(x0, y0, x1, y1, cols, bwidth = 10, twists = 2, min_segments = 50, ...) {
  
  # bwidth is the total width of the braid in *points* (1/72.27 inches)
  
  # --- 0. Convert bwidth from points to inches ---
  width_inches <- bwidth / 72.27
  
  # --- 1. Get aspect-ratio-corrected perpendicular vector ---
  
  # Get plot dimensions
  usr <- par("usr")
  pin <- par("pin")
  
  # Handle zero-size plot dimensions
  if (pin[1] == 0 || pin[2] == 0 || (usr[2] - usr[1]) == 0 || (usr[4] - usr[3]) == 0) {
    x_scale <- 1
    y_scale <- 1
    warning("Plot device not set up or has zero dimensions. Aspect ratio correction may fail.")
  } else {
    # Calculate inches per data unit
    x_scale <- pin[1] / (usr[2] - usr[1]) # inches / x-unit
    y_scale <- pin[2] / (usr[4] - usr[3]) # inches / y-unit
  }
  
  # Vector of the segment in data units
  V_data_x <- x1 - x0
  V_data_y <- y1 - y0
  
  # Length of segment in data units
  L_data <- sqrt(V_data_x^2 + V_data_y^2)
  
  # Handle zero-length segment
  if (L_data == 0) return(invisible(NULL))
  
  # Normalized segment vector (unit vector)
  U_data_x <- V_data_x / L_data
  U_data_y <- V_data_y / L_data
  
  # --- Calculate the aspect-corrected perpendicular vector ---
  P_inch_x <- -V_data_y * y_scale
  P_inch_y <- V_data_x * x_scale
  L_inch_p <- sqrt(P_inch_x^2 + P_inch_y^2)
  
  if (L_inch_p == 0) {
    P_data_x_unit_inch <- -U_data_y
    P_data_y_unit_inch <- U_data_x
  } else {
    P_norm_inch_x <- P_inch_x / L_inch_p
    P_norm_inch_y <- P_inch_y / L_inch_p
    P_data_x_unit_inch <- P_norm_inch_x / x_scale
    P_data_y_unit_inch <- P_norm_inch_y / y_scale
  }
  
  # --- 2. Set up braid parameters ---
  n_strands <- length(cols)
  if (n_strands == 0) return(invisible(NULL))
  
  Width_vec_x <- width_inches * P_data_x_unit_inch
  Width_vec_y <- width_inches * P_data_y_unit_inch
  
  strand_width_vec_x <- Width_vec_x / n_strands
  strand_width_vec_y <- Width_vec_y / n_strands
  
  amplitude_scalar <- (n_strands - 1) / (2 * n_strands)
  if (n_strands <= 1) amplitude_scalar <- 0
  
  Amp_vec_x <- amplitude_scalar * Width_vec_x
  Amp_vec_y <- amplitude_scalar * Width_vec_y
  
  phases <- (seq_len(n_strands) - 1) / n_strands * 2 * pi
  
  # Frequency (cycles per data unit length)
  frequency <- twists / L_data
  
  # --- 3. Generate segment breaks ---
  
  # 3a. Generate base breaks for smoothness
  t_breaks_smooth <- seq(0, L_data, length.out = min_segments + 1)
  
  # 3b. Calculate all topological crossover points
  t_breaks_crossover <- c()
  if (n_strands > 1 && twists > 0) {
    for (i in 1:(n_strands - 1)) {
      for (j in (i + 1):n_strands) {
        # Find t values where cos(A) = cos(B)
        # A = 2*pi*f*t + phase_i, B = 2*pi*f*t + phase_j
        # Solved by A = -B + 2*pi*k
        # t = (2*pi*k - phase_i - phase_j) / (4 * pi * f)
        
        phase_sum <- phases[i] + phases[j]
        
        # Find range of k to check
        # k_min_val = (4*pi*f*0 + phase_sum) / (2*pi)
        # k_max_val = (4*pi*f*L_data + phase_sum) / (2*pi)
        k_min_val <- (phase_sum) / (2 * pi)
        k_max_val <- (4 * pi * frequency * L_data + phase_sum) / (2 * pi)
        
        k_range <- floor(k_min_val):ceiling(k_max_val)
        
        t_cross <- (2 * pi * k_range - phase_sum) / (4 * pi * frequency)
        
        # Keep only those within the segment (and not endpoints)
        t_cross_valid <- t_cross[t_cross > 1e-6 & t_cross < (L_data - 1e-6)]
        t_breaks_crossover <- c(t_breaks_crossover, t_cross_valid)
      }
    }
  }
  
  # 3c. Combine all breaks and sort
  all_breaks <- sort(unique(c(t_breaks_smooth, t_breaks_crossover, 0, L_data)))
  n_seg <- length(all_breaks) - 1
  
  # --- 4. Generate and draw polygons ---
  
  # Capture ... args, and ensure border=NA
  args <- list(...)
  if (is.null(args$border)) {
    args$border <- NA
  }
  
  # Loop over each (now minimal) segment
  for (j in seq_len(n_seg)) {
    t0_calc <- all_breaks[j]
    t1_calc <- all_breaks[j+1]
    
    # Add a tiny overlap to t1 to hide anti-aliasing seams
    overlap <- (t1_calc - t0_calc) * 0.001 # 0.1% overlap
    t0 <- t0_calc
    t1 <- t1_calc + overlap
    
    # Z-order is constant within this segment, so we check midpoint
    t_mid <- (t0 + t1) / 2
    
    # --- Calculate Z-order at the segment midpoint ---
    angle_mid <- 2 * pi * frequency * t_mid + phases
    z_order <- cos(angle_mid)
    draw_order <- order(z_order, decreasing = FALSE) # Draw from back to front
    
    # --- Draw the polygon for each strand in z-order ---
    for (i in draw_order) {
      
      # --- Calculate 4 vertices of the polygon segment ---
      angle0 <- 2 * pi * frequency * t0 + phases[i]
      angle1 <- 2 * pi * frequency * t1 + phases[i]
      
      offset0_scalar <- sin(angle0)
      offset1_scalar <- sin(angle1)
      
      C0_x <- x0 + t0 * U_data_x + offset0_scalar * Amp_vec_x
      C0_y <- y0 + t0 * U_data_y + offset0_scalar * Amp_vec_y
      
      C1_x <- x0 + t1 * U_data_x + offset1_scalar * Amp_vec_x
      C1_y <- y0 + t1 * U_data_y + offset1_scalar * Amp_vec_y
      
      edge_vec_x <- strand_width_vec_x / 2
      edge_vec_y <- strand_width_vec_y / 2
      
      # 4 vertices:
      P1_x <- C0_x - edge_vec_x
      P1_y <- C0_y - edge_vec_y
      
      P2_x <- C1_x - edge_vec_x
      P2_y <- C1_y - edge_vec_y
      
      P3_x <- C1_x + edge_vec_x
      P3_y <- C1_y + edge_vec_y
      
      P4_x <- C0_x + edge_vec_x
      P4_y <- C0_y + edge_vec_y
      
      # Prepare arguments for polygon()
      poly_args <- c(
        list(x = c(P1_x, P2_x, P3_x, P4_x), 
             y = c(P1_y, P2_y, P3_y, P4_y), 
             col = cols[i]),
        args
      )
      
      # Draw the polygon
      do.call("polygon", poly_args)
    }
  }
  
  return(invisible(NULL))
}

# --- Example script for using braid_segment() ---
#
# 1. Make sure the `braid_segment.R` file is in your working directory
#    or source it first:
#    source("braid_segment.R")
#
# 2. This test function will create a 3-panel plot to show that the
#    braid width is consistent across different plot aspect ratios.
#

run_braid_test <- function() {
  
  # Define some color palettes
  cols2 <- c("#E41A1C", "#377EB8") # Red, Blue
  cols3 <- c("#E41A1C", "#377EB8", "#4DAF4A") # Red, Blue, Green
  cols5 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00") # 5 Colors
  
  # Use a wide bwidth for high visibility in the example
  bw <- 20 # 20 points
  
  # Set a base smoothness. The function will add more segments
  # as needed at crossover points.
  min_seg <- 50
  
  # Save old plot parameters and set up a 1x3 layout
  old_par <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  # --- Test 1: Stretched plot window (asp = 0.5) ---
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
       asp = 0.5, main = "Stretched (asp=0.5)",
       xlab = "X", ylab = "Y")
  box()
  # Diagonal line
  braid_segment(1, 1, 5, 5, cols = cols2, bwidth = bw, twists = 3, min_segments = min_seg)
  # Vertical line
  braid_segment(7, 2, 7, 8, cols = cols3, bwidth = bw, twists = 4, min_segments = min_seg)
  # Horizontal line
  braid_segment(2, 8, 8, 8, cols = cols2, bwidth = bw, twists = 3, min_segments = min_seg)
  # Other diagonal
  braid_segment(1, 9, 9, 1, cols = cols5, bwidth = bw, twists = 5, min_segments = min_seg)
  
  # --- Test 2: Square plot window (asp = 1) ---
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
       asp = 1, main = "Square (asp=1)",
       xlab = "X", ylab = "Y")
  box()
  braid_segment(1, 1, 5, 5, cols = cols2, bwidth = bw, twists = 3, min_segments = min_seg)
  braid_segment(7, 2, 7, 8, cols = cols3, bwidth = bw, twists = 4, min_segments = min_seg)
  braid_segment(2, 8, 8, 8, cols = cols2, bwidth = bw, twists = 3, min_segments = min_seg)
  braid_segment(1, 9, 9, 1, cols = cols5, bwidth = bw, twists = 5, min_segments = min_seg)
  
  # --- Test 3: Squished plot window (asp = 2) ---
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
       asp = 2, main = "Squished (asp=2)",
       xlab = "X", ylab = "Y")
  box()
  braid_segment(1, 1, 5, 5, cols = cols2, bwidth = bw, twists = 3, min_segments = min_seg)
  braid_segment(7, 2, 7, 8, cols = cols3, bwidth = bw, twists = 4, min_segments = min_seg)
  braid_segment(2, 8, 8, 8, cols = cols2, bwidth = bw, twists = 3, min_segments = min_seg)
  braid_segment(1, 9, 9, 1, cols = cols5, bwidth = bw, twists = 5, min_segments = min_seg)
  
  mtext("Braid Segment Test (Perceived width should be constant)", 
        outer = TRUE, cex = 1.5)
  
  # Restore old plot parameters
  par(old_par)
}

# --- To Run This Example ---
# 1. Save the `braid_segment` function from `braid_segment.R`
# 2. Make sure it's in your R environment (e.g., by running `source("braid_segment.R")`)
# 3. Run the test function:
#
   run_braid_test()
#
#    (You may need to resize your plot window to see all three panels clearly)