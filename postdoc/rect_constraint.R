
get_rect_info <- function(rectangles) {
  # Compute values for constraints
  properties <- list(
    ar1 = rectangles$h1 / rectangles$w1,
    ar3 = rectangles$h3 / rectangles$w3,
    total_width = max(rectangles$w2 + rectangles$w3 + prop_dist*rectangles$h2, rectangles$w1),
    total_height = rectangles$h1 + rectangles$h2 * (prop_dist * 2 + 1),
    ar2_factor = (rectangles$w2 / rectangles$h2) / mean_iar2f,
    top_vs_bottom_height = rectangles$h1 / rectangles$h2
  )
  return(properties)
}

## A helper function: negative log PDF of a truncated normal( mu, sd ) for x>lower
##  If x <= lower, we return a large penalty (effectively infinite).
##  If x > lower, use standard formula:
##       pdf(x) = dnorm(x, mu, sd) / [1 - pnorm(lower, mu, sd)]
##  => -log(pdf(x)) is what we'll add as penalty. 
truncnorm_neglogpdf <- function(x, mean=1.5, sd=0.1, lower=1) {
  if (x <= lower) {
    return(1e8)  # big penalty if below lower
  }
  denom <- 1 - pnorm(lower, mean, sd)
  dens  <- dnorm(x, mean, sd)
  # PDF for truncated normal (x>lower) 
  pdf  <- dens / denom
  return(-log(pdf + 1e-15))  # add small to avoid log(0)
}

place_rectangles_bfgs <- function(
    prop_dist, 
    total_width,
    total_height,
    ar1,  # aspect ratio for b1
    ar3,  # aspect ratio for b3
    ## Some starting guesses:
    start_w1 = 0.5,
    start_w2 = 0.4,
    start_h2 = 0.6,
    ## Weights controlling how strictly we enforce near "hard" constraints
    big_penalty = 1e6,
    ## Weights for the truncated normal penalty
    ratio_mean_b1b2_w = 1.5,
    ratio_mean_b1b2_h = 1.5,
    ratio_sd_b1b2_w   = 0.1,
    ratio_sd_b1b2_h   = 0.1,
    ratio_mean_b2b3_w = 1,
    ratio_sd_b2b3_w = 0.1,
    ratio_mean_b1_b2b3_w = 0.75,
    ratio_sd_b1_b2b3_w = 0.25,
    mean_iar2f = 0.5,
    sd_iar2f = 0.1,
    iar2f_weight = 1,
    mean_ntwists = 2,
    sd_ntwists = 1,
    two_stage_iar2f = F,
    print_penalties = F,
    nrep = 20,
    ncores = 8
)
{
  # We'll define x = c(w1, w2, h2) as parameters to optimize.
  
  obj_fun <- function(x) {
    w1 <- x[1]
    w2 <- x[2]
    h2 <- x[3]
    
    # If w1, w2, h2 are negative or zero, penalty out
    if (w1 <= 0 || w2 <= 0 || h2 <= 0) {
      return(1e12 + sum(pmax(-x, 0)))  # giant penalty for negative
    }
    
    # 1) Dimensions that follow from x
    h1 <- ar1 * w1            # b1
    w3 <- h2 / ar3            # b3
    # b2 has width w2, height h2
    
    penalties <- c()
    
    # 2) "Hard" constraints as big penalties if violated:
    #    a) total_width:  w2 + prop_dist*h2 + w3 == total_width
    #       We'll penalize squared deviation from target:
    curr_width <- max(w1, w2 + prop_dist*h2 + w3)
    total_width_penalty <- big_penalty * (curr_width - total_width)^2
    penalties <- c(penalties, total_width_penalty = total_width_penalty)
    
    curr_height <- h1 + prop_dist*h2 * 2 + h2
    total_height_penalty <- ifelse(curr_height > total_height, big_penalty, 1) * (curr_height - total_height)^2
    penalties <- c(penalties, total_height_penalty = total_height_penalty)
    
    #    b) b1 bigger than b2 in width: (w1 > w2)
    b1b2_width_penalty <- 0
    if (w1 <= w2) {
      b1b2_width_penalty <- big_penalty * (w2 - w1)^2
    }
    penalties <- c(penalties, b1b2_width_penalty = b1b2_width_penalty)
    
    #       b1 bigger than b2 in height: (h1 > h2)
    b1b2_height_penalty <- 0
    if (h1 <= h2) {
      b1b2_height_penalty <- big_penalty * (h2 - h1)^2
    }
    penalties <- c(penalties, b1b2_height_penalty = b1b2_height_penalty)
    
    #    c) b2 + b3 > b1 in width => (w2 + w3 - w1 > 0)
    slack1 <- w2 + w3 - w1
    b2b3b1_width_penalty <- 0
    if (slack1 <= 0) {
      b2b3b1_width_penalty <- big_penalty * ( -slack1 )^2
    }
    penalties <- c(penalties, b2b3b1_width_penalty = b2b3b1_width_penalty)
    
    # 3) "Softer" constraints: 
    #    a) width_ratio_b1b2 = w1/w2 >1, prefer ~1.5
    width_ratio_b1b2  <- w1 / w2
    width_ratio_b1b2_penalty <- truncnorm_neglogpdf(width_ratio_b1b2, mean=ratio_mean_b1b2_w, sd=ratio_sd_b1b2_w, lower=1)
    penalties <- c(penalties, width_ratio_b1b2_penalty = width_ratio_b1b2_penalty)
    
    #    b) height_ratio = h1/h2 = (ar1*w1)/h2 >1, prefer ~1.5
    height_ratio_b1b2 <- h1 / h2
    height_ratio_b1b2_penalty <- truncnorm_neglogpdf(height_ratio_b1b2, mean=ratio_mean_b1b2_h, sd=ratio_sd_b1b2_h, lower=1)
    penalties <- c(penalties, height_ratio_b1b2_penalty = height_ratio_b1b2_penalty)
    
    #ratio of widths of w2 and w3
    width_b2b3_ratio <- w2 / w3
    width_b2b3_ratio_penalty <- truncnorm_neglogpdf(width_b2b3_ratio, mean=ratio_mean_b2b3_w, sd=ratio_sd_b2b3_w, lower=0)
    penalties <- c(penalties, width_b2b3_ratio_penalty = width_b2b3_ratio_penalty)
    
    #ratio of top row to bottom row widths
    width_b1_b2b3_ratio <- w1 / (w2 + w3 + prop_dist * h2)
    width_b1_b2b3_ratio_penalty <- truncnorm_neglogpdf(width_b1_b2b3_ratio, mean=ratio_mean_b1_b2b3_w, sd=ratio_sd_b1_b2b3_w, lower=0)
    penalties <- c(penalties, width_b1_b2b3_ratio_penalty = width_b1_b2b3_ratio_penalty)
    
    #aspect ratio of b2 -- strong minimum at 1, weak minimum at subsequent multiples
    b2_iar <- w2 / h2
    mod_b2_iar <- b2_iar %% mean_iar2f
    # b2_iar_value <- ifelse(b2_iar < mean_iar2f, 
    #                        mean_iar2f - b2_iar, #how far below the min you are (don't care for closeness to 0)
    #                        min(mod_b2_iar, mean_iar2f - mod_b2_iar)) #how far away from a multiple you are
    # b2_iar_penalty <- -dnorm(b2_iar_value, mean=0, sd=sd_iar2f, log=T) * 
    #   ifelse(b2_iar < mean_iar2f, big_penalty * 5, iar2f_weight) #imposes strong constraint if obs rat < min rat
    #hm this discontinuity makes the optimization always gravitate to just below 0
    b2_iar_penalty <- -dnorm(min(mod_b2_iar, mean_iar2f - mod_b2_iar), mean=0, sd=sd_iar2f, log=T) * iar2f_weight
    penalties <- c(penalties, b2_iar_penalty = b2_iar_penalty)
    
    #big penalty for being shorter than 1 full twist
    b2_iar_sub0_penalty <- 0
    if (b2_iar <= mean_iar2f) {
      b2_iar_sub0_penalty <- big_penalty * (b2_iar - mean_iar2f)^2
    }
    
    #number of twists
    ntwists <- b2_iar / mean_iar2f
    ntwists_penalty <- truncnorm_neglogpdf(ntwists, mean=mean_ntwists, sd=sd_ntwists, lower=0)
    penalties <- c(penalties, ntwists_penalty = ntwists_penalty)
    
    #sum up all costs
    cost <- sum(penalties)
    if(print_penalties){
      print(penalties[which.max(penalties)])
      print(round(penalties, 2))
      # if(which.min(penalties) == length(penalties)){
      #   print(paste0("b2_iar = ", b2_iar))
      #   print(paste0("mod_b2_iar = ", mod_b2_iar))
      #   print(paste0("b2_iar_value = ", b2_iar_value))
      #   print(paste0("b2_iar_penalty = ", b2_iar_penalty))
      #   print(paste0("sd_iar2f = ", sd_iar2f))
      # }
    }
    
    return(cost)
  }
  
  # We'll call optimx with method="BFGS" (unconstrained).
  # The large penalties should steer the solution into feasibility.
  
  if(two_stage_iar2f){
    #using a continuation method for the ar multiple of b2
    #first optimize not caring about iar2f, then initialize to those values
    orig_iar2f_weight <- iar2f_weight
    iar2f_weight <- 0
    res <- do.call(rbind, parallel::mclapply(1:nrep, function(nri){
      if(nri > 1){
        start_par <- c(runif(1,0,1) * total_width, runif(1,0,1) * total_width, runif(1,0,1) * total_height)
      } else {
        start_par <- c(start_w1 * total_width, start_w2 * total_width, start_h2 * total_height)
      }
      optimx::optimx(
        par     = start_par,
        fn      = obj_fun,
        method  = "BFGS",
        control = list(
          maxit = 1000,
          trace = 0   # set >0 for more debugging
        )
      ) 
    }, mc.cores = ncores))
    
    best_row <- which.min(res$value)
    w1_opt <- res[best_row, "p1"]
    w2_opt <- res[best_row, "p2"]
    h2_opt <- res[best_row, "p3"]
    start_par <- c(w1_opt, w2_opt, h2_opt)
    iar2f_weight <- orig_iar2f_weight
  } else {
    start_par <- c(start_w1 * total_width, start_w2 * total_width, start_h2 * total_height)
  }
  
  #perform the final optimization
  res <- optimx::optimx(
    par     = start_par,
    fn      = obj_fun,
    method  = "BFGS",
    control = list(
      maxit = 1000,
      trace = 0
    )
  )
  
  # Extract best solution from optimx result (could have multiple rows for method combos)
  best_row <- which.min(res$value)
  w1_opt <- res[best_row, "p1"]
  w2_opt <- res[best_row, "p2"]
  h2_opt <- res[best_row, "p3"]
  
  # Reconstruct final geometry
  h1_opt <- ar1 * w1_opt
  w3_opt <- h2_opt / ar3
  h3_opt <- h2_opt
  
  # Coordinates: 
  # b1 top-left: (0,0)
  # b2 top-left: (0, h1_opt + prop_dist*h2_opt)
  # b3 top-left: (w2_opt + prop_dist*h2_opt, same y as b2)
  x1 <- 0
  y1 <- 0
  x2 <- 0
  y2 <- h1_opt + prop_dist * h2_opt
  x3 <- w2_opt + prop_dist * h2_opt
  y3 <- y2
  
  # Return a list with geometry & final objective
  out <- list(
    objective  = res[best_row, "value"],
    w1         = w1_opt,
    w2         = w2_opt,
    h2         = h2_opt,
    h1         = h1_opt,
    w3         = w3_opt,
    h3         = h3_opt,
    b1         = list(x = x1, y = y1, width = w1_opt, height = h1_opt),
    b2         = list(x = x2, y = y2, width = w2_opt, height = h2_opt),
    b3         = list(x = x3, y = y3, width = w3_opt, height = h3_opt),
    optimx_out = res
  )
  return(out)
}

plot_rectangles <- function(rectangles){
  
  curr_mar <- par("mar")
  par(mar = c(4,4,2,2))
  
  plot(
    0, 0, type = "n", xlim = c(0, max(rectangles$b3$x + rectangles$b3$width, 
                                      rectangles$b1$x + rectangles$b1$width)), 
    ylim = c(-max(rectangles$b2$y + rectangles$b2$height), 0),
    xlab = "X", ylab = "Y", asp = 1,
    main = "Optimized Rectangles"
  )
  
  # Define colors
  colors <- c("red", "blue", "green")
  
  # Draw each rectangle
  rect(
    rectangles$b1$x, rectangles$b1$y, 
    rectangles$b1$x + rectangles$b1$width, -(rectangles$b1$y + rectangles$b1$height), 
    col = colors[1], border = "black", xpd = NA
  )
  text(rectangles$b1$x + rectangles$b1$width/2, 
       -(rectangles$b1$y + rectangles$b1$height/2), labels = "b1", cex = 2, col = 0)
  
  rect(
    rectangles$b2$x, -rectangles$b2$y, 
    rectangles$b2$x + rectangles$b2$width, -(rectangles$b2$y + rectangles$b2$height), 
    col = colors[2], border = "black", xpd = NA
  )
  text(rectangles$b2$x + rectangles$b2$width/2, 
       -(rectangles$b2$y + rectangles$b2$height/2), labels = "b2", cex = 2, col = 0)
  
  rect(
    rectangles$b3$x, -rectangles$b3$y, 
    rectangles$b3$x + rectangles$b3$width, -(rectangles$b3$y + rectangles$b3$height), 
    col = colors[3], border = "black", xpd = NA
  )
  text(rectangles$b3$x + rectangles$b3$width/2, 
       -(rectangles$b3$y + rectangles$b3$height/2), labels = "b3", cex = 2, col = 0)
  
  par(mar = curr_mar)

  # Example of using the function and constructing the output messages
  rect_info <- get_rect_info(rectangles)
  
  # Print outputs
  print(paste0("ar1: ", round(ar1, 3), " vs. ", round(rect_info$ar1, 3)))
  print(paste0("ar3: ", round(ar3, 3), " vs. ", round(rect_info$ar3, 3)))
  print(paste0("total width: ", round(total_width, 3), " vs. ", round(rect_info$total_width, 3)))
  print(paste0("total height: ", round(total_height, 3), " vs. ", round(rect_info$total_height, 3)))
  print(paste0("ar2 factor: ", round(rect_info$ar2_factor, 3)))
  print(paste0("top vs. bottom height: ", round(rect_info$top_vs_bottom_height, 3)))
  
}

#### example ####
run_example <- F
if(run_example){
  prop_dist   <- 0.2
  total_width <- 1
  total_height <- 0.5
  ar1         <- 0.2   # fixed aspect ratio for b1 (h / w)
  ar3         <- 0.3  # fixed aspect ratio for b3 (h / w)
  mean_iar2f <- 2 # inverse aspect ratio for b2 (w / h)
  
  rectangles <- place_rectangles_bfgs(
    prop_dist   = prop_dist,
    total_width = total_width,
    total_height = total_height,
    ar1         = ar1,
    ar3         = ar3,
    start_w1    = 0.5,
    start_w2    = 0.4,
    start_h2    = 0.5,
    big_penalty = 1e2,
    ratio_mean_b1b2_w  = 1.55,
    ratio_sd_b1b2_w    = 1,
    ratio_mean_b1b2_h  = 2,
    ratio_sd_b1b2_h    = 0.25,
    ratio_mean_b2b3_w = 1,
    ratio_sd_b2b3_w = 0.5,
    ratio_mean_b1_b2b3_w = 0.85,
    ratio_sd_b1_b2b3_w = 0.1,
    mean_iar2f = mean_iar2f,
    sd_iar2f = 0.2,
    iar2f_weight = 1,
    two_stage_iar2f = T, 
    nrep = 5
  )
  
  plot_rectangles(rectangles)
}


