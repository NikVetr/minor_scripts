library(ipumsr)
library(data.table)
library(parallel)

#### plotting functions ####

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


#### data processing ####

if(!exists("dsub")){
  #read in the data
  ddi <- read_ipums_ddi("~/Downloads/usa_00004.xml")
  d <- read_ipums_micro(ddi)
  
  #subset to married couples with spouse present in household
  dsub <- d[d$MARST == 1, c("SERIAL", "PERNUM", "SPLOC", "SEX", "AGE", "BIRTHYR", "YRMARR")]
  
  #subsample d to write code more efficiently
  subsamp <- 1:1E4
  ds <- dsub[subsamp,]
  rm(d); gc()
}

#assign data to work with (for easy switching)
dat <- ds
dat <- data.table(dat)
dat$serial <- dat$SERIAL

#extract households and ID spouse pairs (splitting households with >2 couples)
# h <- split(dat, dat$SERIAL)
# n_parts <- 12
# n_cores <- 12
dat$household_index <- floor(rank(dat$SERIAL))
# part_breaks <- round(quantile(dat$household_index, prob = seq(0, 1, length.out = n_parts + 1)))[-(n_parts+1)]
# dat$part <- Reduce("+", lapply(part_breaks, function(bi) dat$household_index >= bi))
h <- dat[, .(list(.SD)), by = SERIAL]$V1
# h <- unlist(mclapply(h, function(x) {x[, .(list(.SD)), by = household_index]$V1}, mc.cores = n_cores), recursive = F)

if(nrow(dat) != nrow(dsub)){h <- h[-length(h)]}
n_per_household <- sapply(h, nrow)
hc <- h[n_per_household == 2]
hc_gr2 <- h[n_per_household > 2]
hc_gr2 <- unlist(mclapply(hc_gr2, function(x){
  spouse_IDs <- x[,c("PERNUM", "SPLOC")]
  spouse_inds <- cbind(1:nrow(x), match(spouse_IDs$SPLOC, spouse_IDs$PERNUM))
  ind_pairs <- t(apply(spouse_inds, 1, sort))
  uniq_ind_pairs <- unique(ind_pairs)
  out <- lapply(1:nrow(uniq_ind_pairs), function(i){
    sub_out <- x[uniq_ind_pairs[i,],]
    sub_out$household_index <- sub_out$household_index + i / 10^(floor(log10(i))+1)
    return(sub_out)
  })
  
  return(out)
}, mc.cores = 12), recursive = F)

#combine back into h and reorder into new index
h <- do.call(rbind, c(hc, hc_gr2))
rm(list = c("hc", "hc_gr2"))
h <- h[order(h$household_index),]
ci_ranks <- rank(h$household_index, ties.method = "min")
dci_ranks <- diff(ci_ranks)
if(all(names(table(dci_ranks)) %in% c(0,2))){
  dci_ranks[dci_ranks==2] <- 1
  ci_ranks <- cumsum(c(1, dci_ranks))
  h$couple_index <- ci_ranks
} else {
  stop("something is wrong.")
}
max_index <- max(h$couple_index)

#confirm we did this right
if(!all((1:nrow(h) - h$couple_index*2) %in% c(0,-1))){
  stop("something is wrong.")
}
spot_checks <- sample(unique(h$couple_index), 100)
if(!all(sapply(spot_checks, function(ci){
  his <- h[ci*2 - 0:1,]$household_index
  his[1] == his[2]
}))){
  stop("something is wrong.")
}

#process pairs of households into a single household
mcprint <- function(...){
  system(sprintf('printf "%s"', paste0(..., collapse="")))
}

#this is too slow on the full dataset, even with mclapply (also bad memory handling)
# m <- data.table(do.call(rbind, mclapply(1:max_index, function(i){
#   
#   #report progress
#   if(i %% 1000 == 0){
#     print_message <- paste0("(", round(i / max_index * 100, 1), ") ")
#     mcprint(print_message)  
#   }
#   x <- h[ci*2 - 0:1,]
# 
#   #check that these pairs married each other that year
#   if(x$YRMARR[1] != x$YRMARR[2]){
#     return(NULL)
#   }
#   #remove same-sex spouses
#   if(x$SEX[1] == x$SEX[2]){
#     return(NULL)
#   }
#   #remove pairs with unknown sex
#   if(any(x$SEX == 9)){
#     return(NULL)
#   }
# 
# 
#   male <- x[x$SEX == 1,]
#   female <- x[x$SEX == 2,]
#   out <- c(age_m = as.numeric(male$AGE),
#            age_f = as.numeric(female$AGE),
#            byear_m = as.numeric(male$BIRTHYR),
#            byear_f = as.numeric(female$BIRTHYR),
#            year_married = as.numeric(x$YRMARR[1]))
#   return(out)
#   
#   # male <- h[ci*2 - 0:1,][c(h[ci*2 - 0:1,"SEX"] == 1),]
#   # female <- h[ci*2 - 0:1,][c(h[ci*2 - 0:1,"SEX"] == 2),]
#   # out <- c(age_m = as.numeric(male$AGE),
#   #          age_f = as.numeric(female$AGE),
#   #          byear_m = as.numeric(male$BIRTHYR),
#   #          byear_f = as.numeric(female$BIRTHYR),
#   #          year_married = as.numeric(x$YRMARR[1]))
#   # return(out)
#   
# }, mc.cores = 12)))

#match couples with join
h$sex <- h$SEX
sex_mismatch <- which(!(h$sex %in% c(1,2)))
if(length(sex_mismatch) != 0){
  smm_cis <- h$couple_index[sex_mismatch]
  smm_i <- unlist(lapply(smm_cis*2, function(ci) ci - 0:1))
  h <- h[-smm_i,]
}
hxs <- h[, .(list(.SD)), by = sex]$V1
couples <- merge(hxs[[1]], hxs[[2]], by = "couple_index")

#check for errors and remove entries accordingly
couples <- couples[couples$YRMARR.x == couples$YRMARR.y,]
couples <- couples[couples$PERNUM.x == couples$SPLOC.y,]
couples <- couples[couples$PERNUM.y == couples$SPLOC.x,]
couples <- couples[couples$serial.x == couples$serial.y,]

married <- data.table(age_m = as.numeric(couples$AGE.x),
                age_f = as.numeric(couples$AGE.y),
                byear_m = as.numeric(couples$BIRTHYR.x),
                byear_f = as.numeric(couples$BIRTHYR.y),
                year_married = as.numeric(couples$YRMARR.x))

#compute some other helpful variables
married$diff_age_mvf <- married$age_m - married$age_f
married$diff_age_fvm <- married$age_f - married$age_m
married$age_m_at_marriage <- married$year_married - married$byear_m
married$age_f_at_marriage <- married$year_married - married$byear_f
married$mean_age_at_marriage <- (married$age_m_at_marriage + married$age_f_at_marriage) / 2

#### plot preprocessing ####

#reset variable
m <- married

#set plot parameters
use_f <- T
helper_plots <- unlist(lapply(!use_f, function(use_f){
  
  if(use_f){
    var_to_use <- "age_f_at_marriage"
    m$diff_age <- m$diff_age_fvm
    focal_spouse <- "Wife"
    nonfocal_spouse <- "Husband"
  } else {
    var_to_use <- "age_m_at_marriage"  
    m$diff_age <- m$diff_age_mvf
    focal_spouse <- "Husband"
    nonfocal_spouse <- "Wife"
  }
  
  # bin the ages
  age_diff_quantiles <- quantile(m$diff_age, probs = c(1E-3, 1-1E-3))
  m_diff_range <- max(abs(age_diff_quantiles)) * c(-1,1)
  xlabs <- c(rev(seq(from = 0, to = m_diff_range[1], by = -5)), 
             seq(from = 0, to = m_diff_range[2], by = 5)[-1])
  m_bin_all <- m[, .(split_list = list(.SD)), by = var_to_use]
  m_bin <- m_bin_all$split_list
  m_bin_ages <- as.numeric(unlist(m_bin_all[,..var_to_use]))
  m_bin <- m_bin[order(m_bin_ages)]
  m_bin_ages <- sort(m_bin_ages)
  
  #filter out low samples
  min_sample_size <- 20
  to_exclude <- sapply(m_bin, nrow) < min_sample_size
  m_bin <- m_bin[!to_exclude]
  m_bin_ages <- m_bin_ages[!to_exclude]
  m_age_range <- range(m_bin_ages) + c(-1,1)
  
  #compute things to plot
  m_bin_mean_diffs <- sapply(m_bin, function(x) mean(x$diff_age))
  m_bin_SD_diffs <- sapply(m_bin, function(x) sd(x$diff_age))
  m_bin_prop_diffs <- sapply(m_bin, function(x) mean(x$diff_age > 0))
  n_bin_prop_total <- sapply(m_bin, function(x) nrow(x))
  
  return(list(m_bin_ages = m_bin_ages,
              m_bin_mean_diffs = m_bin_mean_diffs,
              m_bin_SD_diffs = m_bin_SD_diffs,
              m_bin_prop_diffs = m_bin_prop_diffs,
              n_bin_prop_total = n_bin_prop_total))
}), recursive = F)

#set actual plot parameters
if(use_f){
  var_to_use <- "age_f_at_marriage"
  m$diff_age <- m$diff_age_fvm
  focal_spouse <- "Wife"
  nonfocal_spouse <- "Husband"
} else {
  var_to_use <- "age_m_at_marriage"  
  m$diff_age <- m$diff_age_mvf
  focal_spouse <- "Husband"
  nonfocal_spouse <- "Wife"
}

# bin the ages
age_diff_quantiles <- quantile(m$diff_age, probs = c(1E-3, 1-1E-3))
m_diff_range <- max(abs(age_diff_quantiles)) * c(-1,1)
xlabs <- c(rev(seq(from = 0, to = m_diff_range[1], by = -5)), 
           seq(from = 0, to = m_diff_range[2], by = 5)[-1])
m_bin_all <- m[, .(split_list = list(.SD)), by = var_to_use]
m_bin <- m_bin_all$split_list
m_bin_ages <- as.numeric(unlist(m_bin_all[,..var_to_use]))
m_bin <- m_bin[order(m_bin_ages)]
m_bin_ages <- sort(m_bin_ages)

#filter out low samples
min_sample_size <- 20
to_exclude <- sapply(m_bin, nrow) < min_sample_size
m_bin <- m_bin[!to_exclude]
m_bin_ages <- m_bin_ages[!to_exclude]
m_age_range <- range(m_bin_ages) + c(-1,1)

#compute things to plot
m_bin_mean_diffs <- sapply(m_bin, function(x) mean(x$diff_age))
m_bin_SD_diffs <- sapply(m_bin, function(x) sd(x$diff_age))
m_bin_prop_diffs <- sapply(m_bin, function(x) mean(x$diff_age > 0))
n_bin_prop_total <- sapply(m_bin, function(x) nrow(x))

#get colors for plot
colgrad <- colorRampPalette(c("brown2", "darkviolet", "deepskyblue4"))(diff(range(m_bin_mean_diffs)) + 1)
if(use_f){
  colgrad <- rev(colgrad)
}
col_seq_vals <- seq(min(round(m_bin_mean_diffs)), max(round(m_bin_mean_diffs)), by = 1)
bin_col_inds <- round(m_bin_mean_diffs) - min(round(m_bin_mean_diffs)) + 1
bin_cols <- colgrad[bin_col_inds]

#obtain the KDEs
ndens <- 256
kern_adj <- 2
m_dens <- lapply(m_bin, function(x){
  if(nrow(x) < min_sample_size) return(NULL) #should be redundant
  #for smoothing, perturb this with a difference of two uniforms
  diff_age <- x$diff_age + runif(length(x$diff_age), 0, 1) - runif(length(x$diff_age), 0, 1)
  out <- density(diff_age, from = min(m$diff_age), to = max(m$diff_age), n = ndens, adjust = kern_adj)
  return(data.frame(x = out$x, y = out$y))
})
scale_vertical <- 16
dens_thres <- 1/ndens/4
mass_thresh <- 0.999
mass_thresh_bounds <- c((1-mass_thresh)/2, mass_thresh+(1-mass_thresh)/2)
hard_bounds_range <- c(-26,35)

#### plotting ###
filename <- paste0("marriage-ages_", tolower(focal_spouse), "-age-at-marriage")
grDevices::cairo_pdf(filename = paste0("~/", filename,".pdf"), 
                     width = 850 / 72, height = 850 / 72)


par(xpd = NA, mar = c(6,6,3,19), fig=c(0,1,0,1))
plot.new()
plot.window(xlim = m_diff_range, ylim = m_age_range)
axis(2, at = seq(from = m_age_range[1], to = m_age_range[2], by = 5))
axis(1, at = xlabs, labels = xlabs)

for(i in length(m_dens):1){
  if(is.null(m_dens[[i]])){
    next()
  }
  ad <- m_dens[[i]]
  
  #subset to non-zero mass regions
  ad_cumul_mass <- cumsum((ad$y[-1] + ad$y[-nrow(ad)]) / 2 * diff(ad$x))
  ad_subset <- ad_cumul_mass > mass_thresh_bounds[1] & ad_cumul_mass < mass_thresh_bounds[2]
  ad_subset <- c(ad_subset[1], ad_subset)
  ad_subset <- ad_subset & ad$x > hard_bounds_range[1] & ad$x < hard_bounds_range[2]
  # ad_subset <- ad_subset & ad$y > dens_thres
  
  ad_subset_dens_thresh <- rep(F, nrow(ad))
  ad_subset_dens_thresh[min(which(ad$y > dens_thres)) : max(which(ad$y > dens_thres))] <- T
  ad_subset <- ad_subset & ad_subset_dens_thresh
  
  ad <- ad[ad_subset,]
  polygon(x = c(ad$x, rev(ad$x)),
          y = m_bin_ages[i] + c(ad$y, rep(0, nrow(ad))) * scale_vertical, 
          border = 1, col = bin_cols[i]
  )
}

#annotate

segments(x0 = 0, x1 = 0, y0 = par("usr")[3], y1 =  m_age_range[2], lty = 2, col = adjustcolor(1, 0.5), lwd = 2)
text(x = 2, y = par("usr")[4] + diff(par("usr")[3:4])/80, pos = 3, 
     labels = paste0("Age Differences Between ", gsub("f", "v", focal_spouse), 
                     "s and ", gsub("f", "v", nonfocal_spouse), "s at Marriage"), 
     cex = 2,
     font = 1, family = "Times")
text(x = mean(xlabs[xlabs<0]), 
     y = m_age_range[2], pos = 3, 
     labels = paste0(nonfocal_spouse, " is Older"), cex = 2.5,
     font = 2, col = colgrad[1])
text(x = mean(xlabs[xlabs>0]) * 1.2, 
     y = m_age_range[2], pos = 3, 
     labels = paste0(focal_spouse, " is Older"), cex = 2.5,
     font = 2, col = colgrad[length(colgrad)])
text(x = 0, y = par("usr")[3] - diff(par("usr")[3:4])/20, cex = 1.5, pos = 1,
     labels = paste0("Difference in Age between ", focal_spouse, " and ", nonfocal_spouse, " (Years)"))
text(x = par("usr")[1] - diff(par("usr")[1:2])/10, y = par("usr")[3] + diff(par("usr")[3:4])/2, cex = 1.5, pos = 3,
     labels = paste0("Age of ", focal_spouse, " at Marriage (Years)"), srt = 90)

add_continuous_legend(rev(colgrad), labels = col_seq_vals, 
                      positions = (col_seq_vals - min(col_seq_vals)) / diff(range(col_seq_vals)), 
                      x = m_diff_range[2] - diff(m_diff_range[1:2])/15, y = 34, h = 20, left_below = F, main = "Mean Age\nDifference")


#add four helper plots
helper_x_prop <- c(0.825, 0.98)
helper_y_prop <- c(0.15, 0.3)
par(fig=c(helper_x_prop, helper_y_prop), mar = c(0,0,0,0), new=TRUE)
plot(m_bin_ages, m_bin_mean_diffs, type = "l", xlab = paste0(focal_spouse, " Age"), ylab = "Mean(Difference)", lwd = 2)
lines(helper_plots$m_bin_ages, helper_plots$m_bin_mean_diffs, 
      col = adjustcolor(1,0.25), lwd = 2, xpd = F)

abline(h=0, xpd = T, lty = 2, col = adjustcolor(1, 0.5))
abline(v=1:4*20, col = adjustcolor(1, 0.3), xpd = T)

helper_y_prop <- helper_y_prop + 0.16
par(fig=c(helper_x_prop, helper_y_prop), mar = c(0,0,0,0), new=TRUE)
plot(m_bin_ages, m_bin_prop_diffs, 
     ylim = c(0, 1),
     type = "l", xlab = "", ylab = "Pr(Difference > 0)", lwd = 2, xaxt = "n")
lines(helper_plots$m_bin_ages, helper_plots$m_bin_prop_diffs, 
      col = adjustcolor(1,0.25), lwd = 2, xpd = F)
abline(h=0.5, xpd = T, lty = 2, col = adjustcolor(1, 0.5))
abline(v=1:4*20, col = adjustcolor(1, 0.3), xpd = T)

helper_y_prop <- helper_y_prop + 0.16
par(fig=c(helper_x_prop, helper_y_prop), mar = c(0,0,0,0), new=TRUE)
plot(m_bin_ages[!is.na(m_bin_SD_diffs)], m_bin_SD_diffs[!is.na(m_bin_SD_diffs)], 
     ylim = c(0, max(m_bin_SD_diffs, na.rm = T)),
     type = "l", xlab = "", ylab = "SD(Difference)", lwd = 2, xaxt = "n")
lines(helper_plots$m_bin_ages[!is.na(helper_plots$m_bin_SD_diffs)], 
      helper_plots$m_bin_SD_diffs[!is.na(helper_plots$m_bin_SD_diffs)], 
      col = adjustcolor(1,0.25), lwd = 2, xpd = F)
abline(v=1:4*20, col = adjustcolor(1, 0.3), xpd = T)

helper_y_prop <- helper_y_prop + 0.16
par(fig=c(helper_x_prop, helper_y_prop), mar = c(0,0,0,0), new=TRUE)
plot(m_bin_ages, n_bin_prop_total / sum(n_bin_prop_total) * length(n_bin_prop_total),
     type = "l", xlab = "", ylab = "Relative Proportion\nSample (vs. Uniform)", lwd = 2, xaxt = "n")
lines(helper_plots$m_bin_ages, 
      helper_plots$n_bin_prop_total / sum(helper_plots$n_bin_prop_total) * length(helper_plots$n_bin_prop_total), 
      col = adjustcolor(1,0.25), lwd = 2, xpd = F)
abline(v=1:4*20, col = adjustcolor(1, 0.3), xpd = T)

helper_y_prop <- helper_y_prop + 0.16
par(fig=c(helper_x_prop, helper_y_prop), mar = c(0,0,0,0), new=TRUE)
plot(m_bin_ages, cumsum(n_bin_prop_total) / sum(n_bin_prop_total), ylim = c(0,1),
     type = "l", xlab = "", ylab = "Cumulative Proportion", lwd = 2, xaxt = "n")
lines(helper_plots$m_bin_ages, 
      cumsum(helper_plots$n_bin_prop_total) / sum(helper_plots$n_bin_prop_total), 
      col = adjustcolor(1,0.25), lwd = 2, xpd = F)
abline(v=1:4*20, col = adjustcolor(1, 0.3), xpd = T)
axis(3)

dev.off()

