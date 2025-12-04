#high level parameters
org <- c("RP", "WAI")[1]
specify_bins <- T
k <- 7
bin_breaks <- c(0, 1E3, 1E4, 5E4, 1E5, 5E5, 1E6, 1E9)
remove_IAPS <- F
fig_dir <- paste0("~/Documents/Documents - nikolai-mba/", org, "-Donation-Figures/")
subfig_dir <- paste0(fig_dir, "subfigs/")
figure_paths <- character(0)
if(!dir.exists(subfig_dir)){dir.create(subfig_dir, recursive = T)}
remove_big_donations <- F
plot_num_indiv <- F
adj_for_inflation <- T
infl_ym <- "2015/01"
infl_ym_str <- paste0("(", month.abb[as.numeric(strsplit(infl_ym, "/")[[1]][2])], 
                      ". ", strsplit(infl_ym, "/")[[1]][1], " equiv.)")
bw <- 1
ndays_to_merge <- 0
fakesim <- F

#set overall figure names
fig_labels <- c(
  "Fig. 1: Binned Cumulative Donations, Resets each Year",
  "Fig. 2: Binned Monthly Donations",
  "Fig. 3: Binned Cumulative Donations, Years are Stacked",
  "Fig. 4: Binned Monthly Donations, Years are Stacked",
  "Fig. 5: Proportion of Funding from Each Bin"
)

#### functions ####
source("~/scripts/minor_scripts/postdoc/1Drepel_v2.R")

fig_lab <- function(xp = 0, yp = 0, cex = 2, label = "", xpd = NA, draw_grid = F, 
                    stepsize = 0.02, grid_sparsity = 3, ...){
  ds <- dev.size("in")
  # xy coordinates of device corners in user coordinates
  xb <- grconvertX(c(0, ds[1]), from="in", to="user")
  yb <- grconvertY(c(0, ds[2]), from="in", to="user")
  text(x = xb[1] + diff(xb) * xp / 100, 
       y = yb[1] + diff(yb) * yp / 100, 
       label = label, cex=cex, xpd = xpd)
  
  if(draw_grid){
    props <- seq(stepsize, 1 - stepsize, by = stepsize)
    ngrid <- length(props)
    x_pos_grid <- xb[1] + diff(xb) * props
    y_pos_grid <- yb[1] + diff(yb) * props
    
    # Draw grid lines
    abline(h = y_pos_grid, col = adjustcolor(1, 0.4), lty = 1, xpd = NA)
    abline(v = x_pos_grid, col = adjustcolor(1, 0.4), lty = 1, xpd = NA)
    
    #label pts in grid
    
    grid_lab_inds <- seq(1, ngrid, by = grid_sparsity)
    ngrid_pts <- length(grid_lab_inds)
    
    #hline labs
    points(x = rep(x_pos_grid[grid_lab_inds] + diff(x_pos_grid[1:2])/2, ngrid), 
           y = rep(y_pos_grid, each = ngrid_pts), 
           pch = 15, 
           xpd = NA, col = adjustcolor("white", 0.8), cex = 2)
    text(x = rep(x_pos_grid[grid_lab_inds] + diff(x_pos_grid[1:2])/2, ngrid), 
         y = rep(y_pos_grid, each = ngrid_pts), 
         labels = rep(props * 100, each = ngrid_pts), 
         xpd = NA, col = adjustcolor("darkred", 0.4), cex = 0.75)
    
    #yline labs
    points(y = rep(y_pos_grid[grid_lab_inds] + diff(y_pos_grid[1:2])/2, ngrid), 
           x = rep(x_pos_grid, each = ngrid_pts), 
           pch = 15, 
           xpd = NA, col = adjustcolor("white", 0.8), cex = 2)
    text(y = rep(y_pos_grid[grid_lab_inds] + diff(y_pos_grid[1:2])/2, ngrid), 
         x = rep(x_pos_grid, each = ngrid_pts), 
         labels = rep(props * 100, each = ngrid_pts), 
         xpd = NA, col = adjustcolor("darkblue", 0.4), cex = 0.75)
    
  }
  
}

format_percent <- function(x, sci_digits = 2) {
  # handle vector input and NAs
  out <- character(length(x))
  is_na <- is.na(x)
  if (any(is_na)) out[is_na] <- NA_character_
  
  small <- !is_na & abs(x) < 0.01
  big   <- !is_na & !small
  
  # small: scientific notation (E) with the % sign
  if (any(small)) {
    out[small] <- paste0(formatC(x[small], format = "e", digits = sci_digits), "%")
  }
  # big: one decimal place with the % sign
  if (any(big)) {
    out[big] <- paste0(format(round(x[big], 1), trim = TRUE, scientific = FALSE), "%")
  }
  out
}

pretty_large_number <- function(num) {
  if (num == 0) return("0")
  
  if(num < 1) return(as.character(round(num, 2)))
  
  suffixes <- c("", "K", "M", "B", "T")
  index <- floor(log10(num) / 3)
  
  divisor <- 10^(index * 3)
  rounded_num <- round(num / divisor, 1)
  
  return(paste0(rounded_num, suffixes[index + 1]))
}

prettier_large_number <- function(num, sig = 0) {
  # vectorized, robust handling
  suffixes <- c("", "K", "M", "B", "T", "P", "E")
  max_i <- length(suffixes) - 1L
  
  fmt_one <- function(n) {
    if (!is.finite(n)) return(as.character(n))
    if (n == 0) return("0")
    
    sign_chr <- if (n < 0) "-" else ""
    a <- abs(n)
    
    # tiny numbers: show up to 2 decimals without suffix
    if (a < 1) {
      s <- sub("\\.?0+$", "", format(round(a, 2), trim = TRUE, scientific = FALSE))
      return(paste0(sign_chr, s))
    }
    
    # choose suffix index by 1000s
    i <- max(0L, min(floor(log10(a) / 3), max_i))
    scaled <- a / (1000 ^ i)
    
    # round to desired significant digits
    r <- signif(scaled, sig)
    
    # rollover case: e.g., 999.5K -> 1000K -> bump to 1.00M
    if (r >= 1000 && i < max_i) {
      i <- i + 1L
      r <- signif(a / (1000 ^ i), sig)
    }
    
    # strip trailing zeros and any dangling dot
    s <- format(r, trim = TRUE, scientific = FALSE)
    s <- sub("\\.0+$", "", s)
    s <- sub("(\\.\\d*?)0+$", "\\1", s)
    
    paste0(sign_chr, s, suffixes[i + 1L])
  }
  
  vapply(num, fmt_one, character(1L))
}

dexp_smooth <- function(x, y, r, reweight_trunc_tail = F, reweight_n_pts = F, fix_endpoints = F, interpolate_at = NA){
  
  nx <- length(x)
  if(is.na(interpolate_at[1])){
    n <- nx
    w <- t(sapply(1:(n-1), function(i) c(dexp((x[i] - x[0:(i-1)]), rate = r), dexp(0, rate = r), dexp(-(x[i] - x[(i+1):n]), rate = r))))
    w <- rbind(w, dexp(x[n] - x, rate = r))
  } else {
    n <- length(interpolate_at)
    w <- t(sapply(1:(n-1), function(i) c(dexp((interpolate_at[i] - x[x<=interpolate_at[i]]), rate = r), 
                                         dexp((x[x>interpolate_at[i]] - interpolate_at[i]), rate = r))))
    w <- rbind(w, dexp(x[nx] - x, rate = r))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_trunc_tail){
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:(n-1), function(i) c(pexp((x[i] - x[1]), rate = r), pexp(-(x[i] - x[n]), rate = r))))
    } else {
      tw <- t(sapply(1:(n-1), function(i) c(pexp((interpolate_at[i] - x[1]), rate = r), pexp(-(interpolate_at[i] - x[nx]), rate = r))))
    }
    tw[1,] <- tw[2,]
    tw <- rbind(tw, tw[n-1,])
    tw <- 1/tw
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
    
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(reweight_n_pts){
    if(is.na(interpolate_at[1])){
      tw <- cbind(0:(n-1), (n-1):0)
      tw <- 1/tw
    } else {
      tw <- sapply(1:n, function(i) sum(interpolate_at[i] >= x))
      tw <- cbind(tw, nx - tw)
    }
    mintw <- apply(tw, 1, min)
    tw[,1] <- tw[,1] / mintw 
    tw[,2] <- tw[,2] / mintw 
    
    tw[1,] <- tw[2,]; tw[n,] <- tw[n-1,]
    
    if(is.na(interpolate_at[1])){
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], ri), rep(1, n-ri)) * c(rep(1, ri), rep(tw[ri,2], n-ri)))) #eh let's give it to the nth pt too
    } else {
      tw <- t(sapply(1:n, function(ri) c(rep(tw[ri,1], sum(interpolate_at[ri] >= x)), rep(1, nx-sum(interpolate_at[ri] >= x))) * 
                       c(rep(1, sum(interpolate_at[ri] >= x)), rep(tw[ri,2], nx-sum(interpolate_at[ri] >= x)))))
    }
    
    w <- t(sapply(1:n, function(ri) w[ri,] * tw[ri,]))
  }
  # dim(w)
  # plot(w[2300,])
  # plot(w[2400,])
  # plot(w[2500,])
  # plot(w[2600,])
  
  if(fix_endpoints){
    w[1,] <- c(1, rep(0,nx-1))
    w[n,] <- c(rep(0,nx-1), 1)
  }
  
  wy <- c(w %*% t(t(y)))
  wy <- wy / apply(w,1,sum)
  return(wy)
}


blend_with_color <- function(base_color, blend_color = "white", alpha) {
  # Split the base color into red, green, and blue components
  base_red <- col2rgb(base_color)["red", ] / 255
  base_green <- col2rgb(base_color)["green", ] / 255
  base_blue <- col2rgb(base_color)["blue", ] / 255
  
  # Split the blend color into red, green, and blue components
  blend_red <- col2rgb(blend_color)["red", ] / 255
  blend_green <- col2rgb(blend_color)["green", ] / 255
  blend_blue <- col2rgb(blend_color)["blue", ] / 255
  
  # Blend each component with the corresponding component of the blend color
  blended_red <- (alpha * base_red) + ((1 - alpha) * blend_red)
  blended_green <- (alpha * base_green) + ((1 - alpha) * blend_green)
  blended_blue <- (alpha * base_blue) + ((1 - alpha) * blend_blue)
  
  # Convert the blended components back to a hex color
  blended_color <- rgb(blended_red, blended_green, blended_blue)
  return(list(rgb = rbind(blended_red, blended_green, blended_blue), hex = blended_color))
}

#### load data ####
if(org == "WAI"){
  d <- read.csv("~/Downloads/WAI_donations.csv")
  d <- d[d$Date != "",]
  d <- d[d$Amount != "",]
  d$dollars <- as.numeric(gsub(",", "", trimws(d$Amount)))
  d$date <- as.Date(d$Date, format = "%m/%d/%Y")
}

if(org == "RP"){
  d <- read.csv("~/data/rp-donations_Oct-2025.csv")
  d$Date <- d$Gift.Date
  d$dollars <- as.numeric(gsub(",|\\$", "", trimws(d$Gift.Amount)))
  d$Name <- paste0(d$Hard.Credit.Name..Legal.Donor., ".",
                   d$Primary.Accredited.Donor)
  if(remove_IAPS){
    d <- d[d$Fund != "IAPS",]
  }
}

#simulate fake data for sharing purposes
if(fakesim){
  d$Date <- paste0(substr(d$Date, start = 1, nchar(d$Date)-2), 
                   sample(15:23, size = nrow(d), T, prob = 15:23-14))
  dy <- split(d, as.numeric(do.call(rbind, strsplit(d$Date, "/"))[,3]))
  # d <- do.call(rbind, lapply(names(dy), function(yi) {
  #   rinds <- sample(x = 1:nrow(dy[[yi]]), 
  #          size = nrow(dy[[yi]]) * (0.7 + (as.numeric(yi) - min(as.numeric(names(dy)))) / 1), 
  #          replace = T)
  #   dy[[yi]][rinds,]
  # }))
  d$dollars <- sample(d$dollars, length(d$dollars), T) * rlnorm(nrow(d), 1, 2)
}

#### process data ####

#clean & process data file
d <- d[d$dollars > 0,] #some adjustments for things are negative, ignore for now
d$date <- as.Date(d$Date, format = "%m/%d/%Y")
d <- d[order(d$date),]
d$week <- as.numeric(format(d$date, "%V"))
d$year <- as.numeric(format(d$date, "%Y"))
d$day <- as.numeric(format(d$date, "%d"))
d$month <- as.numeric(format(d$date, "%m"))
d$year.month <- format(d$date, "%Y/%m")

#skip the days from this month
d <- d[d$year.month < format(Sys.Date(), "%Y/%m"),]

#incorporate inflation from https://www.bls.gov/regions/mid-atlantic/data/consumerpriceindexhistorical_us_table.htm
cpi <- read.csv("~/data/CPI_USA.csv", header = T)
colnames(cpi)[colnames(cpi) == "X"] <- "Year"
cpi <- reshape(cpi, 
               varying = list(names(cpi)[-1]), 
               v.names = "CPI", 
               idvar = "Year", 
               direction = "long", 
               timevar = "Month")
cpi$year.month <- format(as.Date(paste0(cpi$Year, "/", cpi$Month, "/01")), "%Y/%m")

#project CPI to most recent month if missing
cpi <- cpi[order(cpi$year.month),]
cpi <- cpi[cpi$year.month <= max(d$year.month),]
fit <- with(cpi[!is.na(cpi$CPI),], { t <- 12*(Year - min(Year)) + Month; arima(CPI, order = c(1, 1, 0), 
                                                                               xreg = cbind(sin = sin(2*pi*t/12), 
                                                                                            cos = cos(2*pi*t/12))) })
fc  <- with(list(t = seq(with(cpi[!is.na(cpi$CPI),], max(12*(Year - min(Year)) + Month)) + 1, length.out = 12)),
            predict(fit, n.ahead = 12, newxreg = cbind(sin = sin(2*pi*t/12), cos = cos(2*pi*t/12))))
plot(cpi$CPI, type = "l", 
     xlim = c(0, nrow(cpi) + 12), 
     ylim = range(c(cpi$CPI, fc$pred + fc$se * 2, fc$pred - fc$se * 2), na.rm = T))
lines(sum(!is.na(cpi$CPI)) + 0:12, c(tail(cpi$CPI[!is.na(cpi$CPI)], 1), fc$pred), col = 2)
lines(sum(!is.na(cpi$CPI)) + 0:12, c(tail(cpi$CPI[!is.na(cpi$CPI)], 1), fc$pred + fc$se * 2), 
      col = 2, lty = 3)
lines(sum(!is.na(cpi$CPI)) + 0:12, c(tail(cpi$CPI[!is.na(cpi$CPI)], 1), fc$pred - fc$se * 2), 
      col = 2, lty = 3)
cpi$CPI[is.na(cpi$CPI)] <- fc$pred[1:sum(is.na(cpi$CPI))]

#account for inflation
d$dollars_infl_adj <- d$dollars * cpi$CPI[match(infl_ym, cpi$year.month)] / 
  cpi$CPI[match(d$year.month, cpi$year.month)]
d$dollars_nominal <- d$dollars

if(adj_for_inflation){
  d$dollars <- d$dollars_infl_adj
} else {
  d$dollars <- d$dollars_nominal
}

#remove negative values and values that are NA (eg from having no inflation data)
d <- d[!is.na(d$dollars),]
d <- d[d$dollars > 0,]

#use log dollars
d$logdol <- log10(d$dollars)

#merge all donations < 6d apart
merge_donations <- F
if(is.null(d$Name)){
  d$Name <- sample(1:(nrow(d)/3), nrow(d), replace = T)
}
if(merge_donations){
  if(is.null(d$Name)){
    d$Name <- sample(1:(nrow(d)/3), nrow(d), replace = T)
  }
  indivs <- split(d, d$Name)
  d <- do.call(rbind, lapply(indivs, function(x){
    same_donations <- diff(x$date) < ndays_to_merge
    if(any(same_donations)){
      rlesd <- rle(same_donations)
      sd_inds <- c(1,cumsum(rlesd$lengths))
      sd_inds <- cbind(sd_inds[-length(sd_inds)], sd_inds[-1])
      merged_conations <- do.call(rbind, lapply(seq_along(rlesd$values), function(tfi){
        donation_chunk <- x[sd_inds[tfi,1]:sd_inds[tfi,2],]
        if(rlesd$values[tfi]){
          new_donation_chunk <- donation_chunk[1,]
          weights <- donation_chunk$dollars / sum(donation_chunk$dollars)
          new_donation_chunk$date <- new_donation_chunk$date + 
            sum(c(0,diff(donation_chunk$date)) * weights)
          new_donation_chunk$dollars <- 
            new_donation_chunk$dollars_infl_adj <- 
            sum(donation_chunk$dollars)
          new_donation_chunk$dollars_nominal <- sum(donation_chunk$dollars_nominal)
          new_donation_chunk$logdol <- log10(new_donation_chunk$dollars)
          new_donation_chunk$week <- as.numeric(format(new_donation_chunk$date, "%V"))
          new_donation_chunk$year <- as.numeric(format(new_donation_chunk$date, "%Y"))
          new_donation_chunk$day <- as.numeric(format(new_donation_chunk$date, "%d"))
          new_donation_chunk$month <- as.numeric(format(new_donation_chunk$date, "%m"))
          new_donation_chunk$year.month <- format(new_donation_chunk$date, "%Y/%m")
          donation_chunk <- new_donation_chunk
        }
        return(donation_chunk)
      }))
      return(merged_conations)
    } else {
      return(x)
    }
  }))  
}

# overall trend
d <- d[order(d$date),]
plot(d$date, y = d$logdol, type = "l")
head(d[is.na(d$logdol),])

#### set k bins ####

#find breaks / clusters?
try_kmeans <- F
if(try_kmeans){
  kmeans_out <- lapply(2:20, function(k) kmeans(d$logdol, k, iter.max = 50, nstart = 50))
  dbi <- sapply(kmeans_out, function(x) clusterSim::index.DB(d$logdol, x$cluster)$DB)
  kmeans_picked <- which.min(dbi)
}

#hmm too many from kmeans, what does jenks say?
# classInt::classIntervals(d$logdol, 5, "jenks")
#jenks is too slow

#how about we just divide by quantiles? 
breaks <- quantile(d$logdol, 0:k/k)
hist(d$logdol, breaks = seq(min(d$logdol), max(d$logdol), length.out = 50))
abline(v = breaks, col = 2, lwd = 2) #too lumpy at the top
d$group <- sapply(d$logdol, function(x) min(sum(x >= breaks), k))

#modify to contain equal intervals?
disp <- mean(diff(breaks))
breaks <- seq(min(d$logdol), max(d$logdol), by = disp)
d$group <- sapply(d$logdol, function(x) min(sum(x >= breaks), length(breaks) - 1))
hist(d$logdol, breaks = seq(min(d$logdol), max(d$logdol), length.out = 50))
abline(v = breaks, col = 2, lwd = 2) #too lumpy at the top

#split according to total funding amount
# cumtotdol <- cumsum(sort(d$dollars))
# bin_size <- tail(cumtotdol, 1) / k
# cuts <- bin_size * (0:(k-1))
# sapply(cuts, function(ci) min(which(cumtotdol > ci)))
#hmm this is almost all just the largest transactions

disp <- mean(diff(breaks))
breaks <- seq(min(d$logdol), max(d$logdol), by = disp)
d$group <- sapply(d$logdol, function(x) min(sum(x >= breaks), length(breaks) - 1))
hist(d$logdol, breaks = seq(min(d$logdol), max(d$logdol), length.out = 50))
abline(v = breaks, col = 2, lwd = 2) #too lumpy at the top

#or just specify bins a priori
if(specify_bins){
  breaks <- log10(bin_breaks)
  k <- length(bin_breaks) - 1
  d$group <- sapply(d$logdol, function(x) min(sum(x >= breaks), length(breaks) - 1))
  hist(d$logdol, breaks = seq(min(d$logdol), max(d$logdol), length.out = 50))
  abline(v = breaks, col = 2, lwd = 2) #too lumpy at the top
}

#how much money goes in each bin?
if(adj_for_inflation){
  usd_in_each_bin <- sapply(split(d$dollars_infl_adj, d$group), sum)
} else {
  usd_in_each_bin <- sapply(split(d$dollars, d$group), sum)
}
percent_usd_in_each_bin <- c(usd_in_each_bin / sum(usd_in_each_bin) * 100, 100)
usd_in_each_bin <- c(usd_in_each_bin, sum(usd_in_each_bin))
bin_info <- paste0("∑ in Bin: $", prettier_large_number(usd_in_each_bin, 2), 
                   ",\n% All Donations: ",
                   format_percent(percent_usd_in_each_bin, sci_digits = 1))

#split according to these groups
ds <- split(d, d$group)
#add an extra group for all $
ds[[length(ds) + 1]] <- d
plot(x = range(d$date), y = range(d$logdol), type = "l", col = adjustcolor(1,0))
for(i in 1:(k+1)){
  lines(ds[[i]]$date, ds[[i]]$logdol, col = i, lwd = 2)
}

#tells us little -- need # donors per unit time
all_months <- seq(from = as.Date(paste0(min(d$year.month), "/15")), 
                  to = max(d$date), by = "month")
dsc_n_indiv <- lapply(ds, function(x){
  out <- as.data.frame(table(x$year.month))
  out$month <- as.Date(paste0(as.character(out$Var1), "/15"), format = "%Y/%m/%d")
  out <- out[,c("month", "Freq")]
  if(!all(all_months %in% out$month)){
    out <- rbind(out, data.frame(month = all_months[!all_months %in% out$month], Freq = 0))  
  }  
  return(out[order(out$month),])
})

#or the total sm in those bins?
dsc_n_tot <- lapply(ds, function(x){
  out <- sapply(split(x, x$year.month), function(y) sum(y$dollars))
  out <- data.frame(Var1 = names(out), Freq = out)
  out$month <- as.Date(paste0(as.character(out$Var1), "/15"), format = "%Y/%m/%d")
  out <- out[,c("month", "Freq")]
  if(!all(all_months %in% out$month)){
    out <- rbind(out, data.frame(month = all_months[!all_months %in% out$month], Freq = 0))  
  }
  return(out[order(out$month),])
})

if(plot_num_indiv){
  dsc <- dsc_n_indiv
} else {
  dsc <- dsc_n_tot
}

#plot now
plot(x = range(d$date), y = log10(range(do.call(rbind, dsc)$Freq)), 
     type = "l", col = adjustcolor(1,0))
for(i in 1:(k+1)){
  lines(dsc[[i]]$month, log10(dsc[[i]]$Freq), col = i, lwd = 2)
}


#### fig: monthly donations plot ####
pcbys <- c(T, F) #plot cumulative by year options
for(plot_cumulative_by_year in pcbys){
  
  if(plot_num_indiv){
    dsc <- dsc_n_indiv
  } else {
    dsc <- dsc_n_tot
  }
  
  season_cols <- c("#DC3545", 
                   "#FFC107",
                   "#007BFF", 
                   "#28A745")
  
  cols <- blend_with_color(c(viridisLite::plasma(k), 1), 1, 0.7)$hex
  fig_path <- paste0(subfig_dir, "donations_through_time", 
                     ifelse(plot_cumulative_by_year, "_cumulative-by-year", ""),".pdf")
  grDevices::cairo_pdf(filename = fig_path, 
                       width = 1500 / 72 * 0.85, height = 1500 / 72 * (ceiling(k/2)*2)/6 / 2 * 1.5 * 0.85, family="Arial Unicode MS", 
                       pointsize = 19)
  figure_paths <- c(figure_paths, fig_path)
  
  par(mar = c(6,6,1.5,4), 
      oma = c(1,ifelse(plot_cumulative_by_year, 3, 1),5,1), 
      mfrow = c(ceiling(k/2), 2))
  for(i in 1:(k+1)){
    
    #accumulate?
    if(plot_cumulative_by_year){
      #split by year
      dsc[[i]]$year <- as.numeric(format(dsc[[i]]$month, "%Y"))
      dsc[[i]]$actual_month <- as.numeric(format(dsc[[i]]$month, "%m"))
      dscxy <- split(dsc[[i]], dsc[[i]]$year)
      
      dscxy <- lapply(dscxy, function(temp){
        temp <- temp[order(temp$actual_month),]
        temp$Freq <- cumsum(temp$Freq)
        
        #add in initial 0 money month at the start
        temp <- rbind(temp[1,], temp)
        temp$month[1] <- temp$month[1] - 16
        temp$Freq[1] <- 0
        
        return(temp)
      })
      dsc[[i]] <- do.call(rbind, dscxy)
    }
    
    plot(dsc[[i]]$month+15, c(dsc[[i]]$Freq), col = cols[i], lwd = 2, type = "l", 
         frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         xlim = range(d$date), cex.lab = 1.5, ylim = c(0, max((dsc[[i]]$Freq))))
    
    if(i <= k){
      interval_label <- paste0("$", pretty_large_number(10^breaks[i]), " - $", 
             pretty_large_number(10^breaks[i+1]))
    } else {
      interval_label <- paste0("$", pretty_large_number(10^breaks[1]), " - $", 
                               pretty_large_number(10^breaks[length(breaks)]))
    }
    mtext(paste0(interval_label, 
                 paste0(" USD ", ifelse(adj_for_inflation, 
                                        infl_ym_str, 
                                        "(nominal)"))), 
          line = 0.5, col = cols[i])  
    
    #smoothed line
    lines(dsc[[i]]$month+15, dexp_smooth(as.integer(dsc[[i]]$month), (dsc[[i]]$Freq), r = 1/100),
          col = adjustcolor(cols[i], 0.25), lty = 3, lwd = 5)
    
    #vert axis
    vticks <- pretty(c(0, max(dsc[[i]]$Freq)))
    vticks <- vticks[vticks %% 1 == 0]
    if(plot_num_indiv){
      vticks_lab <- vticks
    } else {
      vticks_lab <- paste0("$", sapply(vticks, pretty_large_number))
    }
    axis(2, vticks, vticks_lab, cex.axis = 1.25, las = 2)
    yax_lab <- paste0(ifelse(plot_num_indiv, 
                             "# Donors / Month", 
                             "Total USD / Month"),
                      ifelse(plot_cumulative_by_year, "\n(cumulative by year)", ""))
    mtext(text = yax_lab, side = 2, line = ifelse(plot_num_indiv, 
                                                  3, 
                                                  5))
    
    #horiz axis
    line_height <- (par("usr")[4] - par("usr")[3]) / par("pin")[2] * par("cin")[2] * 0.75
    axis_dates <- seq(from = as.Date(paste0(format(min(d$date), "%Y"), "-01-01")), 
                      to = as.Date(paste0(as.integer(format(max(d$date), "%Y"))+1, "-01-01")), 
                      by = "year")
    axis(1, at = as.integer(axis_dates), labels = format(axis_dates, "%Y"),
         line = 2, cex.axis = 1.5, xpd = NA)
    segments(x0 = as.integer(axis_dates)[-c(1,length(axis_dates))], 
             x1 = as.integer(axis_dates)[-c(1,length(axis_dates))],
             y0 = par("usr")[3] - 2 * line_height, y1 = par("usr")[4], xpd = NA,
             lty = 1, lwd = 1, 
             col = adjustcolor(1, 0.75))
    season_dates <- seq(from = head(axis_dates, 1), 
                        to = tail(axis_dates, 1), by = "month")
    full_season_dates <- season_dates[grepl(pattern = "-01-|-04-|-07-|-10-", season_dates)]
    season_dates <- season_dates[grepl(pattern = "-04-|-07-|-10-", season_dates)]
    axis(1, at = as.integer(season_dates), labels = FALSE, tcl = -0.3, xpd = NA, line = 2)
    
    segments(x0 = as.integer(season_dates), x1 = as.integer(season_dates),
             y0 = par("usr")[3] - 2 * line_height, y1 = par("usr")[4], xpd = NA,
             lty = 3, lwd = 0.75, 
             col = adjustcolor(1, 0.5))
    
    for(j in 2:length(full_season_dates)){
      scol <- season_cols[j%%4+1]
      rect(xleft = full_season_dates[j-1], xright = full_season_dates[j],
           ybottom = -line_height*1.9, ytop = -line_height*1.5, 
           col = adjustcolor(scol, 0.5),
           xpd = NA, border = blend_with_color(scol, 1, 0.5)$hex)
    }
    
    #label horizontal axis
    if(i %in% (ceiling(k/2)*2-c(1,0))){
      mtext("Date", 1, line = 5.5, cex = 1.15, xpd = NA)  
    }
    
    #label panel
    panel_label <- paste0(letters[i], ")")
    pusr <- par("usr")
    text(x = pusr[1] - diff(pusr[1:2])/50, 
         y = pusr[4] + diff(pusr[3:4])/50, 
         pos = 3, label = panel_label, xpd = NA,
         cex = 2, font = 2)
    
    #give bin info
    text(x = pusr[1], 
         y = pusr[4] - diff(pusr[3:4])/10, 
         pos = 4, label = bin_info[i], xpd = NA,
         cex = 0.65)
    
  }
  
  #label figure
  fig_lab(xp = 50, yp = 97, label = fig_labels[length(figure_paths)])
  
  dev.off()
  
}

#### fig: stacking years ####
pcs <- c(T, F) #plot cumulative options
for(plot_cumulative in pcs){
  
  if(plot_num_indiv){
    dsc <- dsc_n_indiv
  } else {
    dsc <- dsc_n_tot
  }
  season_cols <- c("#DC3545", 
                   "#FFC107",
                   "#007BFF", 
                   "#28A745")
  nyears <- length(unique(d$year))
  year_cols <- setNames(viridisLite::viridis(nyears, begin = 0, end = 0.9), sort(unique(d$year)))
  
  cols <- blend_with_color(viridisLite::plasma(k), 1, 0.7)$hex
  fig_path <- paste0(subfig_dir, "donations-by-year_aggregated-by-",
                     ifelse(plot_num_indiv, "num-indiv", "total-amount"),
                     ifelse(plot_cumulative, "_cumulative", ""),
                     ".pdf")
  grDevices::cairo_pdf(filename = fig_path, 
                       width = 1500 / 72 / 2 * 2 * 0.85 , 
                       height = 1500 / 72 * (ceiling(k/2)*2)/6 / 2 * 1.5 * 0.85, family="Arial Unicode MS", 
                       pointsize = 16)
  figure_paths <- c(figure_paths, fig_path)
  
  par(mar = c(5,6,1.5,8), mfrow = c(ceiling(k/2), 2), oma = c(2,1,5,0))
  dsca <- do.call(rbind, dsc)
  
  for(i in 1:(k+1)){
    
    #split by year
    dsc[[i]]$year <- as.numeric(format(dsc[[i]]$month, "%Y"))
    dsc[[i]]$actual_month <- as.numeric(format(dsc[[i]]$month, "%m"))
    dscxy <- split(dsc[[i]], dsc[[i]]$year)
    
    #accumulate?
    if(plot_cumulative){
      dscxy <- lapply(dscxy, function(temp){
        temp <- temp[order(temp$actual_month),]
        temp$Freq <- cumsum(temp$Freq)
        return(temp)
      })
      dsc[[i]] <- do.call(rbind, dscxy)
    }
    
    plot(NULL, col = cols[i], lwd = 2, type = "l", 
         frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         xlim = c(1,12) + c(-0.5, 0.5), cex.lab = 1.5, ylim = c(0, max(dsc[[i]]$Freq)))
    
    for(j in names(dscxy)){
      year_data <- data.frame(month = dscxy[[j]]$actual_month,
                              amount = dscxy[[j]]$Freq)
      lines(c(0.5, year_data$month + 0.5), c(0, year_data$amount), col = year_cols[j], 
            lwd = 2)
    }
    
    #label the years
    final_month_info <- do.call(rbind, lapply(dscxy, function(temp){
      temp[which.max(temp$actual_month), c("actual_month", "year", "Freq")]}))
    minsep <- strheight("A") * 1.25
    xdisp <- strwidth("w")
    final_month_info$Freq
    laby <- minsep_fit(final_month_info$Freq, w = minsep, par("usr")[3:4])
    elbow_l <- 0.1 
    for(xi in 1:length(final_month_info$Freq)){
      bezier_curve(y0 = final_month_info$Freq[xi], 
                   x0 = final_month_info$actual_month[xi] + 0.5, 
                   y1 = laby[xi], 
                   x1 = par("usr")[2] + elbow_l*4, xpd = NA, 
                   ends = "flat", lty = 3, lwd = 2, 
                   col = year_cols[as.character(final_month_info$year[xi])])
    }
    text(x = par("usr")[2] + elbow_l*4 - xdisp / 2, y = laby, 
         labels = final_month_info$year, 
         xpd = NA, pos = 4, col = year_cols[as.character(final_month_info$year)],
         font = 2)
    
    #smoothed line
    # lines(dsc[[i]]$month, dexp_smooth(as.integer(dsc[[i]]$month), (dsc[[i]]$Freq), r = 1/100),
    #       col = adjustcolor(cols[i], 0.25), lty = 3, lwd = 5)
    
    #panel title
    if(i <= k){
      interval_label <- paste0("$", pretty_large_number(10^breaks[i]), " - $", 
                               pretty_large_number(10^breaks[i+1]))
    } else {
      interval_label <- paste0("$", pretty_large_number(10^breaks[1]), " - $", 
                               pretty_large_number(10^breaks[length(breaks)]))
    }
    mtext(paste0(interval_label, 
                 paste0(" USD ", ifelse(adj_for_inflation, 
                                        infl_ym_str, 
                                        "(nominal)"))), 
          line = 0.5, col = cols[i])  
    
    #vert axis
    vticks <- pretty(c(0, max(dsc[[i]]$Freq)))
    vticks <- vticks[vticks %% 1 == 0]
    if(plot_num_indiv){
      vticks_lab <- vticks
    } else {
      vticks_lab <- paste0("$", sapply(vticks, pretty_large_number))
    }
    axis(2, vticks, vticks_lab, cex.axis = 1.25, las = 2)
    yax_lab <- ifelse(plot_num_indiv, 
                      "# Donors / Month", 
                      paste0("Total USD / Month", 
                             ifelse(plot_cumulative, " (cumul.)", "")))
    mtext(text = yax_lab, side = 2, line = ifelse(plot_num_indiv, 3, 4.5))
    
    #horiz axis
    line_height <- diff(par("usr")[3:4])/100
    line_bot <- par("usr")[3] - 5*line_height
    
    axis_dates <- substr(format(ISOdate(2004,1:12,1),"%B"), 1, 3)
    axis_locs <- 1:12
    axis_ticks <- 0:12 + 0.5
    text(x = axis_locs, y = line_bot, labels = axis_dates, 
         xpd = NA, pos = 1, cex = 1.15)
    
    segments(0.5, line_bot, 12.5, line_bot, xpd = NA)
    segments(x0 = axis_ticks, 
             x1 = axis_ticks,
             y0 = line_bot - 2 * line_height, y1 = line_bot, xpd = NA,
             lty = 1, lwd = 1, 
             col = adjustcolor(1, 0.75))
    segments(x0 = axis_ticks, 
             x1 = axis_ticks,
             y0 = line_bot - 2 * line_height, y1 = par("usr")[4], xpd = NA,
             lty = 3, lwd = 1, 
             col = adjustcolor(1, 0.25))
    
    #plot seasonal boxes
    full_season_dates <- c(0,3,6,9,12)
    for(j in 2:length(full_season_dates)){
      scol <- season_cols[j-1]
      rect(xleft = full_season_dates[j-1] + 0.5, 
           xright = full_season_dates[j] + 0.5,
           ybottom = line_bot + line_height * 2, 
           ytop = line_bot + line_height * 5, 
           col = adjustcolor(scol, 0.5),
           xpd = NA, border = blend_with_color(scol, 1, 0.5)$hex)
    }
    
    #label horizontal axis
    if(i %in% (ceiling(k/2)*2-c(1,0))){
      mtext("Date", 1, line = 4.5, cex = 1.15, xpd = NA)  
    }
    
    #label panel
    panel_label <- paste0(letters[i], ")")
    pusr <- par("usr")
    text(x = pusr[1] - diff(pusr[1:2])/50, 
         y = pusr[4] + diff(pusr[3:4])/50, 
         pos = 3, label = panel_label, xpd = NA,
         cex = 2, font = 2)
    
    #give bin info
    text(x = pusr[1], 
         y = pusr[4] - diff(pusr[3:4])/10, 
         pos = 4, label = bin_info[i], xpd = NA,
         cex = 0.75)
    
  }
  
  #label figure
  fig_lab(xp = 50, yp = 97.5, label = fig_labels[length(figure_paths)])
  
  dev.off()
  
}

#### fig: prop total donations ####

#remove donations >$0.5M
if(remove_big_donations){
  d <- d[d$dollars <= 5E5,]  
}

# d <- d[!grepl("- Restricted", d$Account),]

all_year.months <- format(all_months, "%Y/%m")
grouped_donations <- do.call(cbind, lapply(all_year.months, function(ymi){
  subd <- d[d$year.month == ymi, c("dollars", "group")]
  total_donations <- sapply(split(subd$dollars, subd$group), sum)
  missing_groups <- setdiff(1:8, as.numeric(names(total_donations)))
  total_donations <- c(total_donations, 
                       setNames(rep(x = 0,length(missing_groups)), missing_groups))
  total_donations[order(as.numeric(names(total_donations)))]
}))
grouped_donations <- (t(apply(grouped_donations, 1, unlist)))
colnames(grouped_donations) <- all_year.months

#aggregate
grouped_donations_cumu <- t(apply(grouped_donations, 1, cumsum))
grouped_donations_prop <- grouped_donations %*% 
  diag(1 / (apply(grouped_donations, 2, sum) + 1E-9))
grouped_donations_cumu_prop <- grouped_donations_cumu %*% 
  diag(1 / apply(grouped_donations_cumu, 2, sum))

#smooth?
grouped_donations <- t(apply(as.matrix(grouped_donations), 1, function(x) dexp_smooth(1:length(x), x, bw)))
grouped_donations_cumu <- t(apply(as.matrix(grouped_donations_cumu), 1, function(x) dexp_smooth(1:length(x), x, bw)))
grouped_donations_prop <- t(apply(as.matrix(grouped_donations_prop), 1, function(x) dexp_smooth(1:length(x), x, bw)))
grouped_donations_cumu_prop <- t(apply(as.matrix(grouped_donations_cumu_prop), 1, function(x) dexp_smooth(1:length(x), x, bw)))

#plotting params
cols <- viridisLite::plasma(k)
plotmats <- list(rbind(0, apply(grouped_donations, 2, cumsum)),
                 rbind(0, apply(grouped_donations_cumu, 2, cumsum)),
                 rbind(0, apply(grouped_donations_prop, 2, function(mi) cumsum(mi / sum(mi)))),
                 rbind(0, apply(grouped_donations_cumu_prop, 2, cumsum))
)
titles <- c("Monthly Donation Amounts",
            "Cumulative Monthly Donation Amounts",
            "Monthly Donation Proportions",
            "Cumulative Monthly Donation Proportions")
ylabs <- c(paste0("USD Donated ", ifelse(adj_for_inflation, infl_ym_str, "(nominal)")), 
           paste0("USD ", ifelse(adj_for_inflation, infl_ym_str, "(nominal)")), 
           "Proportion Donations", "Proportion Donations")

#actual plotting
fig_path <- paste0(subfig_dir, "donations_through_time_proptotal.pdf")
grDevices::cairo_pdf(filename = fig_path, 
                     width = 1500 / 72 * 0.85, height = 1000 / 72 * 0.85, family="Arial Unicode MS", 
                     pointsize = 22)
figure_paths <- c(figure_paths, fig_path)
par(mar = c(5,5.5,3,3), oma = c(2,1,3.5,1))
layout(rbind(c(1,2,5), c(3,4,6)), widths = c(1,1,0.2))

for(i in 1:4){
  
  plotmat <- plotmats[[i]]
  plot(0,0, col = adjustcolor(1,0), 
       frame.plot = F, xaxt = "n", yaxt = "n", xlab = "Date", ylab = "",
       xlim = range(d$date), cex.lab = 1.25, ylim = c(0, max(plotmat)))
  mtext(ylabs[i], 2, line = 4)
  for(j in 1:k){
    xcoords <- as.Date(paste0(c(all_year.months, rev(all_year.months)), "/15"))
    ycoords <- c(plotmat[j,], rev(plotmat[j+1,]))
    polygon(x = xcoords, y = ycoords, border = blend_with_color(cols[j], "black", 0.8)$hex, 
            col = adjustcolor(cols[j], 0.5))
  }
  
  mtext(titles[i], line = 0.5, col = 1, cex = 0.9)
  
  #vert axis
  vticks <- pretty(c(0, max(plotmat)))
  axis(2, labels = paste0(ifelse(grepl("Amount", titles[i]), "$", ""), 
                          sapply(vticks, pretty_large_number)), vticks, las = 2)
  
  #horiz axis
  axis_dates <- seq(from = as.Date(paste0(format(min(d$date), "%Y"), "-01-01")), 
                    to = as.Date(paste0(as.integer(format(max(d$date), "%Y"))+1, "-01-01")), 
                    by = "year")
  axis(1, at = as.integer(axis_dates), labels = format(axis_dates, "%Y"),
       line = 0, cex.axis = 1.25, xpd = NA)
  segments(x0 = as.integer(axis_dates)[-c(1,length(axis_dates))], 
           x1 = as.integer(axis_dates)[-c(1,length(axis_dates))],
           y0 = par("usr")[3], y1 = par("usr")[4], xpd = NA,
           lty = 1, lwd = 1, 
           col = adjustcolor(1, 0.75))
  season_dates <- seq(from = head(axis_dates, 1), 
                      to = tail(axis_dates, 1), by = "month")
  season_dates <- season_dates[grepl(pattern = "-04-|-07-|-10-", season_dates)]
  axis(1, at = as.integer(season_dates), labels = FALSE, tcl = -0.3, xpd = NA, line = 0)
  
  segments(x0 = as.integer(season_dates), x1 = as.integer(season_dates),
           y0 = par("usr")[3], y1 = par("usr")[4], xpd = NA,
           lty = 3, lwd = 0.75, 
           col = adjustcolor(1, 0.5))
  
  #legend
  if(i==2){
    legend(x = par("usr")[2], y = par("usr")[4], xpd = NA,
           legend = sapply(1:k, function(bi){
             paste0("$", pretty_large_number(10^breaks[bi]), " - $", 
                    pretty_large_number(10^breaks[bi+1]))
           }),
           pt.bg = adjustcolor(cols, 0.5), 
           col = blend_with_color(cols, "black", 0.8)$hex,bty = "n",
           pch = 22, cex = 1.25, pt.cex = 2, title = "Donation Bins")
  }
  
  #label panel
  panel_label <- paste0(letters[i], ")")
  pusr <- par("usr")
  text(x = pusr[1] - diff(pusr[1:2])/50, 
       y = pusr[4] + diff(pusr[3:4])/50, 
       pos = 3, label = panel_label, xpd = NA,
       cex = 2, font = 2)
  
}

#label figure
fig_lab(xp = 50, yp = 95, label = fig_labels[length(figure_paths)])

dev.off()

#### fig: donor-specific analysis ####

#check if donors give more
remove_2023 <- F
if(remove_2023){
  indivs <- split(d[d$year != 2023, c("dollars", "group", "year.month", "year", "date")],
                  d$Name[d$year != 2023]) #remove 2023
} else {
  indivs <- split(d[,c("dollars", "group", "year.month", "year", "date")], d$Name)
}
indivs <- indivs[sapply(indivs, function(x) length(unique(x$year.month))) >= 2]
spearman_corrs <- sapply(indivs, function(x){
  x <- x[order(x$year.month),]
  monthly_donations <- sapply(split(x$dollars, x$year.month), sum)
  cor(1:length(monthly_donations), monthly_donations, method = "spearman") 
})
hist(spearman_corrs, freq = F, main = ifelse(remove_2023, "2023 Removed", ""), 
     xlab = "Spearman Correlation between\nDonation Month and Donation Amount")


going_up <- sapply(indivs, function(x){
  x <- x[order(x$year.month),]
  indices_of_group_12 <- which(x$group %in% c(1,2))
  if(length(indices_of_group_12) == 0){return(NA)}
  indices_of_group_g12 <- which(!(x$group %in% c(1,2)))
  if(length(indices_of_group_g12) == 0){return(0)}
  if(mean(indices_of_group_12) < mean(indices_of_group_g12)){
    return(2)
  }
  if(mean(indices_of_group_12) == mean(indices_of_group_g12)){
    return(1)
  }
  return(-1)
})
table(going_up)
table(going_up) / sum(table(going_up))

#recurring donors?
top_donors <- indivs[order(sapply(indivs, nrow), decreasing = T)]
days_between_donations <- sapply(top_donors, function(x) diff(x$date))
days_between_donations <- days_between_donations[sapply(days_between_donations, length) > 1]
hist(unlist(days_between_donations), breaks = 0:10000, xlim = c(0,100), 
     xaxt = "n", main = "", xlab = "# Days Between Donations")
axis(1, 0:10*10, 0:10*10)

plot(sort(table(d$Name)), xaxt = "n", ylab = "number of donations", xlab = "quantile of donor")
axis(1, length(unique(d$Name)) * 0:10/10, 0:10/10)

#### combine all figures ####

#combine pdfs
all_figures_path <- paste0(fig_dir, org, 
                           ifelse(org == "RP", ifelse(remove_IAPS, "-without-IAPS", "-with-IAPS"), ""),
                           "_fundraising-summary.pdf")
pdftools::pdf_combine(input = figure_paths, output = all_figures_path)

