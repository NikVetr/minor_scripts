#high level parameters
adj_for_inflation <- F
bw <- 1
ndays_to_merge <- 0

#specify functions
pretty_large_number <- function(num) {
  if (num == 0) return("0")
  
  if(num < 1) return(as.character(round(num, 2)))
  
  suffixes <- c("", "K", "M", "B", "T")
  index <- floor(log10(num) / 3)
  
  divisor <- 10^(index * 3)
  rounded_num <- round(num / divisor, 1)
  
  return(paste0(rounded_num, suffixes[index + 1]))
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

#read in data
d <- read.csv("~/Downloads/WAI_donations.csv")

#clean & process data file
d <- d[d$Amount != "",]
d$dollars <- as.numeric(gsub(",", "", trimws(d$Amount)))
d <- d[d$dollars > 1,] #some adjustments for things are negative, ignore for now
d$date <- as.Date(d$Date, format = "%m/%d/%Y")
d <- d[order(d$date),]
d$week <- as.numeric(format(d$date, "%V"))
d$year <- as.numeric(format(d$date, "%Y"))
d$day <- as.numeric(format(d$date, "%d"))
d$month <- as.numeric(format(d$date, "%m"))
d$year.month <- format(d$date, "%Y/%m")
d <- d[d$year.month < "2023/09",]

#incorporate inflation from https://www.bls.gov/regions/mid-atlantic/data/consumerpriceindexhistorical_us_table.htm
cpi <- read.csv("~/Downloads/CPI_USA.csv", header = T)
colnames(cpi)[colnames(cpi) == "X"] <- "Year"
cpi <- reshape(cpi, 
                     varying = list(names(cpi)[-1]), 
                     v.names = "CPI", 
                     idvar = "Year", 
                     direction = "long", 
                     timevar = "Month")
cpi$year.month <- format(as.Date(paste0(cpi$Year, "/", cpi$Month, "/01")), "%Y/%m")
d$dollars_feb2018_equiv <- d$dollars * cpi$CPI[match("2018/02", cpi$year.month)] / 
  cpi$CPI[match(d$year.month, cpi$year.month)] 
d$dollars_nominal <- d$dollars

if(adj_for_inflation){
  d$dollars <- d$dollars_feb2018_equiv
} else {
  d$dollars <- d$dollars_nominal
}

#use log dollars
d$logdol <- log10(d$dollars)

#merge all donations < 6d apart
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
          new_donation_chunk$dollars_feb2018_equiv <- 
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

# overall trend
plot(d$date, y = d$logdol, type = "l")

#find breaks / clusters?
kmeans_out <- lapply(2:20, function(k) kmeans(d$logdol, k, iter.max = 50, nstart = 50))
dbi <- sapply(kmeans_out, function(x) clusterSim::index.DB(d$logdol, x$cluster)$DB)
kmeans_picked <- which.min(dbi)

#hmm too many from kmeans, what does jenks say?
# classInt::classIntervals(d$logdol, 5, "jenks")
#jenks is too slow

#how about we just divide by quantiles? 
k <- 8
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

#split according to these groups
ds <- split(d, d$group)
plot(x = range(d$date), y = range(d$logdol), type = "l", col = adjustcolor(1,0))
for(i in 1:k){
  lines(ds[[i]]$date, ds[[i]]$logdol, col = i, lwd = 2)
}

#tells us little -- need # donors per unit time
all_months <- seq(from = as.Date(paste0(min(d$year.month), "/15")), 
                  to = max(d$date), by = "month")
dsc_n_indiv <- lapply(ds, function(x){
  out <- as.data.frame(table(x$year.month))
  out$month <- as.Date(paste0(as.character(out$Var1), "/15"), format = "%Y/%m/%d")
  out <- out[,c("month", "Freq")]
  out <- rbind(out, data.frame(month = all_months[!all_months %in% out$month], Freq = 0))
  return(out[order(out$month),])
})

#or the total sm in those bins?
dsc_n_tot <- lapply(ds, function(x){
  out <- sapply(split(x, x$year.month), function(y) sum(y$dollars))
  out <- data.frame(Var1 = names(out), Freq = out)
  out$month <- as.Date(paste0(as.character(out$Var1), "/15"), format = "%Y/%m/%d")
  out <- out[,c("month", "Freq")]
  out <- rbind(out, data.frame(month = all_months[!all_months %in% out$month], Freq = 0))
  return(out[order(out$month),])
})

#plot now
plot(x = range(d$date), y = log10(range(do.call(rbind, dsc)$Freq)), 
     type = "l", col = adjustcolor(1,0))
for(i in 1:k){
  lines(dsc[[i]]$month, log10(dsc[[i]]$Freq), col = i, lwd = 2)
}


#### export to pdf ####
plot_num_indiv <- T
if(plot_num_indiv){
  dsc <- dsc_n_indiv
} else {
  dsc <- dsc_n_tot
}
season_cols <- c("#DC3545", 
                 "#FFC107",
                 "#007BFF", 
                 "#28A745")

cols <- blend_with_color(viridisLite::plasma(k), 1, 0.7)$hex
grDevices::cairo_pdf(filename = paste0("~/Documents/WAI_donations.pdf"), 
                     width = 1500 / 72, height = 1500 / 72 * (k+1)/6 / 2, family="Arial Unicode MS", 
                     pointsize = 19)

par(mar = c(5,6,3,5))
layout(rbind(1:2,3:4,5:6,7:8,9:10), heights = c(1,1,1,1,0.1))
for(i in 1:k){
  
  plot(dsc[[i]]$month, (dsc[[i]]$Freq), col = cols[i], lwd = 2, type = "l", 
       frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       xlim = range(d$date), cex.lab = 1.5, ylim = c(0, max((dsc[[i]]$Freq))))
  mtext(paste0("$", pretty_large_number(10^breaks[i]), " - $", 
               pretty_large_number(10^breaks[i+1]), paste0(" USD ", ifelse(adj_for_inflation, "(Feb. 2018 equiv.)", "(nominal)"))), 
        line = 0.5, col = cols[i])
  
  #smoothed line
  lines(dsc[[i]]$month, dexp_smooth(as.integer(dsc[[i]]$month), (dsc[[i]]$Freq), r = 1/100),
        col = adjustcolor(cols[i], 0.25), lty = 3, lwd = 5)
  
  #vert axis
  vticks <- pretty(c(0, max(dsc[[i]]$Freq)))
  vticks <- vticks[vticks %% 1 == 0]
  if(plot_num_indiv){
    vticks_lab <- vticks
  } else {
    vticks_lab <- sapply(vticks, pretty_large_number)
  }
  axis(2, vticks, vticks_lab, cex.axis = 1.25, las = 2)
  yax_lab <- ifelse(plot_num_indiv, 
                    "# Donors / Month", 
                    "Total Donations / Month")
  mtext(text = yax_lab, side = 2, line = ifelse(plot_num_indiv, 
                                                3, 
                                                4))
  
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
  
  if(i %in% c(7,8)){
    mtext("Date", 1, line = 5.25, cex = 1.15)
  }
  
    
  
}
dev.off()

#### all in one graph ####
plot_num_indiv <- F
if(plot_num_indiv){
  dsc <- dsc_n_indiv
} else {
  dsc <- dsc_n_tot
}
season_cols <- c("#DC3545", 
                 "#FFC107",
                 "#007BFF", 
                 "#28A745")

cols <- blend_with_color(viridisLite::plasma(k), 1, 0.7)$hex
par(mar = c(5,6,3,5))
grDevices::cairo_pdf(filename = paste0("~/Documents/WAI_donations_all_together.pdf"), 
                     width = 1500 / 72 / 2, height = 1500 / 72 * (k+1)/6 / 2 / 4, family="Arial Unicode MS", 
                     pointsize = 16)

par(mar = c(5,5,3,7))
dsca <- do.call(rbind, dsc)
i = 1  
plot(dsc[[i]]$month, (dsc[[i]]$Freq), col = cols[i], lwd = 2, type = "l", 
     frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     xlim = range(dsca$month), cex.lab = 1.5, ylim = c(0, max(dsca$Freq)))
#smoothed line
lines(dsc[[i]]$month, dexp_smooth(as.integer(dsc[[i]]$month), (dsc[[i]]$Freq), r = 1/100),
      col = adjustcolor(cols[i], 0.25), lty = 3, lwd = 5)

for(i in 2:k){
  lines(dsc[[i]]$month, (dsc[[i]]$Freq), col = cols[i], lwd = 2)
  #smoothed line
  lines(dsc[[i]]$month, dexp_smooth(as.integer(dsc[[i]]$month), (dsc[[i]]$Freq), r = 1/100),
        col = adjustcolor(cols[i], 0.25), lty = 3, lwd = 5)
}

mtext(paste0("$", pretty_large_number(10^breaks[i]), " - $", 
             pretty_large_number(10^breaks[i+1]), paste0(" (USD, ", ifelse(adj_for_inflation, "Feb. 2018 equiv.)", "nominal)"))), 
      line = 0.5, col = cols[i])



#vert axis
vticks <- pretty(c(0, max(dsca$Freq)))
vticks <- vticks[vticks %% 1 == 0]
if(plot_num_indiv){
  vticks_lab <- vticks
} else {
  vticks_lab <- sapply(vticks, pretty_large_number)
}
axis(2, vticks, vticks_lab, cex.axis = 1.25, las = 2)
yax_lab <- ifelse(plot_num_indiv, 
                  "# Donors / Month", 
                  "Total Donations / Month")
mtext(text = yax_lab, side = 2, line = ifelse(plot_num_indiv, 
                                              3, 
                                              4))

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

if(i %in% c(7,8)){
  mtext("Date", 1, line = 5.25, cex = 1.15)
}

dev.off()

#### prop total donations ####

#remove donations >$0.5M
d <- d[d$dollars <= 5E5,]
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
                 rbind(0, apply(grouped_donations_prop, 2, cumsum)),
                 rbind(0, apply(grouped_donations_cumu_prop, 2, cumsum))
)
titles <- c("Monthly Donation Amounts",
            "Cumulative Monthly Donation Amounts",
            "Monthly Donation Proportions",
            "Cumulative Monthly Donation Proportions")
ylabs <- c(paste0("USD Donated ", ifelse(adj_for_inflation, "(Feb. 2018 equiv.)", "(nominal)")), 
           paste0("USD ", ifelse(adj_for_inflation, "(Feb. 2018 equiv.)", "(nominal)")), 
           "Proportion Donations", "Proportion Donations")

#actual plotting
grDevices::cairo_pdf(filename = paste0("~/Documents/WAI_donations_proptotal.pdf"), 
                     width = 1500 / 72, height = 1000 / 72, family="Arial Unicode MS", 
                     pointsize = 22)
par(mar = c(5,5.5,3,3))
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
  
  mtext(titles[i], line = 0.5, col = 1)
  
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
  
}


dev.off()

#####

#check if donors give more
indivs <- split(d[,c("dollars", "group", "year.month", "year", "date")], d$Name)
indivs <- split(d[d$year != 2023, c("dollars", "group", "year.month", "year", "date")], 
                d$Name[d$year != 2023]) #remove 2023
indivs <- indivs[sapply(indivs, function(x) length(unique(x$year.month))) >= 2]
spearman_corrs <- sapply(indivs, function(x){
  x <- x[order(x$year.month),]
  monthly_donations <- sapply(split(x$dollars, x$year.month), sum)
  cor(1:length(monthly_donations), monthly_donations, method = "spearman") 
})
hist(spearman_corrs, freq = F, main = "2023 Removed", xlab = "Spearman Correlation between Donation Month and Donation Amount")


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
hist(unlist(days_between_donations), breaks = 0:1000, xlim = c(0,100), 
     xaxt = "n", main = "", xlab = "# Days Between Donations")
axis(1, 0:10*10, 0:10*10)

plot(sort(table(d$Name)), xaxt = "n", ylab = "number of donations", xlab = "quantile of donor")
axis(1, length(unique(d$Name)) * 0:10/10, 0:10/10)

     