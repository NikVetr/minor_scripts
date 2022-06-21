#load libraries

#specify functions
log_alt <- function(x, b = 10){
  logvals <- rep(NA, length(x))
  logvals[x <= b & x >= 0] <- x[x <= b & x >= 0] / b
  logvals[x > b] <- log(x[x > b],b)
  return(logvals)
}

#load data & snag metadata
d <- read.csv("~/Downloads/FP Climate Phil Report Data - FP Data (with additions).csv")
d$Amount.ClimateWorks.original.data...Rescaled. <- d$Amount.ClimateWorks.original.data...M. * 10
d <- d[d$Region != "" & d$Sector != "",]
regions <- setNames(unique(d$Region), unique(d$Region))
regions <- regions[-match(c("Global", "All", "Other/ Unknown"), regions)]
sectors <- setNames(unique(d$Sector), unique(d$Sector))
sectors <- sectors[-match(c("All sectors"), sectors)]

#munge data into relative formats
per_region_total_spending <- d[d$Sector == "All sectors" & d$Region %in% regions,c("Region", "Amount.ClimateWorks.original.data...Rescaled.")]
per_region_total_spending <- setNames(per_region_total_spending$Amount.ClimateWorks.original.data...Rescaled., per_region_total_spending$Region)
per_sector_total_spending <- d[d$Region == "All" & d$Sector %in% sectors,c("Sector", "Amount.ClimateWorks.original.data...Rescaled.")]
per_sector_total_spending <- setNames(per_sector_total_spending$Amount.ClimateWorks.original.data...Rescaled., per_sector_total_spending$Sector)

#reorder things
regions <- regions[names(sort(per_region_total_spending))]
sectors <- sectors[names(sort(per_sector_total_spending))]

#set colors
cols <- list(regions = setNames(pals::cols25(length(regions)), regions),
             sectors = setNames(pals::cols25(length(sectors)), sectors))
cols$sectors["Other Climate Change Mitigation Strategies"] <- "brown"

#find cumulative sum for stacked graph?
stacked = T
log_axis = F
if(log_axis){
  base = 10
} else {
  base = 5000
}
d$Amount.ClimateWorks.original.data...Rescaled.CumSum.RxS <- d$Amount.ClimateWorks.original.data...Rescaled.
for(s in sectors){
  d[d$Sector == s,]$Amount.ClimateWorks.original.data...Rescaled.CumSum.RxS[match(regions, d[d$Sector == s,]$Region)] <-
    cumsum(d[d$Sector == s,]$Amount.ClimateWorks.original.data...Rescaled.CumSum.RxS[match(regions, d[d$Sector == s,]$Region)])
}
d$Amount.ClimateWorks.original.data...Rescaled.CumSum.SxR <- d$Amount.ClimateWorks.original.data...Rescaled.
for(r in regions){
  d[d$Region == r,]$Amount.ClimateWorks.original.data...Rescaled.CumSum.SxR[match(sectors, d[d$Region == r,]$Sector)] <-
    cumsum(d[d$Region == r,]$Amount.ClimateWorks.original.data...Rescaled.CumSum.SxR[match(sectors, d[d$Region == r,]$Sector)])
}

  

#### region_by_sector figure ####
grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/RPOP-climate_viz-region_by_sector.pdf"), 
                     width = 1200 / 72, height = 700 / 72, family="Arial Unicode MS")

#specify graphical parameters
par(mar = c(4,5,1,18))
hax <- regions
vax <- sectors
if(stacked){
  d$Amount.ClimateWorks.original.data...Rescaled. <- d$Amount.ClimateWorks.original.data...Rescaled.CumSum.SxR  
}

#draw plot and axes
ylims <- log_alt(b = base, c(0, max(d[d$Region %in% hax & d$Sector %in% sectors,"Amount.ClimateWorks.original.data...Rescaled."])))
xlims <- c(1,length(hax))
plot(0, 0, col = "white", xaxt = "n", yaxt = "n", xlim = xlims, bty="n",
     ylim = ylims, xlab = "", ylab = "")
segments(x0 = 1:length(hax), x1 = 1:length(hax), y0 = ylims[1], y1 = ylims[2], col =  adjustcolor(1,0.2), lwd = 2, lty = 2, xpd = F)

#vertical axis
if(log_axis){
  segments(y0 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
           y1 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
           x0 = xlims[1], x1 = xlims[2],
           col =  adjustcolor(1,0.2), lwd = 2, lty = 2)
  segments(x0 = xlims[1], x1 = xlims[1], y0 = ylims[1], y1 = ylims[2], lwd = 2)
  segments(x0 = xlims[1], x1 = xlims[1] - diff(xlims) / 45, lwd = 2,
           y0 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
           y1 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1))
  oom <- seq(0, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1)
  
  minor_ticks <- (c(replicate(length(oom), (1:9)) %*% diag(10^(oom))))
  minor_ticks <- setdiff(minor_ticks, 10^oom[-1])
  minor_ticks <- minor_ticks[minor_ticks < 10^ylims[2]]
  text(x = xlims[1] - diff(xlims) / 50, labels = 10^seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
       y = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1), pos = 2, xpd = NA,
       cex = 1.25, font = 2)
  segments(y0 = log_alt(b = base, minor_ticks), y1 = log_alt(b = base, minor_ticks), 
           x0 = xlims[1], x1 = xlims[1] - diff(xlims) / 200, lwd = 1)
  text(y = log_alt(b = base, minor_ticks), x = xlims[1], labels = minor_ticks, pos = 2, cex = 0.5, xpd = NA)
  segments(y0 = log_alt(b = base, minor_ticks), y1 = log_alt(b = base, minor_ticks), 
           x0 = xlims[1], x1 = xlims[2], lwd = 1, lty = 3, col = adjustcolor(1,0.2))
} else {
  
  major_ticks <- log_alt(b = base, seq(0, floor(base*ylims[2]), by = 500))
  segments(x0 = xlims[1], x1 = xlims[1] - diff(xlims) / 50, lwd = 2,
           y0 = major_ticks,
           y1 = major_ticks)
  text(y = major_ticks, x = xlims[1]- diff(xlims) / 50, labels = round(major_ticks * base), pos = 2, cex = 1, xpd = NA)
  segments(y0 = major_ticks,
           y1 = major_ticks,
           x0 = xlims[1], x1 = xlims[2],
           col =  adjustcolor(1,0.2), lwd = 2, lty = 2)
  
}


#draw lines
for(x in vax){
  sub <- d[d$Sector == x,]
  if(stacked){
    sub_prev <- d[d$Sector == vax[match(x, vax)-1],]
    if(x != vax[1]){
      polygon(x = c(1:length(hax), length(hax):1), 
              y = c(log_alt(b = base, sub$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub$Region)]),
                    rev(log_alt(b = base, sub_prev$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub_prev$Region)]))),
            col = adjustcolor(cols$sectors[x], 0.75), lwd = 2.5)
    } else{
      polygon(x = c(1:length(hax), length(hax):1), 
              y = c(log_alt(b = base, sub$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub$Region)]), rep(0, length(hax))),
              col = adjustcolor(cols$sectors[x], 0.75), lwd = 2.5)
    }
  }
  else{
    lines(x = 1:length(hax), y = log_alt(b = base, sub$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub$Region)]),
          col = cols$sectors[x], lwd = 2.5)
  }
    
}

#hax labels
text(1:length(hax), y = ylims[1] - diff(ylims)/100, labels = gsub(hax, pattern = " ", replacement = "\n"), pos = 1, xpd = NA, cex = 1.25)

#vax labels
vax_locs <- d[d$Region == tail(hax,1) & d$Sector %in% sectors,]
vax_locs <- log_alt(b = base, vax_locs[match(sectors, vax_locs$Sector),"Amount.ClimateWorks.original.data...Rescaled."])
if(stacked){
  vax_locs <- (c(0, vax_locs[-length(vax_locs)]) + vax_locs) / 2
}
vax_locs_scale <- 100 / max(vax_locs)
vax_names <- FField::FFieldPtRep(coords = cbind(xlims[2] + 0.25, vax_locs * vax_locs_scale + rnorm(length(vax))/10), 
                                 rep.fact = 2, adj.max = 2)
vax_names$y <- vax_names$y / vax_locs_scale
text(labels = sectors, x = vax_names$x, y = vax_names$y, xpd = NA, pos = 4, col = cols$sectors[sectors])
segments(x0 = length(hax), x1 = vax_names$x, y0 = vax_locs, y1 = vax_names$y, col = cols$sectors[sectors], lwd = 2, lty = 3)

#vertical label
text(x = xlims[1] - diff(xlims)/11, y = mean(ylims), srt = 90, pos = 3,
     "Climate Philanthropic Spending ($100,000 USD)", xpd = NA, cex = 2)

segments(x0 = xlims[1], x1 = xlims[1], y0 = ylims[1], y1 = ylims[2], lwd = 2)
if(stacked){
  segments(x0 = xlims[2], x1 = xlims[2], y0 = ylims[1], y1 = ylims[2], lwd = 2)
}

dev.off()




grDevices::cairo_pdf(filename = paste0("~/Documents/Documents - nikolai/RPOP-climate_viz-sector_by_region.pdf"), 
                     width = 1200 / 72, height = 700 / 72, family="Arial Unicode MS")


par(mar = c(4,5,1,7))
if(stacked){
  d$Amount.ClimateWorks.original.data...Rescaled. <- d$Amount.ClimateWorks.original.data...Rescaled.CumSum.RxS  
}
hax <- sectors
vax <- regions
ylims <- log_alt(b = base, c(0, max(d[d$Sector %in% hax & d$Region %in% regions,"Amount.ClimateWorks.original.data...Rescaled."])))
xlims <- c(1,length(hax))
plot(0, 0, col = "white", xaxt = "n", yaxt = "n", xlim = xlims, bty="n",
     ylim = ylims, xlab = "", ylab = "")
segments(x0 = 1:length(hax), x1 = 1:length(hax), y0 = ylims[1], y1 = ylims[2], col =  adjustcolor(1,0.2), lwd = 2, lty = 2, xpd = F)


#vertical axis
if(log_axis){
  segments(y0 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
           y1 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
           x0 = xlims[1], x1 = xlims[2],
           col =  adjustcolor(1,0.2), lwd = 2, lty = 2)
  segments(x0 = xlims[1], x1 = xlims[1], y0 = ylims[1], y1 = ylims[2], lwd = 2)
  segments(x0 = xlims[1], x1 = xlims[1] - diff(xlims) / 45, lwd = 2,
           y0 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
           y1 = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1))
  oom <- seq(0, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1)
  
  minor_ticks <- (c(replicate(length(oom), (1:9)) %*% diag(10^(oom))))
  minor_ticks <- setdiff(minor_ticks, 10^oom[-1])
  minor_ticks <- minor_ticks[minor_ticks < 10^ylims[2]]
  text(x = xlims[1] - diff(xlims) / 50, labels = 10^seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1),
       y = seq(1, max(log_alt(b = base, d$Amount.ClimateWorks.original.data...Rescaled.)), by = 1), pos = 2, xpd = NA,
       cex = 1.25, font = 2)
  segments(y0 = log_alt(b = base, minor_ticks), y1 = log_alt(b = base, minor_ticks), 
           x0 = xlims[1], x1 = xlims[1] - diff(xlims) / 200, lwd = 1)
  text(y = log_alt(b = base, minor_ticks), x = xlims[1], labels = minor_ticks, pos = 2, cex = 0.5, xpd = NA)
  segments(y0 = log_alt(b = base, minor_ticks), y1 = log_alt(b = base, minor_ticks), 
           x0 = xlims[1], x1 = xlims[2], lwd = 1, lty = 3, col = adjustcolor(1,0.2))
} else {
  
  major_ticks <- log_alt(b = base, seq(0, floor(base*ylims[2]), by = 200))
  segments(x0 = xlims[1], x1 = xlims[1] - diff(xlims) / 50, lwd = 2,
           y0 = major_ticks,
           y1 = major_ticks)
  text(y = major_ticks, x = xlims[1]- diff(xlims) / 50, labels = round(major_ticks * base), pos = 2, cex = 1, xpd = NA)
  segments(y0 = major_ticks,
           y1 = major_ticks,
           x0 = xlims[1], x1 = xlims[2],
           col =  adjustcolor(1,0.2), lwd = 2, lty = 2)
}



#draw lines
for(x in vax){
  sub <- d[d$Region == x,]
  if(stacked){
    sub_prev <- d[d$Region == vax[match(x, vax)-1],]
    if(x != vax[1]){
      polygon(x = c(1:length(hax), length(hax):1), 
              y = c(log_alt(b = base, sub$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub$Sector)]),
                    rev(log_alt(b = base, sub_prev$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub_prev$Sector)]))),
              col = adjustcolor(cols$regions[x], 0.75), lwd = 2.5)
    } else{
      polygon(x = c(1:length(hax), length(hax):1), 
              y = c(log_alt(b = base, sub$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub$Sector)]), rep(0, length(hax))),
              col = adjustcolor(cols$regions[x], 0.75), lwd = 2.5)
    }
  }
  else{
    lines(x = 1:length(hax), y = log_alt(b = base, sub$Amount.ClimateWorks.original.data...Rescaled.[match(hax, sub$Sector)]),
        col = cols$regions[x], lwd = 2.5)
  }
}
# segments(x0 = xlims[1], x1 = xlims[2], y0 = ylims[1], y1 = ylims[1], lwd = 2)
new_hax_labels <- gsub(gsub(hax, pattern = " ", replacement = "\n"), pattern = "-", replacement = "-\n")
text(1:length(hax), y = ylims[1] - diff(ylims)/100, labels = new_hax_labels, pos = 1, xpd = NA, cex = 0.9)

#vax labels
vax_locs <- d[d$Sector == tail(hax,1) & d$Region %in% regions,]
vax_locs <- log_alt(b = base, vax_locs[match(regions, vax_locs$Region),"Amount.ClimateWorks.original.data...Rescaled."])
if(stacked){
  vax_locs <- (c(0, vax_locs[-length(vax_locs)]) + vax_locs) / 2
}

vax_locs_scale <- 100 / max(vax_locs)
vax_names <- FField::FFieldPtRep(coords = cbind(xlims[2] + 0.25, vax_locs * vax_locs_scale + rnorm(length(vax))/10), 
                                 rep.fact = 2, adj.max = 2)
vax_names$y <- vax_names$y / vax_locs_scale

text(labels = regions, x = vax_names$x, y = vax_names$y, xpd = NA, pos = 4, col = cols$regions[regions])
segments(x0 = length(hax), x1 = vax_names$x, y0 = vax_locs, y1 = vax_names$y, col = cols$regions[regions], lwd = 2, lty = 3)

#vertical label

text(x = xlims[1] - diff(xlims)/12, y = mean(ylims), srt = 90, pos = 3,
     "Climate Philanthropic Spending ($100,000 USD)", xpd = NA, cex = 2)

segments(x0 = xlims[1], x1 = xlims[1], y0 = ylims[1], y1 = ylims[2], lwd = 2)
if(stacked){
  segments(x0 = xlims[2], x1 = xlims[2], y0 = ylims[1], y1 = ylims[2], lwd = 2)
}


dev.off()

pdftools::pdf_combine(input = c("~/Documents/Documents - nikolai/RPOP-climate_viz-sector_by_region.pdf",
                                "~/Documents/Documents - nikolai/RPOP-climate_viz-region_by_sector.pdf"),
                      output = paste0("~/Documents/Documents - nikolai/RPOP-climate_viz",
                                      ifelse(stacked, "-stacked", ""),
                                      ifelse(log_axis, "", "-raw_axis"),".pdf"))
    