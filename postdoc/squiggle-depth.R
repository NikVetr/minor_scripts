library(sf)
library(igraph)

#functions
source("~/repos/polylines/R/functions.R")

#### more elaborate examples ####
npts <- 200

#initialize plot
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 4))

#overlapping squiggle
x <- rescale01(sin(seq(0, 4 * pi, length.out = npts))) * 2 * xyrat() +
  rescale01(1:npts) * 6
y <- rescale01(cos(seq(0, 4 * pi, length.out = npts))) * 1.5 
lwd <- (rescale01(max(y) - y) * 0.1 + 0.05) * plogis(1:npts/2 - 5) * plogis(npts:1/2 - 5) * 4

polylines(x, y, lwd = lwd, complex = F, xpd= NA, 
          col = adjustcolor(1, 0.2), border = 1)




