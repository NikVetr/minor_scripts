#problem -- svg plots things with a buffering outer box

#set box boundaries
xlim <- c(-3,10)
ylim <- c(-5,17)
npix <- c(w = 1367, h = 3761)

#plot a temporary svg file
temp_filename <- paste0(tempfile(), ".svg")
svg(filename = temp_filename, width = npix["w"] / 72, height = npix["h"] / 72, )
par(mar = c(0,0,0,0))
plot.new()
plot.window(xlim = xlim,
            ylim = ylim,
            xaxs="i", yaxs="i")
par("usr")
rect(xleft = xlim[1], xright = xlim[2], ybottom = ylim[1], ytop = ylim[2])
dev.off()

#read back in and retrieve relative coordinates
svgdf <- svgparser::read_svg(temp_filename, obj_type = 'data.frame')
bbox <- data.frame(x = svgdf$x[1:4], y = svgdf$y[1:4])
xrbb <- range(bbox$x)
yrbb <- range(bbox$y)
coords <- split(data.frame(x = svgdf$x, y = svgdf$y), f = svgdf$elem_idx)[-1]
coords <- do.call(rbind, coords)

#rescale to outer boundary
coords$x <- (coords$x - xrbb[1]) / diff(xrbb)
coords$y <- (coords$y - yrbb[1]) / diff(yrbb)
plot(coords, xlim = c(0,1), ylim = c(0,1), type = "l")
points(coords)
coords
11/12
1/12

#conclusion -- SVGs add 1/12 of the dimension to the outside of the plotting region