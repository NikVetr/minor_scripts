
n <- 1E2
x <- matrix(rnorm(2*n), n, 2)

plot(NULL, ann=FALSE,
     xlim=quantile(x[,1], c(0.1,1)), ylim=quantile(x[,2], c(0.1,1)))
polygon(c(1,2,3), c(1,3,2), col = adjustcolor(1,0.5))

xpd = NA
usr <- par("usr")
upct <- par("plt")
gr_usr <- usr + c(diff(usr[1:2]) / diff(upct[1:2]) * (c(0,1) - upct[1:2]),
                  diff(usr[3:4]) / diff(upct[3:4]) * (c(0,1) - upct[3:4]))



tmp <- tempfile()
# png(tmp, width = par("pin")[1], height = par("pin")[2], units = "in", res = 100, bg = "transparent") #if restricted to plotting window
png(tmp, width = par("din")[1], height = par("din")[2], units = "in", res = 100, bg = "transparent")

par(mar = c(0,0,0,0), xpd = NA)
plot.new()
# plot.window(xlim=usr[1:2], ylim=usr[3:4], xaxs = "i", yaxs = "i")
plot.window(xlim=gr_usr[1:2], ylim=gr_usr[3:4], xaxs = "i", yaxs = "i")


points(x, pch = 19, col = adjustcolor(1,0.5), cex = 3)
rect(usr[1], usr[3], usr[2], usr[4])
dev.off()

# rasterImage(png::readPNG(tmp), usr[1], usr[3], usr[2], usr[4], interpolate=F, xpd = xpd)
rasterImage(png::readPNG(tmp), gr_usr[1], gr_usr[3], gr_usr[2], gr_usr[4], interpolate=F, xpd = xpd)

points(x, pch = 19, col = adjustcolor(2,0.25), cex = 3)

# 
# 


#### clean example ####

#simulate observations
n <- 1E2
x <- matrix(rnorm(2*n), n, 2)

#create a blank plot in preferred graphical device
xpd <- NA
plot(NULL, ann=FALSE,
     xlim=quantile(x[,1], c(0.1,1)), ylim=quantile(x[,2], c(0.1,1)))
polygon(c(1,2,3), c(1,3,2), col = adjustcolor(1,0.5))

#snag plotting parameters
usr <- par("usr")
upct <- par("plt")
gr_usr <- usr + c(diff(usr[1:2]) / diff(upct[1:2]) * (c(0,1) - upct[1:2]),
                  diff(usr[3:4]) / diff(upct[3:4]) * (c(0,1) - upct[3:4]))

#generate a temporary file
tmp <- tempfile()
png(tmp, width = par("din")[1], height = par("din")[2], units = "in", res = 100, bg = "transparent")

#write to temporary file
par(mar = c(0,0,0,0), xpd = NA)
plot.new()
plot.window(xlim=gr_usr[1:2], ylim=gr_usr[3:4], xaxs = "i", yaxs = "i")
points(x, pch = 19, col = adjustcolor(1,0.5), cex = 3)
rect(usr[1], usr[3], usr[2], usr[4])
dev.off()

#read rasterized image back in and plot it
rasterImage(png::readPNG(tmp), gr_usr[1], gr_usr[3], gr_usr[2], gr_usr[4], interpolate=F, xpd = xpd)

#overlay points just to make sure we're doing this right
points(x, pch = 19, col = adjustcolor(2,0.25), cex = 3)
