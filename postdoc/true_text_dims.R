# set.seed(1)
string <- paste0(sample(c(LETTERS, letters, 0:9), 12), collapse = "")
font_name <- "Source Serif Pro"

str_rat <- c(w = strwidth(cex = 1, string, units = "inches", family = font_name),
             h = strheight(cex = 1, string, units = "inches", family = font_name))
str_rat <- str_rat / str_rat[1]
npix_max <- 1E3
npix <- round(str_rat / max(str_rat) * npix_max)

#generate plot
temp_png <- paste0(tempfile(), ".png")
png(filename = temp_png, units = "px", width = npix["w"], height = npix["h"])
par(mar = c(0,0,0,0))
plot.new()
plot.window(xlim = c(0, 1),
            ylim = c(0, 1),
            xaxs="i", yaxs="i")

#initial dimension estimate
est_dims <- c(w = strwidth(string, cex = 1, units = "user", family = font_name),
              h = strheight(string, cex = 1, units = "user", family = font_name))
max_cex <- 1 / max(est_dims)
cex2use <- max_cex * 0.5
est_dims <- cex2use * est_dims

#plot the text
text(x = 0.5, 
     y = 0.5, 
     labels = string, family = font_name, 
     cex = cex2use, xpd = NA, col = 1)

#close the window and device device
dev.off()
dev.off()

arr <- png::readPNG(temp_png)
mat <- round(apply(arr, c(1,2), sum) / 4)
coords <- which(mat==0, arr.ind = T)
yr <- rev(dim(mat)[1] - range(coords[,1]))
xr <- range(coords[,2])
xc <- xr / npix[1]
yc <- yr / npix[2]
true_dims <- c(w = diff(xc), h = diff(yc))
est_dims
true_dims
true_loc <- c(x = mean(xc), y = mean(yc))
true_loc
