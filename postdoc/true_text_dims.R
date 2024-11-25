# set.seed(1)
string <- paste0(sample(c(LETTERS, letters, 0:9), 12), collapse = "")
font_name <- "Source Serif Pro"

using_png <- F
str_rat <- c(w = strwidth(cex = 1, string, units = "inches", family = font_name),
             h = strheight(cex = 1, string, units = "inches", family = font_name))
str_rat <- str_rat / str_rat[1]
npix_max <- 5E3
npix <- round(str_rat / max(str_rat) * npix_max)

#generate plot
temp_filename <- paste0(tempfile(), ifelse(using_png, ".png", ".svg"))
if(using_png){
        png(filename = temp_filename, units = "px", width = npix["w"], height = npix["h"])
} else {
        svg(filename = temp_filename, width = npix["w"] / 72, height = npix["h"] / 72)
        
}

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
text_center <- c(x = 0.5, y = 0.5)
text(x = text_center["x"], 
     y = text_center["y"], 
     labels = string, family = font_name, 
     cex = cex2use, xpd = NA, col = 1)

#close the window and device device
dev.off()
dev.off()

arr <- png::readPNG(temp_filename)
mat <- round(apply(arr, c(1,2), sum) / 4)
coords <- which(mat==0, arr.ind = T)
coords[,1] <- (dim(mat)[1] - coords[,1]) / npix[2]
coords[,2] <- coords[,2] / npix[1]
yr <- range(coords[,1])
xr <- range(coords[,2])
xc <- xr
yc <- yr
true_dims <- c(w = diff(xc), h = diff(yc))
dim_scale <- true_dims / est_dims
true_loc <- c(x = mean(xc), y = mean(yc))
loc_scale_disp <- (true_loc - 0.5) / true_dims
bounding_poly <- concaveman::concaveman(coords, concavity = 2)


#### plot a quick test of converted scale ####

#note -- inter-letter kerning differs between pdf and png devices
png(filename = "~/test.png", units = "px", 
    width = 5000, height = 2500)

# cairo_pdf(filename = "~/test.pdf", 
#     width = 5000/72, height = 2500/72)

cex <- 40
par(mar = c(0,0,0,0))
plot.new()
plot.window(xlim = c(0, 1),
            ylim = c(0, 1),
            xaxs="i", yaxs="i")
est_wh <- c(w = strwidth(string, cex = cex, 
                         units = "user", family = font_name),
            h = strheight(string, cex = cex, 
                          units = "user", family = font_name))
true_wh <- est_wh * dim_scale
disp_center_by <- loc_scale_disp * true_wh
target_center <- c(x = 0.5, y = 0.25)
text_center <- target_center - disp_center_by
bpoly <- t(bounding_poly[,2:1]) - true_loc
bpoly <- t(bpoly / true_dims * true_wh + target_center) 

#plot text
text(x = text_center["x"], 
     y = text_center["y"], 
     labels = string, family = font_name, 
     cex = cex, xpd = NA, col = adjustcolor(1, 0.5))
# lines(bounding_poly[,2:1], type = "l", col = 2)
lines(bpoly, type = "l", col = 2, lwd = 10)
points(target_center["x"],target_center["y"],col=2,pch=19)
dev.off()
