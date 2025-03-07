library(svgparser)
library(concaveman)

# Example text
string <- "HelloWorld123"
font_name <- "sans"

# 1) Estimate device size in 'inches' based on strwidth/strheight (like in your PNG method)
#    Then scale to keep the max dimension at 5000 pixels worth of resolution.
str_rat <- c(
        w = strwidth(string, cex=1, units="inches", family=font_name),
        h = strheight(string, cex=1, units="inches", family=font_name)
)
str_rat <- str_rat / str_rat[1]  # ratio of w:h
npix_max <- 5000
npix <- round(str_rat / max(str_rat) * npix_max)

# 2) Create the SVG device
#    width/height in inches = npix / 72, so it's like a 72 DPI "canvas"
temp_svg <- tempfile(fileext = ".svg")
svg(filename = temp_svg, width = npix["w"]/72, height = npix["h"]/72)

par(mar = c(0,0,0,0))
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1), xaxs="i", yaxs="i")

# 3) Determine a cex that shrinks the text to about half the device space
est_dims <- c(
        w = strwidth(string, cex=1, units="user", family=font_name),
        h = strheight(string, cex=1, units="user", family=font_name)
)
max_cex <- 1 / max(est_dims)
cex2use <- 0.5 * max_cex
est_dims <- cex2use * est_dims

# 4) Plot text in SVG
text_center <- c(x=0.5, y=0.5)
text(x=text_center["x"], y=text_center["y"],
     labels=string, family=font_name,
     cex=cex2use, xpd=NA, col=1)

dev.off()  # close the SVG device

# 5) Parse the SVG with svgparser (data.frame mode)
svgdf <- read_svg(temp_svg, obj_type="data.frame")

# The device coordinates for all elements:
coords <- data.frame(svgdf$x, svgdf$y)
yr <- range(coords[,1])  # or range of x if we invert indexing
xr <- range(coords[,2])
coords <- coords[!(coords$svgdf.x %in% range(coords$svgdf.x)),]

# Some points might be background, line corners, etc.
# But if the text is converted into shapes, you'll see those x,y too.
# Just do a bounding polygon over everything:
bounding_poly <- concaveman(as.matrix(coords), concavity=2)

# We can replicate the logic from your PNG approach:
coords_fixed <- coords
coords_fixed[,2] <- (max(coords[,2]) - coords[,2])  # invert the y if needed

# Recompute bounding polygon on the flipped coords if you prefer that orientation
bpoly <- concaveman(as.matrix(coords_fixed), concavity=2)

# bounding box for the polygon
xr <- range(coords_fixed[,1])
yr <- range(coords_fixed[,2])
true_dims <- c(w=diff(xr), h=diff(yr))
true_loc <- c(x=mean(xr), y=mean(yr))

# We'll compare that to "est_dims" from earlier to see how the text scaling differs.
dim_scale <- true_dims / est_dims
loc_scale_disp <- (true_loc - 0.5) / true_dims

# Now we can replicate your step of re-plotting at a different location
# in a new device with the bounding polygon.

png(filename="~/test_svg_approach.png", units="px", width=5000, height=2500)
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

cex <- 40
est_wh <- c(
        w = strwidth(string, cex=cex, units="user", family=font_name),
        h = strheight(string, cex=cex, units="user", family=font_name)
)
true_wh <- est_wh * dim_scale
target_center <- c(x=0.5, y=0.25)
disp_center_by <- loc_scale_disp * true_wh
text_center2 <- target_center - disp_center_by

# Transform bounding_poly to user coords:
bpoly_centered <- t(bpoly[,c(1,2)])  # from coords_fixed
# shift it so the text center is at target_center
bpoly_centered <- bpoly_centered - true_loc
bpoly_centered <- bpoly_centered / true_dims * true_wh + target_center
bpoly_centered <- t(bpoly_centered)

# (A) Plot text with big cex
text(x=text_center2["x"], y=text_center2["y"],
     labels=string, family=font_name,
     cex=cex, xpd=NA, col=adjustcolor(1, 0.5))

# (B) Draw bounding polygon
lines(bpoly_centered, type="l", col="red", lwd=10)
points(target_center["x"],target_center["y"],col=2,pch=19)

dev.off()

