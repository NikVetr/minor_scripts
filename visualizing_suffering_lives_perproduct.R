makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}
boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA, 
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5, 
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font')){
  
  ## The Character expansion factro to be used:
  theCex <- graphics::par('cex')*cex
  
  ## Is y provided:
  if (missing(y)) y <- x
  
  ## Recycle coords if necessary:    
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]           
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }       
  }
  
  ## Width and height of text
  textHeight <- graphics::strheight(labels, cex = theCex, font = font)
  textWidth <- graphics::strwidth(labels, cex = theCex, font = font)
  
  ## Width of one character:
  charWidth <- graphics::strwidth("e", cex = theCex, font = font)
  
  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)            
    }        
  } else {
    adj <- c(0.5, 0.5)
  }
  
  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }       
  } else {
    offsetVec <- c(0, 0)
  }
  
  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }
  
  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]
  
  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth    
  graphics::rect(xleft = xMid - rectWidth/2, 
                 ybottom = yMid - rectHeight/2, 
                 xright = xMid + rectWidth/2, 
                 ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg)
  
  ## Place the text:
  graphics::text(xMid, yMid, labels, col = col.text, cex = theCex, font = font, 
                 adj = c(0.5, 0.5))    
  
  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                    yMid + rectHeight/2))
  }    
}
addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}
o_i_colors <- c(rgb(  0,   0,   0, maxColorValue = 255),  # black
                rgb(230, 159,   0, maxColorValue = 255),  # orange
                rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
                rgb(240, 228,  66, maxColorValue = 255),  # yellow
                rgb(  0, 158, 115, maxColorValue = 255),  # green
                rgb(  0, 114, 178, maxColorValue = 255),  # blue
                rgb(213,  94,   0, maxColorValue = 255),  # vermillion
                rgb(204, 121, 167, maxColorValue = 255))[-4]   # purple


recolor_image <- T
d <- read.csv("data/animalsuffering_lives_perproduct.csv")

colnames(d)[1] <- "animal"
d$total.lives.per.serving <- as.numeric(gsub(",", "", trimws(as.character(d$total.nanolives.per.serving)))) / 1E9
d$total.suffering.per.serving <- as.numeric(gsub(",", "", trimws(as.character(d$total.nanodays.suffering.per.serving)))) / 1E9
d$log.lives <- log10(d$total.lives.per.serving)
d$log.days <- log10(d$total.suffering.per.serving)

d <- d[order(d$log.days, decreasing = T),]
d$food <- trimws(as.character(d$food))
d$food[grep(pattern = "unbreaded pieces", d$food)]  <- "unbreaded pieces"
d$food[grepl(pattern = "breaded pieces", d$food) & !grepl(pattern = "unbread", d$food)]  <- "unbreaded pieces"
d$food <- sub(d$food, pattern = "/", replacement = " / ")

colmatch <- cbind(RColorBrewer::brewer.pal(max(as.numeric(as.factor(d$animal))), "Dark2"), trimws(levels(as.factor(d$animal))))
colmatch <- cbind(c("#D85128", "#00668E", "#F28A35", "#74a64d", "#676462", "#8C2C0E", "#871f78"), trimws(levels(as.factor(d$animal))))
colmatch <- cbind(ggthemes::colorblind_pal()(max(as.numeric(as.factor(d$animal)))), trimws(levels(as.factor(d$animal))))
colmatch <- cbind(o_i_colors, trimws(levels(as.factor(d$animal))))
#"beef"    "chicken" "dairy"   "egg"     "fish"    "pork"    "turkey" 

# colmatch[,1] <- sample(colmatch[,1])
# colmatch <- cbind(ggsci::pal_jco()(max(as.numeric(as.factor(d$animal)))), trimws(levels(as.factor(d$animal))))
# colmatch <- cbind(wesanderson::wes_palette("Zissou1", n = (max(as.numeric(as.factor(d$animal)))), type = "continuous"), trimws(levels(as.factor(d$animal))))
# colmatch <- cbind(inlmisc::GetTolColors(max(as.numeric(as.factor(d$animal))), scheme = "vibrant"), trimws(levels(as.factor(d$animal))))

#recolor the image
cols <- colmatch[,1][as.numeric(as.factor(d$animal))]
chkns <- png::readPNG("Downloads/rooster_and_hen_silhouette.png")

if(recolor_image){
  henColor_orig <- chkns[round(nrow(chkns) / 2), round(ncol(chkns) / 2),1:3] #plot(1,1,col = rgb(henColor_orig[1], henColor_orig[2], henColor_orig[3]), pch = 16, cex = 10)
  roosterColor_orig <- chkns[round(ncol(chkns) / 11 * 10), round(nrow(chkns) / 2.1),1:3] #plot(1,1,col = rgb(roosterColor_orig[1], roosterColor_orig[2], roosterColor_orig[3]), pch = 16, cex = 10)
  henColor_new <- as.vector(col2rgb(colmatch[colmatch[,2] == "egg", 1]) / 255)
  roosterColor_new <- as.vector(col2rgb(colmatch[colmatch[,2] == "chicken", 1]) / 255)
  for(i in 1:nrow(chkns)){
    if(i %% round(nrow(chkns) / 100) == 0){cat(paste0(i / round(nrow(chkns) / 100), "% "))}
    for(j in 1:ncol(chkns)){
      if(all(chkns[i,j,1:3] == henColor_orig)){
        chkns[i,j,1:3] <- henColor_new
      } else if(all(chkns[i,j,1:3] == roosterColor_orig)){
        chkns[i,j,1:3] <- roosterColor_new
      }
    }
  }
}
chkns_opacity <- 0.9
chkns_opacities <- chkns[,,4]
chkns_opacities[chkns_opacities == 1] <- chkns_opacity
chkns[,,4] <- chkns_opacities
  
png("Desktop/animal_suffering.png", width = 3840, height = 2060)
layout(mat = rbind(t(as.matrix(c(1,2,2))),
                   t(as.matrix(c(2,2,2)))))
layout(mat = rbind(t(as.matrix(c(1,1,1,2,2,2,2,2,2,2))),
                   t(as.matrix(c(1,1,1,2,2,2,2,2,2,2))),
                   t(as.matrix(c(2,2,2,2,2,2,2,2,2,2))),
                   t(as.matrix(c(2,2,2,2,2,2,2,2,2,2)))))

par(mar = c(0,8,6,0))
par(xpd=TRUE)

plot(d$log.days, d$log.lives, col = makeTransparent(cols, alpha = 0.75), pch = 16, cex = 5.5, xaxt = "n", yaxt = "n", bty = "n", 
     xlab = "", ylab = "", xlim = c(-4,1), ylim = c(-7,0))

axis(2, at = -7:0, labels = rep("", 8), line = -2, lwd = 3, lwd.ticks = 3, tck = -0.01, cex.axis = 2)
mtext(side = 2, at = -7:0, text = c(sapply(-7:-1, function(i) as.expression(bquote(10^ .(i)))), 1), cex = 2.5)

axis(3, at = -4:1, labels = rep("", 6), line = -2, lwd = 3, lwd.ticks = 3, tck = -0.01, cex.axis = 2)
mtext(side = 3, at = 1:-4, text = c(10, 1, sapply(-1:-4, function(i) as.expression(bquote(10^ .(i))))), cex = 2.5)
mtext(side = 3, at = -3.5, line = -7, text = "DAYS", font = 1, col = "darkred", cex = 3.5)
mtext(side = 2, at = -0.75, line = -7, text = "LIVES", font = 1, col = "darkred", cex = 3.5)

mtext(side = 1, at = -11:-6, text = c(10, 1, sapply(-1:-4, function(i) as.expression(bquote(10^ .(i))))), cex = 3)

par(mar = c(1,1,5,4))
par(xpd=TRUE)

plot(c(-12, 2.5), c(0,nrow(d)), xaxt = "n", yaxt = "n",type = "n", xlab = "", ylab = "", bty = "n")

lstart <- -6
rstart <- -4
for(i in 1:nrow(d)){
 
  #black border
  segments(x0 = lstart, x1 = lstart - d$log.days[i] - 4, y0 = i, y1 = i, lwd = 13, col = makeTransparent("black", alpha = 0.75)) 
  segments(x0 = rstart, x1 = rstart + d$log.lives[i] + 7, y0 = i, y1 = i, lwd = 13, col = makeTransparent("black", alpha = 0.75))
 
  segments(x0 = lstart, x1 = lstart - d$log.days[i] - 4, y0 = i, y1 = i, lwd = 10, col = makeTransparent(cols, alpha = 1)[i]) 
  segments(x0 = rstart, x1 = rstart + d$log.lives[i] + 7, y0 = i, y1 = i, lwd = 10, col = makeTransparent(cols, alpha = 1)[i]) 
  segments(x0 = lstart + 0.05, x1 = rstart - 0.05, y0 = i, y1 = i, lwd = 4, lty = 3, col = cols[i])
}
boxtext(x = -5, y = 1:nrow(d), labels = d$food, col = cols, cex = 4.5, col.bg = "white")
legend(x = "topright", legend = colmatch[,2], col = colmatch[,1], pch = 16, cex = 4.5, box.lty = 2, box.lwd = 3, pt.cex = 7, ncol = 4)
axis(side = 1, at = -12:3, labels = rep("", 16), line = -7, lwd = 5, lwd.ticks = c(rep(5, 7), 0, rep(5, 8)), tck = -0.01, cex.axis = 3)
axis(side = 1, at = log10(c(1E-4, c(sapply(2:-4, function(mag) c(10^mag*1:10)[-1])))), 
     labels = rep("", length(log10(c(1E-4, c(sapply(-4:2, function(mag) c(10^mag*1:10)[-1])))))), line = -7, lwd = 5, lwd.ticks = 3, tck = -0.005, cex.axis = 1)
axis(side = 1, at = -10 - log10(c(1E-4, c(sapply(1:-4, function(mag) c(10^mag*1:10)[-1])))), 
     labels = rep("", length(log10(c(1E-4, c(sapply(-4:1, function(mag) c(10^mag*1:10)[-1])))))), line = -7, lwd = 5, lwd.ticks = 3, tck = -0.005, cex.axis = 1)
#-11:-6 on left, -4:3 on right
mtext(side = 1, at = -12:-6, text = c(100, 10, 1, sapply(-1:-4, function(i) as.expression(bquote(10^ .(i))))), cex = 3, line = -2)
mtext(side = 1, at = -4:3, text = c(sapply(-7:-1, function(i) as.expression(bquote(10^ .(i)))), 1), cex = 3, line = -2)
mtext(side = 1, at = -5, text = "I", cex = 6.285, col = "white", line = -2)
text(x = lstart - 0.35, y = nrow(d) + 3.75, labels = "DAYS", cex = 8, font = 2, col = "darkred")
text(x = lstart + 0.45, y = nrow(d) + 4.15, labels = "(of suffering)", cex = 4, font = 2, col = "darkred")
text(x = rstart + 0.35, y = nrow(d) + 3.75, labels = "LIVES", cex = 8, font = 2, col = "darkred")
text(x = lstart - 0.45, y = nrow(d) + 1.45, labels = "per serving", cex = 5.5, font = 1, col = "black")
text(x = rstart + 0.425, y = nrow(d) + 1.45, labels = "per serving", cex = 5.5, font = 1, col = "black")
text(x = 0.4, y = 82.5, labels = "data source:", cex = 3.5, font = 1, col = "black")
text(x = 1.9, y = 82.375, labels = "faunalytics.org/animal-product-impact-scales/", cex = 3.5, font = 1, col = "blue")

#add picture of chickens
addImg(chkns, x = -11.55, y = 5.5, width = 1.25)
# axis(side = 1, labels = c(2*2^(-(0:3)), sapply(-(3:9), function(i) as.expression(bquote(2^ .(i))))), at = log10(2*2^(-(0:10))), las = 2)

dev.off()

