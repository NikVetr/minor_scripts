
grad_arrow_curve <- function(arrow_locs, prop_head_width = 2, prop_shaft_length = 0.1,
                             cols = c("white", "black"), nslices = 500, direc = "h", w = 0.2,
                             raster = F, xpd = NA, raster_res = 150, interp_raster = T,
                             col_alpha = 1, col_pow = 0.5, taper_ratio = 0.5, taper_pow = 2){
  
  if(direc == "h"){
    dx <- diff(arrow_locs[1:2]) * (1-prop_shaft_length)
    dy <- diff(arrow_locs[3:4])
  } else if(direc == "v"){
    dx <- diff(arrow_locs[1:2])
    dy <- diff(arrow_locs[3:4]) * (1-prop_shaft_length)
  }
  
  colsgrad <- colorRampPalette(colors = cols)(nslices * 10)
  colsgrad <- colsgrad[round((seq(1, (nslices*10)^col_pow, length.out=nslices))^(1/col_pow))]
  colsgrad <- adjustcolor(colsgrad, col_alpha)
  
  # split the base of the arrow in nslices and fill each progressively
  piseq <- seq(-pi/2,pi/2,length.out=nslices+1)
  
  if(direc == "h"){
    xs <- seq(arrow_locs[1],arrow_locs[1]+dx, length.out=nslices+1)
    ys <- dy*(sin(piseq)+1)/2 + arrow_locs[3]
    m <- cos(piseq)/2 * sign(dy) * sign(dx)
    t <- atan(m)
    dispx <- w * sin(t) / 2
    dispy <- w * cos(t) / 2
  } else if (direc == "v"){
    ys <- seq(arrow_locs[3],arrow_locs[3]+dy, length.out=nslices+1)
    xs <- dx*(sin(piseq)+1)/2 + arrow_locs[1]
    m <- cos(piseq)/2 * sign(dy) * sign(dx)
    t <- atan(m)
    dispx <- w * cos(t) / 2
    dispy <- w * sin(t) / 2
  }
  
  #taper the arrow if desired
  taper <- seq(1, taper_ratio^taper_pow, length.out=nslices+1)^(1/taper_pow)
  dispx <- dispx * taper
  dispy <- dispy * taper
  
  #final coords
  coords <- data.frame(x1 = xs - dispx, y1 = ys + dispy,
                       x2 = xs + dispx, y2 = ys - dispy
  )
  
  
  if(raster){
    
    #get plotting params
    usr <- par("usr")
    upct <- par("plt")
    gr_usr <- usr + c(diff(usr[1:2]) / diff(upct[1:2]) * (c(0,1) - upct[1:2]),
                      diff(usr[3:4]) / diff(upct[3:4]) * (c(0,1) - upct[3:4]))
    
    #write to temporary png
    tmp <- tempfile()
    
    png(tmp, width = par("din")[1], height = par("din")[2], units = "in", res = raster_res, bg = "transparent", type="cairo")
    par(mar = c(0,0,0,0), xpd = NA)
    plot.new(); plot.window(xlim=gr_usr[1:2], ylim=gr_usr[3:4], xaxs = "i", yaxs = "i")
    
    #plot the arrows
    if(length(cols) == 1){
      polygon(x = c(coords$x1, rev(coords$x2)), y = c(coords$y1, rev(coords$y2)), cols = col, border = NA)  
    } else {
      for(i in 1:nslices){
        polygon(x = c(coords$x1[i], coords$x1[i+1], coords$x2[i+1], coords$x2[i]), 
                y = c(coords$y1[i], coords$y1[i+1], coords$y2[i+1], coords$y2[i]), 
                col = colsgrad[i], border = NA) 
      }
    }
    if(direc == "h"){
      polygon(x = c(xs[nslices], xs[nslices], arrow_locs[2]), 
              y = c(ys[nslices] - prop_head_width * w / 2 * taper_ratio, ys[nslices] + prop_head_width * w / 2 * taper_ratio, ys[nslices]),
              col = colsgrad[nslices], border = NA)
    } else if(direc == "v"){
      polygon(x = c(xs[nslices] - prop_head_width * w / 2 * taper_ratio, xs[nslices] + prop_head_width * w / 2 * taper_ratio, xs[nslices]), 
              y = c(ys[nslices], ys[nslices], arrow_locs[4]),
              col = colsgrad[nslices], border = NA)
    }
    dev.off()
    
    #draw to file
    rasterImage(png::readPNG(tmp), gr_usr[1], gr_usr[3], gr_usr[2], gr_usr[4], interpolate = interp_raster, xpd = xpd)
    
    #delete temp file
    rm(tmp)
    
  } else {
    
    #or alternatively draw it in full vector graphics
    if(length(cols) == 1){
      polygon(x = c(coords$x1, rev(coords$x2)), y = c(coords$y1, rev(coords$y2)), cols = col, border = NA, xpd = xpd)  
    } else {
      for(i in 1:nslices){
        polygon(x = c(coords$x1[i], coords$x1[i+1], coords$x2[i+1], coords$x2[i]), 
                y = c(coords$y1[i], coords$y1[i+1], coords$y2[i+1], coords$y2[i]), 
                col = colsgrad[i], border = NA, xpd = xpd) 
      }
    }
    if(direc == "h"){
      polygon(x = c(xs[nslices], xs[nslices], arrow_locs[2]), 
              y = c(ys[nslices] - prop_head_width * w / 2 * taper_ratio, ys[nslices] + prop_head_width * w / 2 * taper_ratio, ys[nslices]),
              col = colsgrad[nslices], border = NA)
    } else if(direc == "v"){
      polygon(x = c(xs[nslices] - prop_head_width * w / 2 * taper_ratio, xs[nslices] + prop_head_width * w / 2 * taper_ratio, xs[nslices]), 
              y = c(ys[nslices], ys[nslices], arrow_locs[4]),
              col = colsgrad[nslices], border = NA)
    }
  }
  
}

plot(NULL, xaxt="n",yaxt="n",bty="n",pch="",ylab="",xlab="", main="", sub="", 
     xlim = c(0,9), ylim = c(0.5,3.5))
grad_arrow_curve(c(1,5,1,3))
grad_arrow_curve(arrow_locs = c(1,5,3,1))
grad_arrow_curve(c(6,2,3,1))
grad_arrow_curve(c(6,2,1,3))
