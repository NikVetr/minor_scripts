#packages
library(png)
library(magick)

#functions


grad_arrow_curve <- function(arrow_locs, prop_head_width = 2, prop_shaft_length = 0.1,
                             cols = c("white", "black"), nslices = 500, direc = "h", w = 0.25,
                             raster = T, xpd = NA, raster_res = 51, interp_raster = T,
                             col_alpha = 1, col_pow = 0.5, taper_ratio = 0.45, taper_pow = 1.5,
                             outline = T, outline_lwd = 2.5, outline_col = "black"){
  
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
  if(direc == "h"){
    head_coords <- data.frame(x = rev(c(xs[nslices], arrow_locs[2], xs[nslices])), 
                              y = rev(c(ys[nslices] - prop_head_width * w / 2 * taper_ratio, ys[nslices], ys[nslices] + prop_head_width * w / 2 * taper_ratio)))
  } else if(direc == "v"){
    head_coords <- data.frame(x = c(xs[nslices] - prop_head_width * w / 2 * taper_ratio, xs[nslices], xs[nslices] + prop_head_width * w / 2 * taper_ratio), 
                              y = c(ys[nslices], arrow_locs[4], ys[nslices]))
  }
  
  
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
      polygon(x = c(coords$x1, rev(coords$x2)), y = c(coords$y1, rev(coords$y2)), cols = cols, border = NA)  
    } else {
      for(i in 1:nslices){
        polygon(x = c(coords$x1[i], coords$x1[i+1], coords$x2[i+1], coords$x2[i]), 
                y = c(coords$y1[i], coords$y1[i+1], coords$y2[i+1], coords$y2[i]), 
                col = colsgrad[i], border = NA) 
      }
    }
    if(direc == "h"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    } else if(direc == "v"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    }
    dev.off()
    
    #draw to file
    rasterImage(png::readPNG(tmp), gr_usr[1], gr_usr[3], gr_usr[2], gr_usr[4], interpolate = interp_raster, xpd = xpd)
    
    #delete temp file
    rm(tmp)
    
  } else {
    
    #or alternatively draw it in full vector graphics
    if(length(cols) == 1){
      polygon(x = c(coords$x1, rev(coords$x2)), y = c(coords$y1, rev(coords$y2)), cols = cols, border = NA, xpd = xpd)  
    } else {
      for(i in 1:nslices){
        polygon(x = c(coords$x1[i], coords$x1[i+1], coords$x2[i+1], coords$x2[i]), 
                y = c(coords$y1[i], coords$y1[i+1], coords$y2[i+1], coords$y2[i]), 
                col = colsgrad[i], border = NA, xpd = xpd) 
      }
    }
    if(direc == "h"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    } else if(direc == "v"){
      polygon(x = head_coords$x, y = head_coords$y, col = colsgrad[nslices], border = NA)
    }
  }
  
  if(outline){
    polygon(x = c(coords$x1, head_coords$x, rev(coords$x2)), 
            y = c(coords$y1, head_coords$y, rev(coords$y2)), 
            border = outline_col, xpd = xpd, lwd = outline_lwd)
  }
  
}

inverse_hex_color <- function(hex_color) {
  # Convert the hex color to RGB
  rgb_color <- col2rgb(hex_color)
  # Invert each RGB channel
  inverted_rgb <- bitwXor(rgb_color, 255)
  # Convert the inverted RGB color back to hex
  inverted_hex <- rgb(inverted_rgb[1], inverted_rgb[2], inverted_rgb[3], maxColorValue = 255)
  return(inverted_hex)
}

find_factors <- function(A) {
  sqrt_A <- sqrt(A)
  B <- ceiling(sqrt_A)
  while (A %% B != 0) {
    B <- B + 1
  }
  C <- A / B
  return(c(B = B, C = C))
}

dexp_smooth <- function(x, y, r, cyclic = T){
  n <- length(x)
  w <- t(sapply(1:(n-1), function(i) c(dexp((x[i] - x[0:(i-1)]), rate = r), 
                                       dexp(0, rate = r), dexp(-(x[i] - x[(i+1):n]), rate = r))))
  w <- rbind(w, dexp(x[n] - x, rate = r))
  wy <- c(w %*% t(t(y)))
  wy <- wy / apply(w,1,sum)
  return(wy)
}

interpolate <- function(p1, p2) {
  max.length = max(abs(unlist(c(p1[1] - p2[1], p1[2] - p2[2]))))
  return(cbind(x = seq(unlist(p1[1]), unlist(p2[1]), length.out = max.length+1),
               y = seq(unlist(p1[2]), unlist(p2[2]), length.out = max.length+1)))
}

wrap <- function(x, mx){
  xw <- x %% mx
  xw[xw == 0] <- mx
  xw
}

remove_points <- function(pts, rpts){
  
  if(is.null(nrow(pts)) & length(pts) != 0){
    pts <- matrix(pts, nrow = 1)
  }
  
  if(is.null(nrow(rpts)) & length(rpts) != 0){
    rpts <- matrix(rpts, nrow = 1)
  }
  
  if(nrow(pts) == 0 | nrow(rpts) == 0){
    return(pts)
  }
  
  if(nrow(rpts) == 1){
    pts[!apply(pts, 1, function(ri) all(ri == rpts)),]
  } else {
    remove_points(pts[!apply(pts, 1, function(ri) all(ri == rpts[1,])),], rpts[-1,])
  }
  
}
# remove_points(pts, rpts)

remove_points_no_dupes <- function(pts, rpts){
  all_mat <- rbind(as.matrix(rpts), as.matrix(pts))
  pts[!duplicated(all_mat)[-(1:nrow(rpts))],]
}

find_valid_interior_point <- function(shape, si){
  
  #basic constants
  dirs <- t(expand.grid(-1:1, -1:1))
  dirs <- dirs[,!apply(dirs == c(0,0), 2, all)]
  nesw <- cbind(c(0,1), c(1,0), c(0,-1), c(-1,0))
  
  #initialize navigation list
  interior_pts <- list()
  
  #first find *one* point inside the polygon
  found_first_point <- F
  while(!found_first_point){
    if(si >= nrow(shape)){return(NULL)}
    adj_pts <- unlist(shape[si,]) + dirs
    adj_pt_vals <- sp::point.in.polygon(point.x = adj_pts[1,], point.y = adj_pts[2,],
                                        pol.x = shape$x, pol.y = shape$y)
    proposed_inter_points_inds <- which(adj_pt_vals == 1)
    if(length(prop_inter_points_inds) != 0){
      #first point also needs to have somewhere to go, not just be interior to the shape
      proposed_inter_pts <- as.matrix(adj_pts[,proposed_inter_points_inds])
      nearby_line_pts <- shape[wrap(si + -20:20, nrow(shape)),]
      proposed_points_on_line <- apply(proposed_inter_pts, 2, function(prop_pt) 
        sum(duplicated(rbind(t(c(prop_pt) + nesw), as.matrix(nearby_line_pts)))) == 4)
      
      if(all(proposed_points_on_line)){
        si <- si + 1
      } else {
        found_first_point <- T
        inter_dir <- min(which(!proposed_points_on_line))
        interior_pts[[1]] <- proposed_inter_pts[,inter_dir] #x,y
        previously_visited <- rbind(proposed_inter_pts[,inter_dir])  
      }
    } else {
      si <- si + 1
    }
  }
  
  return(list(si = si, interior_pts = interior_pts, previously_visited = previously_visited))
}

find_interior_shape <- function(shape, plot_result = F, debug = F){
  
  #basic constants
  dirs <- t(expand.grid(-1:1, -1:1))
  dirs <- dirs[,!apply(dirs == c(0,0), 2, all)]
  nesw <- cbind(c(0,1), c(1,0), c(0,-1), c(-1,0))
  
  #find our first point
  isi <- i <- 1
  first_point <- find_valid_interior_point(shape, 1)
  if(is.null(first_point)){return(NULL)}
  interior_pts <- first_point$interior_pts
  si <- first_point$si
  previously_visited <- first_point$previously_visited
  
  #iterate over outer polygon pts, following it on the inside
  returned_to_start <- F
  max_iter <- round(3 * nrow(shape))
  buffer <- 0
  max_buffer <- 5
  buffer_reset <- 2 
  previously_visited_length <- 4
  pt_of_return <- T
  isi <- isi + 1
  last_buffer_reset_point <- c(Inf, Inf)
  n_opportunities_to_finish <- 5
  shape_dists <- as.matrix(dist(shape))
  rerun_finish_opportunities <- T
  
  #record locs
  chain_history <- data.frame(matrix(NA, nrow = max_iter, ncol = 8, 
                                     dimnames = list(1:max_iter, c("i", "si", "isi", "sx", "sy", "isx", "isy", "keep"))))
  chain_history$i <- 1:max_iter
  chain_history[i,c("si", "isi")] <- c(si-1, isi-1)
  chain_history[i,c("sx", "sy")] <- shape[wrap(si-1, nrow(shape)),]
  chain_history[i,c("isx", "isy")] <- interior_pts[[isi-1]]
  chain_history$keep <- F
  chain_history$keep[1] <- T
  
  bad_ilocs <- numeric(0)
  while(!returned_to_start & i < max_iter){

    i <- i + 1
    si <- si + 1
    
    if(si > nrow(shape) | si < 1){
      si <- wrap(si, nrow(shape))
    }
      
    #all of these are either interior pts or on the border
    adj_inter_pts_im1 <- t(interior_pts[[isi-1]] + nesw)
    
    #find point around and remove points on outer line from candidacy
    focal_shape_inds <- wrap(si + -(5+buffer):(5+buffer), nrow(shape))
    #this finds everything in the neighborhood of the neighborhood of our point
    focal_shape_inds <- sort(unique(c(focal_shape_inds, 
                                     which(apply(shape_dists[focal_shape_inds,], 2, 
                                                 function(fsi_dists) any(fsi_dists < 4.1))))))
    pts_on_line <- shape[focal_shape_inds,]
    
    
    #points adjacent to border
    adj_pts_i <- do.call(rbind, lapply(1:nrow(pts_on_line), function(pti) t(unlist(pts_on_line[pti,]) + dirs)))
    adj_pts_i <- adj_pts_i[!duplicated(adj_pts_i),]
    adj_pts_i <- rbind(as.matrix(pts_on_line), adj_pts_i)    

    adj_pts_i <- adj_pts_i[-(1:nrow(pts_on_line)),][!duplicated(adj_pts_i)[-(1:nrow(pts_on_line))],]
    combined <- rbind(adj_inter_pts_im1, adj_pts_i)
    
    #remove points already visited in interior shape
    combined <- remove_points(combined, previously_visited)
    
    #remove points that are "bad" ie in corners where you get stuck
    if(length(bad_ilocs) > 0){
      combined <- remove_points(combined, bad_ilocs)  
    }
    
    #find our next point and record it
    proposed_next_point <- combined[duplicated(combined),]
    if(!is.null(nrow(proposed_next_point))){
      if(nrow(proposed_next_point) > 1){
        proposed_next_point <- proposed_next_point[1,] #sometimes more than one direction to go  
      }
    }
    
    #check for legitimacy
    next_interior_pt <- proposed_next_point
    if(length(next_interior_pt) == 0){ #sometimes the next point is empty bc movement around is smaller, reuse previous pt if so
      
      buffer <- buffer + 1
      
    } else {
      
      #find distances of this new interior point to existing points
      nearby_interior_pts <- which(apply(abs(next_interior_pt - t(chain_history[chain_history$keep, c("isx", "isy")])), 2, sum) < 1.1)
      previously_visited <- chain_history[as.numeric(names(nearby_interior_pts)), c("isx", "isy")]
      
      #add to running list
      interior_pts[[isi]] <- next_interior_pt
      chain_history$keep[i] <- T
      isi <- isi + 1
      
    }
    
    #if buffer has grown large, things can be inefficient
    #so find closest point to reset it
    if(buffer > max_buffer){
      #if you hit the buffer within (max_buffer - buffer_reset) moves
      #you are stuck in a corner and need to back out (bc you haven't moved anywhere)
      #so go back a point
      if(all(last_buffer_reset_point == interior_pts[[isi-1]])){
        
        #check we're not just stuck back at the start or at first opportunity
        if(isi > (n_opportunities_to_finish + 5)){
          
          #find the first opportunity to exist the while loop
          if(rerun_finish_opportunities){
            n_opportunities_to_finish <- 5
            isi_to_check <- 1
            opportunities_to_finish <- list()
            while(length(opportunities_to_finish) < n_opportunities_to_finish & isi_to_check < length(interior_pts)){
              poss_neighbors <- interior_pts[[isi_to_check]] + nesw
              
              poss_innershape_pts <- chain_history$isi[chain_history$keep][apply(abs(interior_pts[[isi_to_check]] - 
                                            t(chain_history[chain_history$keep, c("isx", "isy")])), 2, sum) == 1]
              poss_innershape_pts <- setdiff(poss_innershape_pts, isi-1)
              poss_innershape_pts <- do.call(rbind, interior_pts[poss_innershape_pts])
              poss_neighbors <- remove_points(pts = t(poss_neighbors), rpts =  poss_innershape_pts)
              
              poss_outershape_pts <- shape[apply(abs(interior_pts[[isi_to_check]] - t(shape)), 2, sum) == 1,]
              poss_neighbors <- remove_points(pts = poss_neighbors, rpts = poss_outershape_pts)
              
              if(length(poss_neighbors) > 0.1){
                opportunities_to_finish[[length(opportunities_to_finish) + 1]] <- isi_to_check
              }
              
              isi_to_check <- isi_to_check + 1  
            }
            
          }
          opportunities_to_finish <- unlist(opportunities_to_finish)
          rerun_finish_opportunities <- F
          
          #ORRR we could allow the chain to fully retract, 
          #and then take progessively longer chains 
          #and check for their termination?
          
          #maybe we could clarify that the termination is from having run into oneself?
          
          dist_from_return <- unlist(lapply(interior_pts[opportunities_to_finish], 
                                            function(point) sum(abs(point - interior_pts[[isi-1]]))))
          pt_of_return <- dist_from_return < 1.1
          returned_to_start <- any(pt_of_return)
        }
        
        #if chain is not able to make it around a corner with this buffer length
        if(!returned_to_start){
          bad_ilocs <- rbind(bad_ilocs, interior_pts[[isi-1]])
          interior_pts <- interior_pts[-(isi-1)]
          isi <- isi - 1
          chain_history$keep[tail(which(chain_history$keep), 1)] <- F
        }
        
        #if we have returned to start, roll a new starting point from a 
        #position larger than all of the preceding si
        if(length(interior_pts) == 0){
          
          new_first_point <- find_valid_interior_point(shape, max(chain_history$si, na.rm = T)+1)
          if(is.null(new_first_point)){return(NULL)}
          interior_pts <- new_first_point$interior_pts
          isi <- 2
          si <- new_first_point$si
          previously_visited <- new_first_point$previously_visited
          chain_history$keep <- F
          chain_history$keep[i] <- T
          rerun_finish_opportunities <- T
          
        } else {
          
          #this gets us to the closest si to our previous location, but maybe we want 
          #the closest period? or maybe only do if si is far away from isi
          hist_dists_is <- chain_history$isi[!is.na(chain_history$isi)] - isi
          history_reset_ind <- max(which(abs(hist_dists_is) == min(abs(hist_dists_is))))
          si <- chain_history$si[history_reset_ind]
          
          #need to also reset the previously_visited variable
          nearby_interior_pts <- which(apply(abs(interior_pts[[isi-1]] - t(chain_history[chain_history$keep, c("isx", "isy")])), 2, sum) < 1.1)
          previously_visited <- chain_history[as.numeric(names(nearby_interior_pts)), c("isx", "isy")]
                  
        }
        
        #equivalent to just snipping off pouches in the outer shape that have
        #length > (max_buffer - buffer_reset)
        
      }
      
      #find relevant shape points
      if(sum(abs(shape[si,] - interior_pts[[isi-1]])) > 5){
        #use entire shape
        si_inds <- 1:nrow(shape)
      } else {
       #find closest points on outer shape within buffer
        si_inds <- wrap(si + -(2+buffer):(2+buffer), nrow(shape))
      }
      pts_on_shape <- shape[si_inds,]
      
      #reset si to that closest point
      dist_pts <- apply((t(pts_on_shape) - interior_pts[[isi-1]])^2, 2, sum)
      min_dist_pt <- max(which(dist_pts == min(dist_pts)))
      si <- wrap(si_inds[min_dist_pt], nrow(shape))
      
      #reset buffer and record when you did it
      buffer <- buffer_reset
      last_buffer_reset_point <- interior_pts[[isi-1]]

    }
    
    #record history
    chain_history[i,c("si", "isi")] <- c(si-1, isi-1)
    chain_history[i,c("sx", "sy")] <- shape[wrap(si-1, nrow(shape)),]
    chain_history[i,c("isx", "isy")] <- interior_pts[[isi-1]]
    
    if(debug){
      #debug mode -- laboriously check each new point for internality to shape
      if(sp::point.in.polygon(point.x = interior_pts[[isi-1]][1], point.y = interior_pts[[isi-1]][2],
                           pol.x = shape$x, pol.y = shape$y) != 1){
        break
      }
    }
    
    }
  
  inter_pts <- do.call(rbind, interior_pts)
  
  #only keep up to where we circle around back to start
  # 
  if(any(pt_of_return)){
    inter_pts <- inter_pts[min(opportunities_to_finish[which(pt_of_return)]):nrow(inter_pts),]  
  } else {
    inter_pts <- chain_history[chain_history$keep,c("isx", "isy")]
  }
  
  if(debug){
    #check in progress
    par(mfrow = c(2,1))
    curr_inter_pts <- do.call(rbind, interior_pts)
    xr <- range(tail(curr_inter_pts[,1], 5)) + c(-4,4)
    yr <- range(tail(curr_inter_pts[,2])) + c(-4,4)
    plot(rbind(shape, shape[1,]), type = "l", 
         xlim = xr, ylim = yr)
    abline(v=xr[1]:xr[2], lty = 3, col = adjustcolor(1, 0.5), lwd = 0.75)
    abline(h=yr[1]:yr[2], lty = 3, col = adjustcolor(1, 0.5), lwd = 0.75)
    lines(curr_inter_pts, col = 2, lty = 2)
    
    #label interior points
    text(x = curr_inter_pts[,1] + 0.2, curr_inter_pts[,2] + 0.2, labels = 1:nrow(curr_inter_pts), col = 2)
    segments(x0 = curr_inter_pts[,1], x1 = curr_inter_pts[,1] + 0.15, 
             y0 = curr_inter_pts[,2], y1 = curr_inter_pts[,2] + 0.15, lty = 3, col = 2)
    
    #label outer points
    text(x = shape[,1] + 0.2, shape[,2] + 0.2, labels = 1:nrow(shape))
    segments(x0 = shape[,1], x1 = shape[1:nrow(shape),1] + 0.15, 
             y0 = shape[,2], y1 = shape[,2] + 0.15, lty = 3)
    
    if(!is.null(dim(interior_pts[[isi-1]]))){
      points(interior_pts[[isi-1]][,1], interior_pts[[isi-1]][,2], pch = 19)
    } else {
      points(interior_pts[[isi-1]][1], interior_pts[[isi-1]][2], pch = 19)  
    }
    
    points(adj_inter_pts_im1, pch = 19, col = adjustcolor(2, 0.2), cex = 2)
    points(adj_pts_i, pch = 19, cex = 3, col = adjustcolor(4, 0.2))
    points(bad_ilocs, pch = 19, col = 2)
    
  }
  
  #check final output
  if(plot_result | debug){
    plot(rbind(shape, shape[1,]), type = "l", col = 2, lwd = 2, xlim = range(c(shape[,1], inter_pts[,1])),
         ylim = range(c(shape[,2], inter_pts[,2])))
    abline(v=floor(par("usr")[1]):ceiling(par("usr")[2]), lty = 3, col = adjustcolor(1, 0.5), lwd = 0.75)
    abline(h=floor(par("usr")[3]):ceiling(par("usr")[4]), lty = 3, col = adjustcolor(1, 0.5), lwd = 0.75)
    polygon(inter_pts, col = adjustcolor(3, 0.25))
    points(bad_ilocs, pch = 19, col = 2)
  }
  
  return(data.frame(x = inter_pts[,1], y = inter_pts[,2]))
  
}


find_interior_shape_2 <- function(shape, plot_result){
  
  #compass constants
  dirs <- t(expand.grid(-1:1, -1:1))
  dirs <- dirs[,!apply(dirs == c(0,0), 2, all)]
  nesw = cbind(c(0,1), c(1,0), c(0,-1), c(-1,0))
  
  #find a point inside the shape
  first_point <- find_valid_interior_point(shape, 1)$interior_pts[[1]]
  
  #if the polygon is a square
  max_pts <- (ceiling(nrow(shape) / 4) - 1)^2
  inner_pts <- data.frame(matrix(NA, ncol = 4, nrow = max_pts, dimnames = list(1:max_pts, c("x", "y", "f", "r"))))
  inner_pts$f <- F
  inner_pts$r <- F
  
  npts <- 1
  inner_pts[npts, c("x", "y")] <- first_point
  inner_pts$f[npts] <- T
  
  #flood fill from that first point
  i <- 2
  focal_i <- 1
  focal_point <- unlist(inner_pts[focal_i, c("x", "y")])
  
  possible_points <- t(focal_point + nesw)
  possible_points <- remove_points_no_dupes(possible_points, shape)
  possible_points <- remove_points_no_dupes(possible_points, inner_pts[inner_pts$f, c("x", "y")])
  possible_points <- matrix(possible_points, ncol = 2)
  npts <- nrow(possible_points)
  
  inner_pts[i:(i+npts-1), c("x", "y")] <- possible_points
  inner_pts$f[i:(i+npts-1)] <- T
  inner_pts$r[focal_i] <- T
  i <- i + npts
  
  while(sum(inner_pts$r) != sum(inner_pts$f)){
    focal_i <- setdiff(which(inner_pts$f), which(inner_pts$r))[1]
    focal_point <- unlist(inner_pts[focal_i, c("x", "y")])
    
    possible_points <- t(focal_point + nesw)
    possible_points <- remove_points_no_dupes(possible_points, shape)
    possible_points <- remove_points_no_dupes(possible_points, inner_pts[inner_pts$f, c("x", "y")])
    possible_points <- matrix(possible_points, ncol = 2)
    npts <- nrow(possible_points)
    inner_pts$r[focal_i] <- T
    
    if(npts > 0.1){
      inner_pts[i:(i+npts-1), c("x", "y")] <- possible_points
      inner_pts$f[i:(i+npts-1)] <- T
      i <- i + npts  
    }
    
  }
  
  inner_pts <- inner_pts[inner_pts$f,c("x", "y")]
  
  tshape <- t(shape)
  on_border <- apply(inner_pts, 1, function(ipi){
    dists <- abs(tshape - unlist(ipi))
    any(dists[1,] + dists[2,] == 1)
    any(dists[1,]^2 + dists[2,]^2 < 2.1)
    })
  border_pts <- inner_pts[on_border,]
  finterior_border_pts <- t(inner_pts[!on_border,])
  
  #valid border pts need to be touching an interior point in the full 6-way direction
  touching <- apply(border_pts, 1, function(ipi){
    dists <- abs(finterior_border_pts - unlist(ipi))
    any(dists[1,]^2 + dists[2,]^2 < 2.1)
  })
  border_pts <- border_pts[touching,]
  
  #now to sort -- iterate through distances and pick next closest, bouncing around
  dbp <- as.matrix(dist(border_pts, method = "manhattan"))
  
  #now bounce along the rows of dbp, following a path of adjacent border_pts
  
  #initialize
  border_order <- rep(NA, nrow(border_pts))
  i <- 1
  border_order[i] <- 1
  
  #pick a direction
  i <- i + 1
  border_order[i] <- which(dbp[border_order[i-1],] == 1)[1]
  
  for(i in 3:length(border_order)){
    border_order[i] <- setdiff(which(dbp[border_order[i-1],] == 1), border_order[i-2])
  }
  
  border_pts <- border_pts[border_order,]
  
  if(plot_result){
    plot(rbind(shape, shape[1,]), type = "l", col = 2, lwd = 2, xlim = range(c(shape[,1], inter_pts[,1])),
         ylim = range(c(shape[,2], inter_pts[,2])))
    abline(v=floor(par("usr")[1]):ceiling(par("usr")[2]), lty = 3, col = adjustcolor(1, 0.5), lwd = 0.75)
    abline(h=floor(par("usr")[3]):ceiling(par("usr")[4]), lty = 3, col = adjustcolor(1, 0.5), lwd = 0.75)
    polygon(border_pts, col = adjustcolor(3, 0.25))
  }
  
  return(border_pts)
  
}

#ask the shape for directions
find_interior_shape_3 <- function(shape){
  
}

mat2img <- function(x){
  image_read(1-abind::abind(x, x, x, along = 3))
}

img2mat <- function(img){
  mat <- image_data(img)
  mat <- Reduce("+", lapply(1:3, function(i) matrix(as.integer(mat[i,,]), nrow(mat[i,,]), ncol(mat[i,,]))))
  mat <- t(scale2(abs(mat - median(mat))) >= 0.5)
  mat  
}

discretize_image <- function(img, ncol_range = 3:10, offset = 0){
  
  #simplify colorspace
  img <- round(img, 2)
  
  #retrieve and transform colors
  cols <- do.call(rbind, lapply(1:nrow(img), function(i) do.call(rbind, lapply(1:ncol(img), function(j) img[i,j,1:4]))))
  # transcols <- cbind(cols[,1:3] * cols[,4], 
  #                    rep(1:nrow(img), each = ncol(img)),
  #                    rep(1:ncol(img), times = nrow(img)))
  transcols <- cols
  transcols[transcols[,4] < 0.1,1:3] <- -1
  
  #perform clustering
  kmeans_out <- parallel::mclapply(ncol_range, function(k) 
    kmeans(transcols, k, iter.max = 50, nstart = 50), mc.cores = 12
  )
  
  #decide on k
  dbi <- sapply(kmeans_out, function(x) clusterSim::index.DB(transcols, x$cluster)$DB)
  kmeans_picked <- which.min(dbi) + offset
  kmeans_picked <- min(kmeans_picked, length(ncol_range))
  
  #retrieve clustering output
  cluster_inds <- kmeans_out[[kmeans_picked]]$cluster
  centers <- round(clusterSim::index.DB(cols, cluster_inds)$centers, 3)
  mean_cols <- apply(centers, 1, function(x) rgb(x[1], x[2], x[3], x[4]))
  
  return(list(cluster_inds = cluster_inds, mean_cols = mean_cols, centers = centers))

}

sqdists_to_centers <- function(mat, cents){
  if(is.null(nrow(cents))){
    cents <- matrix(cents, 1, length(cents))
  }
  sqdists <- do.call(cbind, lapply(1:nrow(cents), function(ki){
    colSums((t(mat) - cents[ki,])^2)
  }))
}

wkmeans <- function(x, w, k, kmeans.iter.max = NA, n = NA, eps = 1E-6){
  
  #hyperparameters
  if(is.na(n)){
    n <- k * 10
  }
  
  if(is.na(kmeans.iter.max)){
    kmeans.iter.max <- k * 5  
  }
  
  kmeans_out <- parallel::mclapply(1:n, function(kmeans_run){
    
    #initialize centers w/ kmeans++
    #TODO need to weigh probs by count!
    kmeans_iter <- 1
    converged <- F
    cents <- matrix(NA, k, 4)
    cents[1,] <- x[sample(1:nrow(x), 1),]
    for(ki in 2:k){
      sqdists <- sqdists_to_centers(x, cents[1:(ki-1),])
      cls <- apply(sqdists, 1, which.min)
      smallest_sqdists <- sqdists[cbind(1:length(cls), cls)]
      cents[ki,] <- x[sample(1:length(cls), 1, prob = smallest_sqdists),]
    }
    
    #evaluate initial cluster assignments
    sqdists <- sqdists_to_centers(x, cents)
    cls <- apply(sqdists, 1, which.min)
    
    #run kmeans iteratively
    while(!converged & kmeans_iter < kmeans.iter.max){
      
      #increment iter
      kmeans_iter <- kmeans_iter + 1
      
      #check that all clusters present and reassign cents if not
      prev_cents <- cents
      if(!all(1:k %in% cls)){
        empty_cls <- (1:k)[!(1:k %in% cls)]
        for(ki in empty_cls){
          sqdists <- sqdists_to_centers(x, cents)
          cls <- apply(sqdists, 1, which.min)
          smallest_sqdists <- sqdists[cbind(1:length(cls), cls)]
          cents[ki,] <- x[sample(1:length(cls), 1, prob = smallest_sqdists),]
        }
        sqdists <- sqdists_to_centers(x, cents)
        cls <- apply(sqdists, 1, which.min)
      }
      
      #find new centers
      cents <- do.call(rbind, lapply(1:k, function(ki){
        clusmems <- which(cls == ki) #members of the cluster
        clusvals <- x[clusmems,] 
        cluscounts <- w[clusmems]
        apply(t(clusvals) * rep(cluscounts, each = 4), 1, sum) / sum(cluscounts)
      }))
      
      #check for convergence
      converged <- all(abs(cents - prev_cents) < eps)
      
      #evaluate new cluster assignments
      if(!converged){
        sqdists <- sqdists_to_centers(x, cents)
        cls <- apply(sqdists, 1, which.min)
      }
      
    }
    
    #calculate total variance / sum of squares
    total_ss <- sum(sqdists[cbind(1:length(cls), cls)] * w)
    
    #return value
    list(converged = converged, k = k, kmeans_iter = kmeans_iter, cents = cents, cls = cls, total_ss = total_ss, obs_counts = w)
    
  }, mc.cores = 15)
  
  #pick best run
  total_sss <- sapply(kmeans_out, function(kmeans_indiv) kmeans_indiv$total_ss)
  optimal_cluster_index <- which.min(total_sss)
  optimal_cluster <- kmeans_out[[optimal_cluster_index]]
  return(optimal_cluster)
  
}

#more efficient attempt
discretize_image_2 <- function(img, ncol_range = 3:10, slow_clustering = F, n_top = 300, offset = 0){
  
  #simplify colorspace
  img <- round(img, 2)
  
  #retrieve and transform colors
  cols <- apply(img, 3, function(slice) t(slice))
  # transcols <- cbind(cols[,1:3] * cols[,4], 
  #                    rep(1:nrow(img), each = ncol(img)),
  #                    rep(1:ncol(img), times = nrow(img)))
  transcols <- cols
  transcols[transcols[,4] < 0.1,1:3] <- -1
  
  #rescale obs for more efficient kmeans
  colstrings <- apply(transcols, 1, paste0, collapse = "_")
  tcolstrings <- table(colstrings)
  fcolstrings <- do.call(rbind, lapply(strsplit(names(tcolstrings), "_"), as.numeric))
  
  #take n most popular colors and reassign the rest to those
  n_top <- min(n_top, length(tcolstrings))
  
  # plot(cumsum(sort(tcolstrings, T)) / sum(tcolstrings), type = "l")
  top_cols <- sort(tcolstrings, T)[1:n_top]
  prop_cols_represented <- sum(top_cols) / sum(tcolstrings)
  top_cols_mat <- do.call(rbind, lapply(strsplit(names(top_cols), "_"), as.numeric))
  
  #perform weighted clustering slowly (may need to rewrite own fast version)
  clust_out <- lapply(ncol_range, function(k){
    
    optimal_cluster <- wkmeans(x = top_cols_mat, w = top_cols, k)
    
    #reproduce original data
    sqdists <- sqdists_to_centers(fcolstrings, optimal_cluster$cents)
    cls <- apply(sqdists, 1, which.min)
    exploded_cls <- cls[match(colstrings, names(tcolstrings))]
    
    #should be equivalent
    exploded_cents <- do.call(rbind, lapply(1:k, function(ki){
      clusmems <- which(exploded_cls == ki) #members of the cluster
      clusvals <- transcols[clusmems,] 
      cluscounts <- rep(1, length(clusmems))
      apply(t(clusvals) * rep(cluscounts, each = 4), 1, sum) / sum(cluscounts)
    }))
    exploded_cents <- round(exploded_cents, 3)
    
    # yep, and faster too, weirdly
    # microbenchmark::microbenchmark(
    #   clusterSim::index.DB(transcols, exploded_cls)$centers,
    #   do.call(rbind, lapply(1:k, function(ki){
    #     clusmems <- which(exploded_cls == ki) #members of the cluster
    #     clusvals <- transcols[clusmems,] 
    #     cluscounts <- rep(1, length(clusmems))
    #     apply(t(clusvals) * rep(cluscounts, each = 4), 1, sum) / sum(cluscounts)
    #   }))
    # )
    
    #remove empty clusters so index.DB doesn't complain
    pres_cl <- sort(unique(exploded_cls))
    dbi_exploded_cl <- (1:length(pres_cl))[match(exploded_cls, pres_cl)]
    
    davies_bouldin_index <- clusterSim::index.DB(transcols, dbi_exploded_cl)$DB
    
    #return values
    return(list(cluster_inds = exploded_cls, centers = exploded_cents, dbi = davies_bouldin_index, k = k))
    
  })
  
  #decide on k
  dbi <- sapply(clust_out, function(x) x$dbi)
  k_picked <- which.min(dbi) + offset
  k_picked <- min(k_picked, length(ncol_range))
  
  #retrieve clustering output
  cluster_inds <- clust_out[[k_picked]]$cluster
  
  centers <- clust_out[[k_picked]]$centers
  centers[centers < 0] <- 0
  mean_cols <- apply(centers, 1, function(x) rgb(x[1], x[2], x[3], x[4]))
  
  return(list(cluster_inds = cluster_inds, mean_cols = mean_cols, centers = centers, prop_cols_represented = prop_cols_represented))
  
}

#get rid of noisy pixels easily
simplify_image <- function(cluster_inds, img, min_vl = 3, n_passes = 2){
  clmat <- matrix(cluster_inds, nrow(img), ncol(img), byrow = T)
  for(pass in 1:(2*n_passes)){
    clmat <- apply(clmat, 1, function(x){
      rlex <- rle(x)
      ok <- which(rlex$lengths >= min_vl)
      not_ok <- setdiff(1:length(rlex$lengths), ok)
      if(length(not_ok) == 0){
        return(x)
      } else {
        where_to_add <- sapply(not_ok, function(i) min(sum(i > ok) + 1, length(ok)))
        ok_vals <- rlex$values[ok]
        ok_lengths <- rlex$lengths[ok]
        not_ok_lengths <- rlex$lengths[not_ok]
        uniq_w2a <- sort(unique(where_to_add))
        n2a <- sapply(uniq_w2a, function(i) sum(not_ok_lengths[where_to_add == i]))
        ok_lengths[uniq_w2a] <- ok_lengths[uniq_w2a] + n2a
        rep(ok_vals, times = ok_lengths)
      }
    })
  }
  return(clmat)
}

n_nbors <- function(locs, mat){
  sum(mat[t(locs + cbind(c(0,1), c(1,0), c(0,-1), c(-1,0)))])
}

deconvolve_1kern <- function(boolmat, dimk){
  boolmat_inds <- which(boolmat, T)
  dimkr <- c(ceiling(-dimk/2), floor(dimk/2)) - c(0, dimk %% 2 == 0) #our window around each convolved point
  dimkr_combos <- expand.grid(dimkr[1]:dimkr[2], dimkr[1]:dimkr[2]) #vals to add
  truevals <- do.call(rbind, lapply(1:nrow(dimkr_combos), function(dkri) 
    t(t(boolmat_inds) + unlist(dimkr_combos[dkri,])))) #do the addition
  
  #filter out impossible or duplicated values
  xr <- c(1, ncol(boolmat))
  yr <- c(1, nrow(boolmat))
  truevals <- truevals[truevals[,1] >= yr[1] &
                         truevals[,1] <= yr[2] &
                         truevals[,2] >= xr[1] &
                         truevals[,2] <= xr[2], ]
  #this is costlier than just doing the redundant assignment
  # truevals <- truevals[!duplicated(truevals),]
  
  #create and return new matrix of bools
  xdc <- boolmat
  xdc[truevals] <- T
  xdc
}

separate_layers <- function(onehot, clmat, use_fourier = F, use_hadamard = F, use_wavelet = F, kconv_init = 1){
  
  #useful constants
  if(all(!c(use_wavelet, use_fourier, use_hadamard))){use_wavelet <- T}
  kern <- matrix(c(1,2,1,0,0,0,-1,-2,1), 3, 3)
  nesw <- cbind(c(0,1), c(1,0), c(0,-1), c(-1,0))
  dirs <- t(expand.grid(-1:1, -1:1))
  dirs <- dirs[,!apply(dirs == c(0,0), 2, all)]
  
  #initialize variables
  new_shapes <- old_shapes <- list()
  
  #iterate through layers
  for(layer in 1:length(onehot)){
    
    #get cluster matrix
    cmat <- onehot[[layer]]
    
    #check to make sure pixels exist (didn't get eliminated earlier)
    if(sum(cmat) == 0) {
      old_shapes[[layer]] <- list(integer(0))
      next()
    }
    
    #PROBLEM -- since inner and outer shapes are affected differently by this convolution,
    #they no longer overlap
    #also, you lose information on fine shapes
    #ANOTHER IDEA -- maybe make all 2x2 blocks 3x3? no problems with overlaps then
    
    #now start processing the image
    img_temp <- mat2img(cmat)
    
    #make sure all pixels are at least 4x4 blobs (or else edge detection fails)
    img_temp <- image_convolve(img_temp, kernel = matrix(1,kconv_init,kconv_init), bias = '50%')
    
    #deconvolve the resultant image
    img_temp <- deconvolve_1kern(boolmat = img2mat(img_temp), dimk = kconv_init + 1)
    
    #get the two transposes for edge detection
    imgm1 <- mat2img(img_temp)
    imgm2 <- mat2img(t(img_temp))
    
    #find edges
    # skern <- matrix(c(1,2,1,0,0,0,-1,-2,-1), 3, 3)
    skern <- "sobel"
    out1 <- image_data(image_convolve(imgm1, skern, bias = '50%'))
    out2 <- image_data(image_convolve(imgm2, skern, bias = '50%'))
    out1 <- Reduce("+", lapply(1:3, function(i) matrix(as.integer(out1[i,,]), nrow(out1[i,,]), ncol(out1[i,,]))))
    out2 <- Reduce("+", lapply(1:3, function(i) matrix(as.integer(out2[i,,]), nrow(out2[i,,]), ncol(out2[i,,]))))
    # plot(which(scale2(abs(out1 - median(out1))) > 0.5, T))
    # plot(which(scale2(abs(t(out2) - median(t(out2)))) > 0.5, T))
    
    #find edges
    edges <- scale2(abs(out1 - median(out1))) > 0.5 | scale2(abs(t(out2) - median(t(out2)))) > 0.5
     
    #now pass over these with a 2x2 convolution of all 1s to simplify edge data
    out_edges <- image_convolve(mat2img(edges), kernel = matrix(1,2,2), bias = '50%')
    edges <- img2mat(out_edges)
    
    #top-right edge gets shifted, so gotta correct by shrinking
    # inds <- which(edges, T)
    # new_inds <- lapply(1:nrow(inds), function(ei){
    #   
    #   edge <- inds[ei,c(2,1)]
    #   nbors <- t(nesw + edge)
    #   nbors <- nbors[nbors[,1] >= 1 &
    #                    nbors[,1] <= nrow(cmat) &
    #                    nbors[,2] >= 1 &
    #                    nbors[,2] <= ncol(cmat),]
    #   
    #   # xr <- edge[1] + c(-20,20)
    #   # yr <- edge[2] + c(-20,20)
    #   # plot(inds[,2], inds[,1], xlim = xr, ylim = yr, cex = 1.1)
    #   # abline(v = xr[1]:xr[2], lwd = 0.2)
    #   # abline(h = yr[1]:yr[2], lwd = 0.2)
    #   # points(which(t(img_temp), T) %*% matrix(c(0,1,1,0), 2), pch = 19, col = adjustcolor(4, 0.2))
    #   # points(edge[1], edge[2], col = 2, pch = 19)
    #   # # text(inds[,2], inds[,1], labels = 1:nrow(inds), col = 2, cex = 0.5)
    #   # points(nbors, col = 2)
    #   # text(nbors[,1], nbors[,2], labels = 1:4, col = 2, cex = 0.5)
    #   
    #   interior_nbor <- which(img_temp[nbors])
    #   
    #   extended_nbors <- t(dirs + edge)
    #   extended_nbors <- extended_nbors[extended_nbors[,1] >= 1 &
    #                                      extended_nbors[,1] <= nrow(cmat) &
    #                                      extended_nbors[,2] >= 1 &
    #                                      extended_nbors[,2] <= ncol(cmat),]
    #   interior_enbor <- which(img_temp[extended_nbors])
    #   
    #   
    #   if(img_temp[edge[1], edge[2]]){
    #     return(edge)
    #   } else { #ok, this is hacky, but works for now
    #     if(length(interior_nbor) == 0){
    #       if(img_temp[edge[1] - 1, edge[2] - 1]){
    #         return(edge + c(-1,-1))
    #       } else {
    #         return(edge + c(-1,-1))
    #       } 
    #     }
    #     
    #     if(length(interior_nbor) == 1){return(nbors[interior_nbor,])}
    #     
    #     #means we're in a corner, and also need to account for the entire corner
    #     if(length(interior_nbor) == 2){
    #       if(all(sort(interior_nbor) == c(1,2))){return(rbind(edge + c(1,1), nbors[interior_nbor,]))}
    #       if(all(sort(interior_nbor) == c(2,3))){return(rbind(edge + c(1,-1), nbors[interior_nbor,]))}
    #       if(all(sort(interior_nbor) == c(3,4))){return(rbind(edge + c(-1,-1), nbors[interior_nbor,]))}
    #       if(all(sort(interior_nbor) == c(1,4))){return(rbind(edge + c(-1,1), nbors[interior_nbor,]))}
    #     }
    #     
    #     if(length(interior_nbor) > 2){
    #       return(extended_nbors[interior_enbor,])
    #     }
    #   }
    # })
    # 
    # new_inds <- do.call(rbind, new_inds)[,2:1]
    # new_inds <- new_inds
    # new_inds <- new_inds[!duplicated(new_inds),]
    # inds <- new_inds
    # plot(inds)
    # plot(inds[,2], inds[,1], xlim = c(200,260), ylim = c(300, 350))
    # points(new_inds[,2], new_inds[,1], col = 2, cex = 1.2)
    # points(which(t(img_temp), T) %*% matrix(c(0,1,1,0), 2), pch = 19, col = adjustcolor(4, 0.4))
    # #take only first entry of every listed new_inds for matching to work
    # #and don't remove duplicated entries
    # segments(x0 = inds[,2], y0 = inds[,1], x1 = new_inds[,2], y1 = new_inds[,1])
    #actually, this is a much faster way to shrink shapes lol
    #and then just do old uni-directional chain
    
    
    #hmm lemme try rewriting this...
    #we want every ind on the outer border's inside
    #which is to say all the current inds that are in the space
    #but only those with adjacent values outside the space
    #so can first expand each existing ind by dir
    #and then filter out these new inds by adjacency to the outside
    #maybe not super efficient, but still linear in ind
    inds <- which(edges, T)
    new_inds <- lapply(1:nrow(inds), function(ei){
      
      edge <- inds[ei,c(2,1)]
      nbors <- t(dirs + edge)
      nbors <- nbors[nbors[,1] >= 1 &
                       nbors[,1] <= nrow(cmat) &
                       nbors[,2] >= 1 &
                       nbors[,2] <= ncol(cmat),]
      interior_nbors <- matrix(nbors[which(img_temp[nbors]),], ncol = 2)
      
      borders_outside <- sapply(1:nrow(interior_nbors), function(ini){
        nbors_nbors <- t(dirs + interior_nbors[ini,] )
        any(!img_temp[nbors_nbors])
      })
      
      return(interior_nbors[borders_outside,])
      
    })
    
    new_inds <- do.call(rbind, new_inds)[,2:1]
    new_inds <- new_inds[!duplicated(new_inds),]
    inds <- new_inds
    # plot(inds)
    # points(which(t(img_temp), T), pch = 19, col = adjustcolor(4, 0.2))
    
    #recreate edge matrix
    edges <- matrix(F, nrow(edges), ncol(edges))
    edges[inds] <- T
    
    #potentially still some 2x2 blocks left? -- check and correct if so
    #earlier processing may have fixed this though
    out_edges_2 <- image_convolve(mat2img(edges), kernel = matrix(1,2,2), bias = '50%')
    edges_extra <- img2mat(out_edges_2)
    max_2x2_conv <- 10
    convi <- 1
    while(any(edges_extra) & convi <= max_2x2_conv){
      convi <- convi + 1
      # print(convi)
      
      #erase offending values
      edges[which(edges_extra,T)] <- F #hmm, this leaves nubbins & gaps
      
      #find any resultant gaps & nubbins
      inds <- which(edges, T)
      # plot(inds, xlim = c(320,350), ylim = c(200, 250))
      # points(which(t(cmat), T), pch = 19, col = adjustcolor(1, 0.1)) #top right edge is shifted up and right
      # points(inds)
      # points(which(edges_extra, T), col = 2, pch = 19)
      
      ind_neighbors <- apply(inds, 1, n_nbors, edges)
      gaps_and_nubbins <- matrix(inds[ind_neighbors == 1,], ncol = 2)
      # points(gaps_and_nubbins, col = 3, pch = 19)
      # text(x = gaps_and_nubbins[,1], y = gaps_and_nubbins[,2], labels = 1:nrow(gaps_and_nubbins), col = 1)
      
      #gaps and nubbins whose neighbor has 2 connections are gaps (not 1, imposs w/ earlier filtering)
      #but if neighbor has 3 connections they're nubbins
      if(nrow(gaps_and_nubbins) != 0){
        nbors_nnbors <- apply(gaps_and_nubbins, 1, function(locs){
          nbor <- t(locs + nesw)
          nbor <- nbor[edges[nbor],]
          n_nbors(nbor, edges)
        })
        
        #delete nubbins
        edges[gaps_and_nubbins[nbors_nnbors == 3,]] <- F
        
        #identify gaps
        gaps <- gaps_and_nubbins[nbors_nnbors == 2,]
        
        #gaps always have another gap on a diagonal (aka manh dist of 2), so they pair off 
        gap_dist_pairs <- which(as.matrix(dist(gaps, method = "man")) == 2, T)
        gap_dist_pairs <- t(apply(gap_dist_pairs, 1, sort))
        gap_dist_pairs <- matrix(gap_dist_pairs[!duplicated(gap_dist_pairs),], ncol = 2)
        
        #find gaps to fill
        if(nrow(gap_dist_pairs) != 0){
          gap_empty_locs <- t(apply(gap_dist_pairs, 1, function(gis){
            g1 <- gaps[gis[1],]
            g1n <- t(g1 + nesw)
            
            # points(g1[1], g1[2], col = 2, pch = 19)
            g2 <- gaps[gis[2],]
            g2n <- t(g2 + nesw)
            # points(g2[1], g2[2], col = 2, pch = 19)
            
            #find a shared neighbor
            gn <- rbind(g1n, g2n)
            shared_pts <- matrix(gn[duplicated(gn),], ncol = 2)
            shared_pts[1,]
          }))
          
          #fill gaps
          edges[gap_empty_locs] <- T
        }
      }
      
      # # find 3-way junctions... maybe fixed this earlier!
      # trois <- inds[ind_neighbors == 3,]
      # 
      # #combine these with gaps to find shapes
      # points(trois, pch = 19, col = 2)
      # points(gaps, pch = 19, col = 3)
      # gaps
      
      #reprocess
      out_edges_2 <- image_convolve(mat2img(edges), kernel = matrix(1,2,2), bias = '50%')
      edges_extra <- img2mat(out_edges_2)
    }
    
    #retrieve indices
    inds <- which(edges, T)
    # plot(inds, xlim = c(1, ncol(cmat)), ylim = c(1, nrow(cmat)))
    # plot(inds, xlim = c(320,350), ylim = c(200, 250))
    
    #do one last check for gaps and nubbins
    ind_neighbors <- apply(inds, 1, n_nbors, edges)
    gaps_and_nubbins <- matrix(inds[ind_neighbors == 1,], ncol = 2)
    
    #gaps and nubbins whose neighbor has 2 connections are gaps (not 1, imposs w/ earlier filtering)
    #but if neighbor has 3 connections they're nubbins
    if(nrow(gaps_and_nubbins) != 0){
      nbors_nnbors <- apply(gaps_and_nubbins, 1, function(locs){
        nbor <- t(locs + nesw)
        nbor <- nbor[edges[nbor],]
        n_nbors(nbor, edges)
      })
      
      #delete nubbins
      edges[gaps_and_nubbins[nbors_nnbors == 3,]] <- F
      
      #identify gaps
      gaps <- gaps_and_nubbins[nbors_nnbors == 2,]
      
      #gaps always have another gap on a diagonal (aka manh dist of 2), so they pair off 
      gap_dist_pairs <- which(as.matrix(dist(gaps, method = "man")) == 2, T)
      gap_dist_pairs <- t(apply(gap_dist_pairs, 1, sort))
      gap_dist_pairs <- matrix(gap_dist_pairs[!duplicated(gap_dist_pairs),], ncol = 2)
      
      #find gaps to fill
      if(nrow(gap_dist_pairs) != 0){
        gap_empty_locs <- t(apply(gap_dist_pairs, 1, function(gis){
          g1 <- gaps[gis[1],]
          g1n <- t(g1 + nesw)
          
          # points(g1[1], g1[2], col = 2, pch = 19)
          g2 <- gaps[gis[2],]
          g2n <- t(g2 + nesw)
          # points(g2[1], g2[2], col = 2, pch = 19)
          
          #find a shared neighbor
          gn <- rbind(g1n, g2n)
          shared_pts <- matrix(gn[duplicated(gn),], ncol = 2)
          shared_pts[1,]
        }))
        
        #fill gaps
        edges[gap_empty_locs] <- T
      }
      
      #any remaining values are nubbins and need to be deleted
      ind_neighbors <- apply(inds, 1, n_nbors, edges)
      nnc_i <- which(ind_neighbors != 2)
      
      # while(length(gnc_i) != 0){
        
        nubbins_and_cores <- matrix(inds[nnc_i,], ncol = 2)
        nnc_pairs <- matrix(which(as.matrix(dist(nubbins_and_cores, "manh")) == 1, T), ncol = 2)
        
        if(nrow(nnc_pairs) != 0){
          nnc_pairs <- t(apply(nnc_pairs, 1, sort))
          nnc_pairs <- matrix(nnc_pairs[!duplicated(nnc_pairs),], ncol = 2)
          if(nrow(nnc_pairs) != 0){
            for(nncp_i in 1:nrow(nnc_pairs)){
              nnc_pair <- nnc_i[nnc_pairs[nncp_i,]]
              nnc_to_del <- inds[nnc_pair[which(ind_neighbors[nnc_pair] == 1)],]
              edges[nnc_to_del[1], nnc_to_del[2]] <- F
            }
          }
        } else if(nrow(nnc_pairs) / 2 != length(nnc_i) & length(nnc_i) > 0){
          accounted_for <- unique(c(nnc_pairs))
          if(length(accounted_for) != 0){
            leftovers <- nubbins_and_cores[-accounted_for,]  
          } else {
            leftovers <- nubbins_and_cores
          }
          edges[leftovers] <- F
        }
        
        inds <- which(edges, T)
        ind_neighbors <- apply(inds, 1, n_nbors, edges)
        nnc_i <- which(ind_neighbors != 2)
        
      # }
      
      
      
      
    }
    
    
    #now let's hack together a chain to go around the edges
    edges2process <- edges
    
    shapes <- list()
    while(any(edges2process)){
      inds2process <- which(edges2process, T)
      foc_ind <- inds2process[1,]
      curr_ind_in_chain <- 1
      running_inds <- list(x = numeric(0), y = numeric(0))
      
      has_adjacent_vals <- T
      while(has_adjacent_vals){
        edges2process[foc_ind[1], foc_ind[2]] <- F
        running_inds$x[curr_ind_in_chain] <- foc_ind[1]
        running_inds$y[curr_ind_in_chain] <- foc_ind[2]
        nesw_inds <- t(foc_ind + nesw)
        adj_vals <- edges2process[nesw_inds]
        has_adjacent_vals <- any(adj_vals)
        
        if(has_adjacent_vals){
          connected_ind_ind <- which(adj_vals)[1]
          connected_ind <- nesw_inds[connected_ind_ind,]
          foc_ind <- connected_ind
          curr_ind_in_chain <- curr_ind_in_chain + 1
        }
      }
      
      #do quick checks
      tentative_shape <- data.frame(x = running_inds$x, y = running_inds$y)
      returned_to_start <- sum(abs(head(tentative_shape, 1) - tail(tentative_shape, 1))) < 1.01
      long_enough <- nrow(tentative_shape) > 8.9
      
      if(returned_to_start & long_enough){
        shapes[[length(shapes)+1]] <- tentative_shape
      } else {
        if(!returned_to_start){print("failed returned_to_start")}
        if(!long_enough){print("failed long_enough")}
      }
      
    }
    
    old_shapes[[layer]] <- shapes
    
  }
  
  #now we can smooth this... can either 
  #1) fit a parametric shape to it, 
  #2) apply non-parametric smoothing, eg laplace smoothing or a spline
  #3) compute the fourier / hadamard / wavelet transform and discard extra information
  #4) use a statistical method or else smoothly blend multiple transforms in proportion to their fit goodness
  for(layer in 1:length(onehot)){ 
    
    cmat <- onehot[[layer]]
    if(sum(cmat) == 0) {
      new_shapes[[layer]] <- list(integer(0))
      next()
    }
    
    shapes <- old_shapes[[layer]]
    
    if(use_fourier){
      
      fshapes <- lapply(shapes, function(running_inds){
        
        #weird edge case
        if(nrow(running_inds) == 1){return(running_inds)}
        
        #transform and extract magnitudes
        fftout <- fft(as.matrix(running_inds), F)
        mag <- Mod(fftout)
        
        #magnitude-based filtering
        noise_sigma <- median(mag[(nrow(fftout) %/% 2 + 1):nrow(fftout), (ncol(fftout) %/% 2 + 1):ncol(fftout)], na.rm = TRUE) / 0.6745
        threshmag <- noise_sigma * sqrt(2 * log(nrow(running_inds))) #universal threshold
        threshmag <- noise_sigma * sqrt(2 * log(nrow(running_inds))^2 / log(2)) #visushrink threshold
        to_discard <- lapply(1:2, function(i){
          x <- mag[,i]
          which(x < threshmag)
        })
        to_discard <- intersect(to_discard[[1]], to_discard[[2]])
        
        #recompose
        fftout[to_discard,] <- 0 + 0i
        filtered_inds <- Re(fft(fftout, T)) / nrow(running_inds) / ncol(running_inds)
        data.frame(filtered_inds)
      })
      
    } else if(use_hadamard){
      
      fshapes <- lapply(shapes, function(running_inds){
        
        #weird edge case
        if(nrow(running_inds) == 1){return(running_inds)}
        
        #transform and extract magnitudes
        hout <- gsignal::fwht(as.matrix(running_inds))
        mag <- Mod(hout)
        
        # Use a subset of coefficients to estimate the noise standard deviation
        n_coef_sd <- 20
        noise_sigma_1 <- sd(hout[,1][(length(hout[,1]) - n_coef_sd + 1):length(hout[,1])])
        noise_sigma_2 <- sd(hout[,2][(length(hout[,2]) - n_coef_sd + 1):length(hout[,2])])
        
        # Set the threshold level using a multiplier and the noise standard deviation
        threshold_multiplier <- 1
        threshold_1 <- threshold_multiplier * noise_sigma_1
        threshold_2 <- threshold_multiplier * noise_sigma_2
        
        # Apply soft thresholding to the Hadamard coefficients
        thresh_coeffs_1 <- sign(hout[,1]) * (abs(hout[,1]) - threshold_1) * (abs(hout[,1]) > threshold_1)
        thresh_coeffs_2 <- sign(hout[,2]) * (abs(hout[,2]) - threshold_2) * (abs(hout[,2]) > threshold_2)
        
        denoised_signal_1 <- gsignal::ifwht(thresh_coeffs_1)
        denoised_signal_2 <- gsignal::ifwht(thresh_coeffs_2)
        
        # Combine the denoised signals into a 2D time series
        filtered_inds <- cbind(denoised_signal_1, denoised_signal_2)
        
        #magnitude-based filtering
        # noise_sigma <- median(mag[(nrow(hout) %/% 2 + 1):nrow(hout), (ncol(hout) %/% 2 + 1):ncol(hout)], na.rm = TRUE) / 0.6745
        # threshmag <- noise_sigma * sqrt(2 * log(nrow(running_inds))) #universal threshold
        # threshmag <- noise_sigma * sqrt(2 * log(nrow(running_inds))^2 / log(2)) #visushrink threshold
        # to_discard <- lapply(1:2, function(i){
        #   x <- mag[,i]
        #   which(x < mean(x) * threshmag)
        # })
        # to_discard <- intersect(to_discard[[1]], to_discard[[2]])
        # 
        # #recompose
        # hout[to_discard,] <- 0
        # filtered_inds <- Re(gsignal::ifwht(hout)) / nrow(running_inds) / ncol(running_inds)
        filtered_inds <- data.frame(filtered_inds[1:nrow(running_inds),])
        filtered_inds
      })
      
    } else if(use_wavelet){
      
      fshapes <- lapply(shapes, function(running_inds){
        
        #weird edge case
        if(nrow(running_inds) == 1){return(running_inds)}
        
        #transform and extract magnitudes
        wavelet_filter <- c("haar", "la8", "d4", "d6", "mb4")[5]
        n_levels <- max(1, floor(log2(nrow(running_inds)) / 2))
        
        waveout.1 <- waveslim::modwt(running_inds[,1], wf = wavelet_filter, n.levels = n_levels)
        waveout.2 <- waveslim::modwt(running_inds[,2], wf = wavelet_filter, n.levels = n_levels)
        
        # Get the detail coefficients from the first level
        detail_coeffs_1 <- waveout.1$d1
        detail_coeffs_2 <- waveout.2$d1
        
        # Estimate noise_sigma using the Median Absolute Deviation (MAD)
        noise_sigma_1 <- sd(waveout.1$d1)
        noise_sigma_2 <- sd(waveout.2$d1)
        
        # Set the threshold level using a multiplier and the noise standard deviation
        threshold_multiplier <- 3
        threshold_1 <- threshold_multiplier * noise_sigma_1
        threshold_2 <- threshold_multiplier * noise_sigma_2
        
        # Apply soft thresholding to the detail coefficients
        thresholded_coeffs_1 <- waveout.1
        thresholded_coeffs_2 <- waveout.2
        
        for (i in 1:n_levels) {
          level_name <- paste0("d", i)
          thresholded_coeffs_1[[level_name]] <- sign(thresholded_coeffs_1[[level_name]]) * 
            (abs(thresholded_coeffs_1[[level_name]]) - threshold_1) * (abs(thresholded_coeffs_1[[level_name]]) > threshold_1)
          thresholded_coeffs_2[[level_name]] <- sign(thresholded_coeffs_2[[level_name]]) * 
            (abs(thresholded_coeffs_2[[level_name]]) - threshold_2) * (abs(thresholded_coeffs_2[[level_name]]) > threshold_2)
        }
        
        denoised_signal_1 <- waveslim::imodwt(thresholded_coeffs_1)
        denoised_signal_2 <- waveslim::imodwt(thresholded_coeffs_2)
        
        # Combine the denoised signals into a 2D time series
        filtered_inds <- cbind(denoised_signal_1, denoised_signal_2)
        data.frame(filtered_inds)
      })
      
    }
    
    new_shapes[[layer]] <- fshapes
    
  }
  
  #find outermost box and add to shapes
  bounding_box <- apply(do.call(rbind, lapply(old_shapes, do.call, what = rbind)), 2, range) + matrix(c(-1,1,-1,1),2,2)
  bounding_box <- expand.grid(bounding_box[,1], bounding_box[,2])[c(1,2,4,3),]
  bounding_box <- do.call(rbind, lapply(1:4, function(i) interpolate(bounding_box[i,], bounding_box[wrap(i+1, 4),])))
  outer_box <- matrix(c(1, nrow(clmat), 1, ncol(clmat)), 2, 2)
  outer_box <- expand.grid(outer_box[,1], outer_box[,2])[c(1,2,4,3),]
  outer_box <- do.call(rbind, lapply(1:4, function(i) interpolate(outer_box[i,], outer_box[wrap(i+1, 4),])))
  outer_cluster <- as.numeric(names(which.max(table(clmat[outer_box]))))
  old_shapes[[outer_cluster]][[length(old_shapes[[outer_cluster]]) + 1]] <- outer_box[,c(2,1)]
  new_shapes[[outer_cluster]][[length(new_shapes[[outer_cluster]]) + 1]] <- outer_box[,c(2,1)]
  
  #re-orient so output is upright
  new_shapes <- lapply(new_shapes, function(layer) lapply(layer, function(shape){
    if(length(layer[[1]]) == 0){return(layer)}
    out <- data.frame(as.matrix(shape) %*% matrix(c(1,0,0,-1),2,2))
    colnames(out) <- c("x", "y")
    out
  }))
  
  #re-orient so output is upright
  old_shapes <- lapply(old_shapes, function(layer) lapply(layer, function(shape){
    if(length(layer[[1]]) == 0){return(layer)}
    out <- data.frame(as.matrix(shape) %*% matrix(c(1,0,0,-1),2,2))
    colnames(out) <- c("x", "y")
    out
  }))
  
  return(list(old_shapes = old_shapes, new_shapes = new_shapes))
}


find_shape_plotting_order <- function(new_shapes){
  
  #check if shape contains > 5% of volume with color -- if it does, plot it, if it doesn't, remove it
  #check if point on border is inside other polygon -- generate s x s matrix where s = # of shapes
  
  #Find matrix of contained within values (-1 if on the outside, 0 if neither, 1 if on the inside), 
  #sum the rows to get priority, take order priority to get order of painting
  
  #get useful variables
  new_shape_inds <- do.call(rbind, lapply(1:length(new_shapes), function(i) 
    cbind(layer = i, shape = 1:length(new_shapes[[i]]))))
  n_shapes <- nrow(new_shape_inds)
  
  #identify proportion of shape overlap
  prop_overlap <- lapply(2:n_shapes, function(i){
    lapply(1:(i-1), function(j){
      
      s1 <- new_shape_inds[i,]
      s2 <- new_shape_inds[j,]
      
      p1 <- sf::st_polygon(list(as.matrix(rbind(new_shapes[[s1[1]]][[s1[2]]], new_shapes[[s1[1]]][[s1[2]]][1,]))))
      p2 <- sf::st_polygon(list(as.matrix(rbind(new_shapes[[s2[1]]][[s2[2]]], new_shapes[[s2[1]]][[s2[2]]][1,]))))
      
      sf_p1 <- sf::st_sf(geometry = sf::st_sfc(p1))
      sf_p2 <- sf::st_sf(geometry = sf::st_sfc(p2))
      
      isec <- sf::st_intersection(sf_p1, sf_p2)
      
      a_p1 <- sf::st_area(sf_p1)
      a_p2 <- sf::st_area(sf_p2)
      a_isec <- sf::st_area(isec)
      
      if(length(a_isec) == 0){
        return(c(0,0))
      } else {
        return(c(a_isec / a_p1, a_isec / a_p2))
      }
      
    })
  })
  
  #construct matrix of proportion overlaps between shapes
  poverlap <- diag(n_shapes)
  for(i in 1:length(prop_overlap)){
    for(j in 1:length(prop_overlap[[i]])){
      poverlap[i+1,j] <- prop_overlap[[i]][[j]][1]
      poverlap[j,i+1] <- prop_overlap[[i]][[j]][2]
    }
  }
  
  #identify if shape contains or is contained by other shapes
  contains <- which(poverlap > 0.95 & t(poverlap) < 0.95, T)
  is_contained <- which(poverlap < 0.95 & t(poverlap) > 0.95, T)
  same_shapes <- which(poverlap > 0.95 & t(poverlap) > 0.95, T)
  same_level <- which(poverlap < 0.95 & t(poverlap) < 0.95, T)
  
  #construct containment matrix
  shape_hierarchy = list(contains = contains,
                         is_contained = is_contained,
                         same_shapes = same_shapes,
                         same_level = same_level)
  shape_hierarchy_vals = c(contains = -1,
                           is_contained = 1,
                           same_shapes = 0,
                           same_level = 0)
  containment_matrix <- diag(nrow(new_shape_inds))
  for(i in names(shape_hierarchy)){
    if(length(shape_hierarchy[[i]]) == 0){next()}
    for(j in 1:nrow(shape_hierarchy[[i]])){
      containment_matrix[shape_hierarchy[[i]][j,1], 
                         shape_hierarchy[[i]][j,2]] <- shape_hierarchy_vals[i]
      
    }
  }
  
  
  #identify impostor shapes (that define external, not internal colors)
  n_contains <- apply(containment_matrix, 1, function(x) sum(x == 1))
  impostor_pairs <- same_shapes[apply(same_shapes, 1, function(x) x[1] != x[2]),]
  impostor_pairs <- impostor_pairs[duplicated(t(apply(impostor_pairs, 1, sort))),]
  impostor_pairs <- impostor_pairs[order(n_contains[impostor_pairs[,1]], decreasing = T),] #resolve top disagreements first
  
  impostors <- rep(0, nrow(impostor_pairs))
  if(nrow(impostor_pairs) > 0){
    for(i in 1:nrow(impostor_pairs)){ #can also vectorize within rank
      pair <- impostor_pairs[i,]
      contained_by <- which(containment_matrix[pair[1],] == -1) #should always have at least 1
      next_shapes_up <- setdiff(contained_by[n_contains[contained_by] == min(n_contains[contained_by])], impostors)
      layers_nsu <- new_shape_inds[next_shapes_up,1]
      impostors[i] <- pair[new_shape_inds[pair,1] %in% layers_nsu]
    }  
  }
  
  #compile final plotting order
  plotting_ranks <- ceiling(rank(apply(-1 * containment_matrix, 1, sum)))
  plotting_order <- order(plotting_ranks)
  plotting_order <- setdiff(plotting_order, impostors)
  
  #return result
  return(list(plotting_order = plotting_order, new_shape_inds = new_shape_inds))
  
}



# centers[,4] <- as.numeric(centers[,4] > 0.95)

slow_method <- F
if(slow_method){
  #hmm this is too slow :/
  shapes_contain_color <- lapply(1:length(old_shapes), function(layer){
    print(layer)
    layer_cols <- onehot[[layer]]
    shapes <- old_shapes[[layer]]
    sapply(shapes, function(shape){
      
      #find elements just inside the shape
      interior_shape <- find_interior_shape(shape, T)
      second_interior_shape <- find_interior_shape(interior_shape, T)
      
      debug <- F
      if(debug){
        plot(shape, type = "l")
        lines(interior_shape, col = 2)
        lines(second_interior_shape, col = 3)
      }
      
      #test for color presence
      mean(layer_cols[as.matrix(interior_shape)]) > 0.2 |
        mean(layer_cols[as.matrix(second_interior_shape)]) > 0.2
      
    })
  })
  
  # dirs <- t(expand.grid(-1:1, -1:1))
  # dirs <- dirs[,!apply(dirs == c(0,0), 2, all)]
  # shapes_contain_color_fast <- lapply(1:length(old_shapes), function(layer){
  #   print(layer)
  #   layer_cols <- onehot[[layer]]
  #   shapes <- old_shapes[[layer]]
  #   sapply(shapes, function(shape){
  #     
  #     #find elements just inside the shape
  #     n_random_pts <- 20
  #     n_random_pts <- min(nrow(shape), n_random_pts)
  #     random_pts_on_shape <- shape[sample(1:nrow(shape), n_random_pts, replace = F),]
  #     
  #     random_pts_by_shape <- do.call(rbind, lapply(1:n_random_pts, function(i) t(unlist(random_pts_on_shape[i,]) + dirs)))
  #     random_pts_by_shape <- random_pts_by_shape[!duplicated(random_pts_by_shape),]
  #     random_pts_by_shape <- remove_points_no_dupes(random_pts_by_shape, random_pts_on_shape)
  #     
  #     
  #     random_pts_in_shape <- ?sp::point.in.polygon(point.x = random_pts_by_shape[,1], point.y = random_pts_by_shape[,2],
  #                                            pol.x = shape$x, pol.y = shape$y)
  #     random_pts_in_shape <- random_pts_by_shape[random_pts_in_shape == 1,]
  #     
  #     #test for color presence
  #     mean(layer_cols[as.matrix(random_pts_by_shape)]) > 0.4 |
  #       mean(layer_cols[as.matrix(random_pts_by_shape)]) > 0.4
  #     
  #   })
  # })
}


pad <- function(x, n=2, val=0, center=T){
  if(center){
    pl <- round(n / 2 - 1E-3)
    pr <- round(n / 2 + 1E-3) 
    pb <- round(n / 2 - 1E-3)
    pt <- round(n / 2 + 1E-3)  
  } else {
    pl <- 0
    pr <- n
    pb <- n
    pt <- 0
  }
  
  cbind(matrix(val, nrow(x) + pb + pt, pl), 
        rbind(matrix(val, pt, ncol(x)), x, matrix(val, pb, ncol(x))),
        matrix(val, nrow(x) + pb + pt, pr))
}

trim <- function(x, n=2, center=T){
  if(center){
    pl <- round(n / 2 - 1E-3)
    pr <- round(n / 2 + 1E-3) 
    pb <- round(n / 2 - 1E-3)
    pt <- round(n / 2 + 1E-3)  
  } else {
    pl <- 0
    pr <- n
    pb <- n
    pt <- 0
  }
  x[(pl+1):(ncol(x)-pr), (pt+1):(ncol(x)-pb)]
}

scale2 <- function(x) {(x - min(x)) / ifelse(diff(range(x)) != 0, diff(range(x)), Inf)}

#### process an example image ####

#retrieve image
img_name <- c("narutomaki_smol", "openai", "stephen_southparlk_head", "deathrats1", 
              "stan_logo", "rat_anatomy", "treadmill_rat", "smile", "stephen_southparlk_head_smol")[8]
path_to_image <- paste0("~/Pictures/", img_name, ".png")
img <- readPNG(source = path_to_image)
resmul <- 1
new_dims <- dim(img) * c(resmul, resmul, 1)
big_img <- array(rep(0, prod(new_dims)), dim = new_dims)
for(i in 1:dim(img)[1]){
  for(j in 1:dim(img)[2]){
    for(k in 1:dim(img)[3]){
      big_img[((i-1)*resmul+1):(i*resmul),
              ((j-1)*resmul+1):(j*resmul),
              k] <- img[i,j,k]
    }
  }
}
img <- big_img

#buffer image (helps with later simplification)
nbuffers <- 3
for(i in 1:nbuffers){
  img <- abind::abind(img[1,,], img, along = 1)
  img <- abind::abind(img, img[nrow(img),,], along = 1)
  img <- abind::abind(img[,1,], img, along = 2)
  img <- abind::abind(img, img[,ncol(img),], along = 2)
}

#discreetize and simplify colors
#TODO just write own version of weighted kmeans
# discrete_img <- discretize_image(img, ncol_range = 3:10)
discrete_img <- discretize_image_2(img, ncol_range = 3:10, offset = 0)

cluster_inds <- discrete_img$cluster_inds
mean_cols <- discrete_img$mean_cols
centers <- discrete_img$centers
#paint transparent colors white (for now)
mean_cols[centers[,4] < 0.2] <- rgb(col2rgb("white")[1], col2rgb("white")[2], col2rgb("white")[3], maxColorValue = 255)

processed_img <- abind::abind(lapply(1:4, function(i) matrix(centers[cluster_inds,i], nrow(img), ncol(img), byrow = T)), along = 3)

#simplify colorspace
#this isn't actually removing blocks of 4 pixels
#TODO can rewrite routefinder to turn around if necessary
#and only terminate if it returns to the start
clmat <- simplify_image(cluster_inds, processed_img, n_passes = 1)
cluster_inds <- c(t(clmat))
clusters_left <- sort(unique(cluster_inds))

#adjust in case we lost a color
mean_cols <- mean_cols[clusters_left]
centers <- centers[clusters_left,]
ncols <- length(clusters_left)
for(i in 1:ncols){
  cluster_inds[cluster_inds == clusters_left[i]] <- i
  clmat[which(clmat == clusters_left[i])] <- i
}

#specify colors for shape borders
border_cols <- rep("black", length(mean_cols))

#plot early result
processed_img <- abind::abind(lapply(1:4, function(i) 
  matrix(centers[cluster_inds,i], nrow(processed_img), ncol(processed_img), byrow = T)), along = 3)

#retrieve shape data and simplify
onehot <- lapply(1:ncols, function(i) (clmat == i))
processed_shapes <- separate_layers(onehot, clmat, use_wavelet = T, kconv_init = 1)
old_shapes <- processed_shapes$old_shapes
new_shapes <- processed_shapes$new_shapes

#find shape plotting order
shape_plotting_order <- find_shape_plotting_order(new_shapes = new_shapes)
plotting_order <- shape_plotting_order$plotting_order
new_shape_inds <- shape_plotting_order$new_shape_inds

#plot reconstructed shapes
par(xpd = NA)
plot(0,0,col = 0, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1), 
     axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n')

#plot initial image
rasterImage(img,0,0,1,1)

#plot color-simplified image
rasterImage(processed_img,0,-1,1,0)

#plot old shapes
s2p <- new_shape_inds[plotting_order[1],]
xyr <- apply(new_shapes[[s2p[1]]][[s2p[2]]], 2, range)
for(i in 2:length(plotting_order)){
  s2p <- new_shape_inds[plotting_order[i],]
  xs <- old_shapes[[s2p[1]]][[s2p[2]]][,1]
  ys <- old_shapes[[s2p[1]]][[s2p[2]]][,2]
  polygon(x = (xs - xyr[1,1]) / diff(xyr[,1]) - 1, y = (ys - xyr[1,2]) / diff(xyr[,2]) - 1, 
          col = mean_cols[new_shape_inds[plotting_order[i],1]],
          border = border_cols[new_shape_inds[plotting_order[i],1]])
}

#plot new shapes
for(i in 2:length(plotting_order)){
  s2p <- new_shape_inds[plotting_order[i],]
  xs <- new_shapes[[s2p[1]]][[s2p[2]]][,1]
  ys <- new_shapes[[s2p[1]]][[s2p[2]]][,2]
  polygon(x = (xs - xyr[1,1]) / diff(xyr[,1]) - 1, y = (ys - xyr[1,2]) / diff(xyr[,2]), 
          col = mean_cols[new_shape_inds[plotting_order[i],1]],
          border = border_cols[new_shape_inds[plotting_order[i],1]])
}

#add in some arrows
arrow_cols <- table(cluster_inds)[mean_cols != "#FFFFFF"]
arrow_cols <- mean_cols[as.numeric(names(arrow_cols)[order(arrow_cols, decreasing = T)[1:2]])]
arrow_font <- "Courier"
grad_arrow_curve(c(1.1,1.1,0.25,-0.25), cols = arrow_cols, direc = "v", w = 0.05, 
                 prop_head_width = 3, prop_shaft_length = 0.15,
                 nslices = 50)
text(x = 1.15, y = 0 + 0.15 * 0.5 / 2, srt = 270, labels = "simplify colors", cex = 0.75,
     family = arrow_font)

xyrat <- diff(range(par("usr")[1:2])) / par("pin")[1] / diff(range(par("usr")[3:4])) * par("pin")[2]
grad_arrow_curve(c(0.25 * xyrat, -0.25 * xyrat, -1.1, -1.1), cols = arrow_cols, direc = "h", w = 0.05 / xyrat, 
                 prop_head_width = 3, prop_shaft_length = 0.15,
                 nslices = 50)
text(y = -1.1 - 0.05 / xyrat, x = 0 + 0.15 * 0.5 * xyrat / 2, labels = "extract shapes", cex = 0.75,
     family = arrow_font)

grad_arrow_curve(c(-1.1,-1.1,-0.25,0.25), cols = arrow_cols, direc = "v", w = 0.05, 
                 prop_head_width = 3, prop_shaft_length = 0.15,
                 nslices = 50)
text(x = -1.15, y = 0 - 0.15 * 0.5 / 2, srt = 90, labels = "smooth shapes", cex = 0.75,
     family = arrow_font)


#### debugging and other exploratory code code ####

# par(mfrow = find_factors(length(new_shapes)) + c(0,1))
plot(1:length(mean_cols), rep(1, length(mean_cols)), pch = 19, col = 1, cex = 5.1)
points(1:length(mean_cols), rep(1, length(mean_cols)), pch = 19, col = mean_cols, cex = 5.0)
text(1:length(mean_cols), rep(1, length(mean_cols)), lwd = 2, 
     col = sapply(mean_cols, inverse_hex_color), cex = 2.0)

img_dims <- apply(do.call(rbind, lapply(new_shapes, function(layer) 
  do.call(rbind, lapply(layer, function(shp) apply(shp, 2, range))))), 2, range)
for(li in 1:length(new_shapes)){
  plot(NULL, xlim = img_dims[,1], ylim = img_dims[,2])
  for(si in 1:length(new_shapes[[li]])){
    polygon(as.matrix(new_shapes[[li]][[si]]), 
            col = adjustcolor(mean_cols[li], 0.3))
  top_pt <- which.max(new_shapes[[li]][[si]]$y)
    text(new_shapes[[li]][[si]]$x[top_pt], 
         new_shapes[[li]][[si]]$y[top_pt] + diff(img_dims[,2]) / 30, 
            labels = paste0("(l", li, ", s", si, ", i", 
            which(new_shape_inds[,1] == li & new_shape_inds[,2] == si), ")"), cex = 1, xpd = NA)
  }
}

for(li in 1:length(old_shapes)){
  plot(NULL, xlim = img_dims[,1], ylim = img_dims[,2])
  for(si in 1:length(old_shapes[[li]])){
    polygon(as.matrix(old_shapes[[li]][[si]]), 
            col = adjustcolor(mean_cols[li], 0.3))
    top_pt <- which.max(old_shapes[[li]][[si]]$y)
    text(old_shapes[[li]][[si]]$x[top_pt], 
         old_shapes[[li]][[si]]$y[top_pt] + diff(img_dims[,2]) / 30, 
         labels = paste0("(l", li, ", s", si, ", i", 
                         which(new_shape_inds[,1] == li & new_shape_inds[,2] == si), ")"), cex = 1, xpd = NA)
  }
}

plot(which(onehot[[1]], T) %*% matrix(c(0,1,-1,0), 2, 2), col = adjustcolor(mean_cols[1], 0.5), pch = 19)
for(layer in 2:length(onehot)){
  points(which(onehot[[layer]], T) %*% matrix(c(0,1,-1,0), 2, 2), col = adjustcolor(mean_cols[layer], 0.5), pch = 19)
}


plot(old_shapes[[1]][[1]]$x, -old_shapes[[1]][[1]]$y, type = "l")
for(layer in sample(1:length(onehot))){
  for(i in 1:length(old_shapes[[layer]])){
    polygon(old_shapes[[layer]][[i]]$x, -old_shapes[[layer]][[i]]$y, col = adjustcolor(mean_cols[layer], 1))
  }
}

plot(new_shapes[[1]][[1]]$x, -new_shapes[[1]][[1]]$y, type = "l")
for(layer in sample(1:length(onehot))){
  for(i in 1:length(new_shapes[[layer]])){
    polygon(new_shapes[[layer]][[i]]$x, -new_shapes[[layer]][[i]]$y, col = adjustcolor(mean_cols[layer], 1))
  }
}

par(mfrow = c(2,2))
plot(shapes[[1]]$x, -shapes[[1]]$y, type = "l")
for(i in 1:length(shapes)){
  polygon(shapes[[i]]$x, -shapes[[i]]$y, col = adjustcolor(mean_cols[layer], 1))
}

plot(shapes[[1]]$x, -shapes[[1]]$y, type = "l")
for(i in 1:length(shapes)){
  polygon(shapes[[i]]$x, -shapes[[i]]$y, col = adjustcolor(2, 0.2))
}
for(i in 1:length(fshapes)){
  polygon(fshapes[[i]]$x, -fshapes[[i]]$y, col = adjustcolor(3, 0.2))
}

plot(fshapes[[1]]$x, -fshapes[[1]]$y, type = "l")
for(i in 1:length(fshapes)){
  polygon(fshapes[[i]]$x, -fshapes[[i]]$y, col = adjustcolor(mean_cols[layer], 1))
}

#####

#diy

dfft <- nrow(kern) + nrow(cmat) - 1
pcmat <- pad(x=cmat, n=dfft-nrow(cmat))
pkern1 <- pad(x=kern, n=dfft-nrow(kern))
pkern2 <- pad(x=t(kern), n=dfft-nrow(kern))
cmat_conv_list <- list(Re(fft(fft(pcmat) * fft(pkern1), inverse = T)), 
                  Re(fft(fft(pcmat) * fft(pkern2), inverse = T)))
cmat_conv_list <- lapply(cmat_conv_list, scale2)


cmat_conv_spatial_list <- replicate(2, cmat * 0, F)
for(i in 1:nrow(cmat_sdc)){
  for(j in 1:ncol(cmat_sdc)){
    input_cmat <- pcmat[i:(i+nrow(kern)-1), j:(j+ncol(kern)-1)]
    cmat_conv_spatial_list[[1]][i, j] <- sum(kern * input_cmat) 
    cmat_conv_spatial_list[[2]][i, j] <- sum(t(kern) * input_cmat)
  }
}

heatmap(cmat_conv_spatial_list[[1]], NA, NA)
heatmap(scale2(trim(cc(cmat_conv_list[[1]]), 2)), NA, NA)

cc <- function(x){ #for center corner
  cn <- ceiling(dim(x) / 2)
  y <- rbind(cbind(x,x),
        cbind(x,x))
  y[(cn[1]+1):(3*cn[1]), (cn[2]+1):(3*cn[2])]
}

center_corner(cmat_conv_list[[2]])



