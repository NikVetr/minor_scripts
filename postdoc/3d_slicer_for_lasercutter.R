library(rgl)
options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)

#top level parameters
frame_dim <- c(w = 2, h = 1)
buffer_x <- 0.5
xlim <- c(0, frame_dim["w"])
ylim <- c(0, frame_dim["h"])
plot_out_of_frame <- F
plot_3D <- T

#wavy lines
yf <- function(x){
 list(
   sin(x * 3  - 0) / 2 + (x - 1 / 4) / 2,
   x * sin(x * 3) / 4 + 1/2 + x / 2,
   (x + 1) ^ (1 / 2) * sin(x * 3 + 1) / 3 - 1/2 + (max(x) - x / 4) / 5
 ) 
}

#spiral shape
use_spiral <- F
xyf <- function(x){
  rx <- range(x)
  midp_x <- min(rx) + diff(rx) * 0.65
  midp_y <- 0.75
  rad_x <- diff(rx) / 2
  t <- seq(0, 8 * pi, length.out = length(x))
  rads <- seq(0, rad_x, length.out = length(x))
  x <- sin(t) * rads + midp_x
  y <- cos(t) * rads + midp_y
  return(list(x = x, y = y))
}

#return a matrix of z values for all points in grid given by limits

#first specify bivariate density functions about lines
invlogit <- function(p) exp(p) / (1 + exp(p))

z_min_prop <- 0.1
z_max_prop <- 0.5
logdens <- list(
  function(x, y, xl, yl){dnorm(x = x, mean = xl, sd = 0.1) * 
      dnorm(x = y, mean = yl, sd = 0.1) * x * (sin(x * 10) / 100 + 1)},
  function(x, y, xl, yl){dnorm(x = x, mean = xl, sd = 0.1) * 
      dnorm(x = y, mean = yl, sd = 0.1 * x ^ 0.25)},
  function(x, y, xl, yl){dnorm(x = x, mean = xl, sd = 0.1) * 
      dnorm(x = y, mean = yl, sd = 0.1)}
)

logdens <- list(
  function(x, y, xl, yl){dexp(x = abs(x - xl), rate = 10) * 
      dexp(x = abs(y - yl), rate = 10) / 3 * x * (sin(x * 10) / 100 + 1) *
      invlogit(-abs(y - 1 / 2) * 4 *
      abs(x - 2 / 2))},
  function(x, y, xl, yl){dnorm(x = x, mean = xl, sd = 0.1) * 
      dnorm(x = y, mean = yl, sd = 0.1 * x ^ 0.25) *
      invlogit(-abs(y - 1 / 2) * 4 *
      abs(x - 2 / 2))},
  function(x, y, xl, yl){dnorm(x = x, mean = xl, sd = 0.1) * 
      dnorm(x = y, mean = yl, sd = 0.1) *
      invlogit(-abs(y - 1 / 2) * 4 *
      abs(x - 2 / 2))}
)

#then write function for getting grid densities
zf <- function(xy, xlim, ylim, logdens, res = 1E2, z_max = NULL){
  
  n_lines <- length(xy)
  
  #construct grid over range of values
  xy_rat <- diff(range(xlim)) / diff(range(ylim))
  if(xy_rat > 1){
    x_grid <- seq(xlim[[1]], xlim[[2]], length.out = ceiling(res * xy_rat))
    y_grid <- seq(ylim[[1]], ylim[[2]], length.out = res)
  } else {
    x_grid <- seq(xlim[[1]], xlim[[2]], length.out = res)
    y_grid <- seq(ylim[[1]], ylim[[2]], length.out = ceiling(res / xy_rat))
  }
  xy_grid <- expand.grid(x = x_grid, y = y_grid)
  xy_grid_inds <- expand.grid(xi = 1:length(x_grid), yi = 1:length(y_grid))
  
  #compute densities at points
  ptdens <- Reduce("+", lapply(1:n_lines, function(li){
    lf <- logdens[[li]]
    apply(apply(xy_grid, 1, function(xy_pt){
      lf(x = xy_pt["x"], y = xy_pt["y"], 
         xl = xy[[li]][,"x"], yl = xy[[li]][,"y"])
    }), 2, sum)
  }))
  xw <- diff(x_grid[1:2])
  yh <- diff(y_grid[1:2])
  if(is.null(z_max)){
    total_volume <- sum(xw * yh * ptdens)
    z <- ptdens / total_volume  
  } else {
    total_volume <- sum(xw * yh * ptdens)
    z <- ptdens / max(ptdens) * z_max
  }
  
  xyz_df <- cbind(xy_grid, z = z)
  
  #get heights into matrix format?
  x_mat <- matrix(data = xy_grid[,"x"], ncol = length(x_grid),
                  nrow = length(y_grid), byrow = T)
  y_mat <- matrix(data = xy_grid[,"y"], ncol = length(x_grid),
                  nrow = length(y_grid), byrow = T)
  z_mat <- matrix(data = z, ncol = length(x_grid),
                  nrow = length(y_grid), byrow = T)
  
  return(list(xyz_df = xyz_df, 
              x_grid = x_grid, 
              y_grid = y_grid,
              xw = xw, yh = yh,
              xy_grid_inds = xy_grid_inds,
              x_mat = x_mat,
              y_mat = y_mat,
              z_mat = z_mat))
}

#compute 2d line, equal intervals in horizontal dimension
res_line <- 3E2
x <- seq(-buffer_x, frame_dim["w"] + buffer_x, length.out = res_line * 1E2)
if(use_spiral){
  xy <- xyf(x)
  x <- xy[["x"]]
  y <- xy["y"]
} else {
  y <- yf(x)
}
n_lines <- length(y)

#respecify line to have equal intervals in length
xy <- lapply(1:n_lines, function(li){
  yl <- y[[li]]
  seg_lens <- sqrt(diff(x)^2 + diff(yl)^2)
  line_len <- sum(seg_lens)
  interval_len <- line_len / (res_line - 1)
  seg_inds <- floor(cumsum(seg_lens) / interval_len)
  n_inds <- sapply(split(seg_inds, seg_inds), length)[1:(res_line - 1)]
  ind_locs <- c(1, cumsum(n_inds))
  cbind(x = x[ind_locs],
        y = yl[ind_locs])
})

#2d view of lines
par(mfrow = c(2,1), mar = c(0,0,0,0))

#set up window
if(plot_out_of_frame){
  ylim_window <- range(unlist(y))  
} else {
  ylim_window <- ylim
}

plot2d_blank <- function() plot(NA, NA, xlim = xlim, ylim = ylim_window, asp = 1, 
     axes = F, xlab = "", ylab = "", frame.plot	= F, xpd = NA)
plot2d_blank()
for(i in 1:n_lines){
  if(plot_out_of_frame){
    valid_inds <- rep(T, res_line)
  } else {
    valid_inds <- xy[[i]][,"y"] > ylim[1] & xy[[i]][,"y"] < ylim[2] &
                  xy[[i]][,"x"] > xlim[1] & xy[[i]][,"x"] < xlim[2]  
  }
  
  x_plot <- xy[[i]][valid_inds,"x"]
  y_plot <- xy[[i]][valid_inds,"y"]
  lines(x_plot, y_plot, lwd = 2, col = i)
}
rect(xlim[1], ylim[1], xlim[2], ylim[2], lwd = 4)

#compute densities
z_min <- -min(frame_dim) * z_min_prop
z_max <- min(frame_dim) * z_max_prop
xyz_info <- zf(xy, xlim, ylim, logdens, z_max = z_max, res = 40)
xyz_df <- xyz_info$xyz_df
z_mat <- xyz_info$z_mat
n_cells <- nrow(xyz_df)

#plot heatmap of densities?
ncol <- 100
col_inds <- ceiling((xyz_df$z - min(xyz_df$z)) / diff(range(xyz_df$z)) * (ncol - 1)) + 1
xyz_df$col <- viridisLite::rocket(ncol)[col_inds]
plot2d_blank()
rect(xleft = xyz_df$x - xyz_info$xw / 2,
     xright = xyz_df$x + xyz_info$xw / 2,
     ybottom = xyz_df$y - xyz_info$yh / 2,
     ytop = xyz_df$y + xyz_info$yh / 2,
     border = NA, col = xyz_df$col)

#now plot 3D slices of the manifold (in x direction for now)
if(plot_3D){
  
#first compute the slices
n_slices <- 40
slice_inds <- round(seq(1, ncol(z_mat), length.out = n_slices))
slices <- lapply(1:n_slices, function(si){
  sind <- slice_inds[si]
  slice_top <- cbind(x = xyz_info$x_mat[,sind],
                     y = xyz_info$y_mat[,sind],
                     z = xyz_info$z_mat[,sind])
  slice_start <- slice_top[1,] * c(1,1,0) + c(0,0,z_min)
  slice_end <- tail(slice_top, 1) * c(1,1,0) + c(0,0,z_min)
  slice_verts <- rbind(slice_start, slice_top, slice_end, slice_start)
  return(slice_verts)
})

#then plot them (first as lines)
open3d()
bg3d("white")
for(si in 2:n_slices){
  lines3d(slices[[si]], col = "black", lwd = 1.5)  
}
lines3d(slices[[1]], col = "black", lwd = 1.5)

#then as polygons
poly2quads <- function(poly, z_min = 0){
  nverts <- nrow(poly)
  poly_bottom <- poly_top <- poly[-c(1, nverts - 0:1),]
  poly_bottom[,"z"] <- z_min
  nquads <- nverts - 4
  quads <- lapply(1:nquads, function(qi){
    rbind(poly_top[qi,],
          poly_top[qi+1,],
          poly_bottom[qi+1,],
          poly_bottom[qi,]
          )
  })
  return(quads)
}

plot_polyquads <- function(quads){
  for(qi in 1:length(quads)){
    quads3d(quads[[qi]], col = "grey40", alpha = 1)  
  }
}

clear3d()
for(si in 2:n_slices){
  plot_polyquads(poly2quads(slices[[si]], z_min = z_min))
  lines3d(slices[[si]], col = "white", lwd = 1.5)  
}
plot_polyquads(poly2quads(slices[[1]], z_min = z_min))
lines3d(slices[[1]], col = "white", lwd = 0)

#give them some thickness?



}
