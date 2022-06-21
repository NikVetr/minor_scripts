#packages
library(Biostrings)
library(graphics)
library(foreach)
library(doParallel)
library(parallel)


#functions
intersecting_point <- function(xy0a, xy0b, xy1a, xy1b){
  mb0 <- solve(cbind(rbind(xy0a[1], xy0b[1]), c(1,1))) %*% rbind(xy0a[2], xy0b[2])
  mb1 <- solve(cbind(rbind(xy1a[1], xy1b[1]), c(1,1))) %*% rbind(xy1a[2], xy1b[2])
  intspt <- rev(as.vector(solve(cbind(rbind(1,1), rbind(-mb0[1], -mb1[1]))) %*% rbind(mb0[2], mb1[2])))
  return(intspt)
}
closestAngle <- function(focal, thetas){
  #return index of thetas that is closest to focal
  focal <- coerceTo2Pi(focal)
  thetas <- coerceTo2Pi(thetas)
  signed_distances <- rbind(focal - thetas, 2*pi + focal - thetas, -2*pi + focal - thetas)
  distances <- abs(signed_distances)
  
  which_mins <- apply(distances, 2, which.min)
  directions <- -sign(sapply(1:length(thetas), function(md) signed_distances[which_mins[md], md]))
  distances <- apply(distances, 2, min) 
  
  return(list(index = which.min(distances), distance = min(distances), direction = directions[which.min(distances)]))
}
coerceTo2Pi <- function(theta){
  if(length(theta) > 1){
    return(sapply(theta, function(x) coerceTo2Pi(x)))
  }
  if(theta > 0){
    return(theta %% (2*pi))
  } else {
    num_2pis_to_add <- ceiling(abs(theta / (2*pi)))
    return(theta + num_2pis_to_add * 2 * pi)
  }
}
smooth0to0.5 <- function(x){
  (sin(x*pi*2-pi/2) + 1) / 4
}
logistic <- function(x) 1 / (1+exp(-x))
inv_logistic <- function(x) log(x / (1-x))
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
polar2cart <- function(t, r){
  return(c(r*cos(t), r * sin(t)))
}
cart2polar <- function(x,y){
  return(c(atan2(x = x, y = y), sqrt(x^2 + y^2)))
}
rot <- function(x,y,t){
  return(matrix(c(cos(t), -sin(t), sin(t), cos(t)),2,2, byrow = T) %*% t(t(c(x,y))))
}
draw_circle <- function(x,y,r,...){
  angles <- seq(0,2*pi,length.out = 200)
  verts <- t(sapply(angles, function(a) polar2cart(a, r) + c(x,y)))
  polygon(verts[,1], verts[,2], ...)
}

# exploring smoothing functions
# buffer <- 0
# x <- c(rep(1, buffer), 1:50, 49:1, rep(1,buffer))
# plot(x[(buffer+1):(buffer+99)], type = "l")
# r <- 0.5
# #exponentially weighted smoothing
# y <- sapply(1:99, function(i) 
#   sum(c(x[i:99] * dexp(0:99, r)[1:(99-i+1)], 
#         x[(i-1):1] * dexp(1:99, r)[1:(i-1)])) / 
#     sum(c(dexp(0:99, r)[1:(99-i+1)], dexp(0:99, r)[1:(i-1)])))
# lines(y, col = "purple")
# #sliding window
# buffer <- 10
# x <- c(rep(1, buffer), 1:50, 49:1, rep(1,buffer))
# w <- 10
# y <- sapply((buffer+1):(buffer+99), function(i) mean(x[(i-w):(i+w)]))
# lines(y, col = "blue")


n_vert <- 8
r <- 2
init_vert <- c(0,r)
angles <- seq(cart2polar(init_vert[1], init_vert[2])[1], cart2polar(init_vert[1], init_vert[2])[1] + 2*pi, 
              length.out = n_vert + 1)
verts <- t(sapply(angles, function(a) polar2cart(a, r)))
verts <- t(sapply(1:nrow(verts), function(v) rot(verts[v,1], verts[v,2], pi/16)))
verts <- t(sapply(1:nrow(verts), function(v) rot(verts[v,1], verts[v,2], -pi/16)))

#specify colors
cols <- viridis::cividis(n_vert)

#generate basic polygon
plot(1,1,col = "white", xlim = c(-2,2), ylim = c(-2,2), xaxt = "n", yaxt = "n", frame.plot = F, xlab = "", ylab = "")
polygon(verts[,1], verts[,2], col = cols[1])

#inscribe smaller polygon inside
for(i in 2:(n_vert-2)){
  n_vert <- n_vert - 1
  init_angle <- mean(c(angles[1], angles[2]))
  r <- cos(abs(angles[1] - init_angle)) * r
  angles <- seq(init_angle, init_angle + 2*pi, length.out = n_vert + 1)
  verts <- t(sapply(angles, function(a) polar2cart(a, r)))
  polygon(verts[,1], verts[,2], col = cols[i])
}


#### generate inscribed same polygons ####
n_vert <- 6
n_inscr_polygons <- 10
cols <- viridis::cividis(n_inscr_polygons)
init_r <- r <- 1.75
init_vert <- polar2cart(pi/8, init_r)
angles <- seq(cart2polar(init_vert[1], init_vert[2])[1], cart2polar(init_vert[1], init_vert[2])[1] + 2*pi, 
              length.out = n_vert + 1)
verts <- t(sapply(angles, function(a) polar2cart(a, r)))
verts <- t(sapply(1:nrow(verts), function(v) rot(verts[v,1], verts[v,2], pi/16)))
verts <- t(sapply(1:nrow(verts), function(v) rot(verts[v,1], verts[v,2], -pi/16)))

par(mar = c(0,0,0,0), xpd = NA)
plot(1,1,col = "white", xlim = c(-2,2), ylim = c(-2,2), xaxt = "n", yaxt = "n", frame.plot = F, xlab = "", ylab = "")

#draw background with stars
rect(-1E6, -1E6, 1E6, 1E6, col = "#000a1a")
n_stars <- 4E2
star_pos <- matrix(runif(n_stars*2, -2.5, 2.5), n_stars, 2)
star_size <- runif(n_stars,0,1.5)
star_cols <- sample(cols, n_stars, replace = T)
star_alpha <- star_size <- runif(n_stars,0,1)
star_cols_alpha <- sapply(1:n_stars, function(s) adjustcolor(star_cols[s], alpha.f = star_alpha[s]))
points(star_pos[,1], star_pos[,2], pch = 19, col = star_cols_alpha, cex = star_size, xpd = NA)

#draw initial polygon
polygon(verts[,1], verts[,2], col = cols[1])
verts_init <- verts

#put circle on border
angle_circle_init <- angle_circle <- runif(1,0,2*pi)
angle_circle_init <- angle_circle <-  0 

n_chains <- 3
# for(j in 1:n_chains){
for(j in 1:1){
if(j > 1){angle_circle <- angle_circle_init <- (angle_circle_init + 2*pi/n_chains) %% (2*pi)}
n_circles <- 5
prop_angle_skip_edge <- 0.1
cols_circle <- rev(viridis::cividis(n_circles, begin = 0.5))
r_circles <- 0.25 * 1 / 1.25^(0:n_circles)
# for(i in 1:n_circles){
for(i in 1:1){
  r_circle <- r_circles[i]
  if(i > 1){
    sum_rads <- sum(r_circles[c(i, i-1)])
    touch_r_next <- sqrt(touch_r^2 + sum_rads^2 - 2 * touch_r * sum_rads * cos(pi - side_angle - (max(angles[adj_verts]) - angle_circle)))
    angle_circle <- angle_circle + asin(sin(pi - side_angle - (max(angles[adj_verts]) - angle_circle)) / touch_r_next * sum_rads)
  }
  angles_std <- angles %% (2 * pi)
  side_angle <- (pi - 2*pi/n_vert) / 2
  weight_angle_thresh <- prop_angle_skip_edge * 2*pi/n_vert
  adj_verts <- c(max(which(angles <= angle_circle)), 
                 min(which(angles >= angle_circle)))
  if(any(adj_verts == Inf | adj_verts == -Inf)){
    adj_verts <- c(max(which(angles <= (2*pi + angle_circle))), 
                  min(which(angles >= (2*pi + angle_circle))))
  }
  # angle_circle
  # angles[adj_verts]
  # angles_std[adj_verts]
  
  touch_r <- init_r * sin(side_angle) / abs(sin(pi - side_angle - (max(angles[adj_verts]) - angle_circle)))
  touch_pt <- polar2cart(angle_circle, touch_r)
  
  adj_angle <- angles_std[adj_verts[which.min(abs(angles_std[adj_verts] - angle_circle))]]
  adj_angle <- adj_angle + sign(adj_angle - angle_circle) * pi/n_vert
  adj_angle <- adj_angle %% (2 * pi)
  
  prop_adj_angle <- max(weight_angle_thresh - abs(mean(c(adj_angle, mean(angles_std[adj_verts]) %% (2*pi))) - angle_circle), 0) / weight_angle_thresh / 2
  
  if(sign(diff(angles_std[adj_verts]))){
    angle_out <- mean(c(angles_std[adj_verts[1]]-2*pi, angles_std[adj_verts[2]])) * (1 - prop_adj_angle) + adj_angle * prop_adj_angle
  } else {
    angle_out <- mean(angles_std[adj_verts]) * (1 - prop_adj_angle) + adj_angle * prop_adj_angle  
  }
  
  # cent_circ <- polar2cart(angle_out - ifelse(sign(diff(angles_std[adj_verts])) == -1, pi, 0), r_circle) + touch_pt
  cent_circ <- polar2cart(angle_out, r_circle) + touch_pt
  
  draw_circle(cent_circ[1], cent_circ[2], r_circle, col = cols_circle[i], lwd = 0.1, border = cols[1])
  points(touch_pt[1], touch_pt[2])
  points(cent_circ[1], cent_circ[2], pch = 19, col = cols[1])
  if(i > 1){
    segments(x0 = cent_circ[1], y0 = cent_circ[2], x1 = prev_cent_circ[1], y1 = prev_cent_circ[2], col = cols[1], lwd = 1)
  }
  
  prev_cent_circ <- cent_circ
  
  
}

print(angle_circle)
text(x = cent_circ[1], y = cent_circ[2], labels = j, col = "white")
}

# #inscribe smaller polygon inside
for(i in 2:(n_inscr_polygons)){
  n_vert <- n_vert
  init_angle <- mean(c(angles[1], angles[2]))
  r <- cos(abs(angles[1] - init_angle)) * r
  angles <- seq(init_angle, init_angle + 2*pi, length.out = n_vert + 1)
  verts <- t(sapply(angles, function(a) polar2cart(a, r)))
  # polygon(verts[,1], verts[,2], col = cols[i], border = cols[n_inscr_polygons - i], lwd = 0.1)
  polygon(verts[,1], verts[,2], col = cols[i], border = ifelse(i != n_inscr_polygons, NA, cols[1]), lwd = 2)
  #draw lines connecting vertices to center
  if(i != n_inscr_polygons){
    for(i in 1:nrow(verts)){
      segments(x0 = 0, y0 = 0, x1 = verts[i,1], y1 = verts[i,2], col = cols[1], lwd = 2)
    }
  }

}

# ~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~ #
#### make animation ####
# ~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~ #

set.seed(777)

file.remove(paste0("~/Documents/geometric_animation/frames/", list.files(path = "~/Documents/geometric_animation/frames/", pattern = ".png")))

n_runs_reset <- 10
if(!exists("n_runs")){
  n_runs <- 0
}

if(!exists("cl")){
  cl <- makeCluster(16, outfile="")
  registerDoParallel(cl)
  n_runs <- 0
}
getDoParWorkers()

#specify animation-specific parameters
debug_animation <- F
render_video <- T
recalculate_lights <- T

pwidth <- 1080
pheight <- 1080
framethin <- 1
nfps <- 60
nsec = 30
nf <- (nsec * nfps)

#specify star properties
n_stars <- 2E3
star_pos_init <- matrix(runif(n_stars*2, -3, 3), n_stars, 2)
star_size <- runif(n_stars,0,1.5)
star_cols <- sample(cols, n_stars, replace = T)
star_size <- runif(n_stars,0,1)
star_alphas <- matrix(NA, ncol = nf, nrow = n_stars)
star_alphas_optima <- star_alphas[,1] <- inv_logistic(runif(n_stars,0,1))
OU_var <- 1 / nfps * framethin
OU_RB <- 1 / nfps * framethin
for(frame in 2:nf){
  star_alphas[,frame] <- star_alphas[,frame-1] + rnorm(n_stars, 0, sqrt(OU_var)) + (star_alphas_optima - star_alphas[,frame-1]) * OU_RB
}
star_alphas <- logistic(star_alphas)
  
#specify polygon properties
r <- 2
n_vert <- 6
starting_angle <- pi/8
n_inscr_polygons <- 10
light_amplitude <- 0.5
light_period_in_sec <- 1.5
light_period_in_frames <- light_period_in_sec * nfps
max_light_val <- logistic(light_amplitude * sin(2*pi/light_period_in_frames*1:nf) + 4 - logistic((1:nf - nf/2)/nf*7) * 3)
max_opacity_val <- sqrt(max_light_val)
min_light_val <- logistic(light_amplitude * sin(2*pi/light_period_in_frames*1:nf) + -2 - logistic((1:nf - nf/2)/nf*7) * 2)
plot(max_light_val, type = "l", col = "red", ylim = c(0,1)); lines(min_light_val, col = "blue")
cols_matrix <- t(sapply(1:nf, function(f1) viridis::cividis(n_inscr_polygons, begin = min_light_val[f1], end = max_light_val[f1])))

init_r <- 1.75
polygon_opacities <- 0.9 + 1:n_inscr_polygons / 100

#specify rotation speeds
worm_speed <- -0.05 #rotations per second
worm_speed <- worm_speed * 2 * pi / nfps #radians per frame
# polygon_speeds <- 0.05 #rotations per second
polygon_speeds <- runif(n_inscr_polygons, 0.03, 0.15) #rotations per second
polygon_speeds <- polygon_speeds * 2 * pi / nfps #radians per frame
star_speed <- 0.0005

##specify worm properties
n_chains <- 3
n_circles <- 5
prop_angle_skip_edge <- 0.2
cols_circle <- rev(viridis::cividis(n_circles, begin = 0.5))
r_circles <- 0.25 * 1 / 1.25^(0:n_circles)

#specify line properties
n_extra_cols <- 5
conduit_cols <- c(cols[1], viridis::plasma(n_extra_cols, begin = 0.75, end = 1))
transfer_threshold = 0.02
recharge_rate <- 0.1 * nfps
if(recalculate_lights | !exists("angle_matrix")){
  conduit_cols_array <- array(data = 1, dim = c(n_inscr_polygons, n_vert, nf)) #need to pre-roll in order to parallelize
  wormfood_array <- array(data = 0, dim = c(n_chains, n_circles, nf)) #need to pre-roll in order to parallelize
  for(f in 2:nf){
    
    if((f*100) %% nf == 0){cat(paste0(f*100/nf, " "))}
    
    
    #work out angles of vertices
    angle_matrix <- matrix(data = NA, nrow = n_inscr_polygons, ncol = n_vert)
    init_vert <- polar2cart(pi/8 + polygon_speeds[1] * f, init_r)
    angles <- seq(cart2polar(init_vert[1], init_vert[2])[1], cart2polar(init_vert[1], init_vert[2])[1] + 2*pi, 
                  length.out = n_vert + 1)
    angle_matrix[1,] <- angles[-length(angles)]
    for(i in 2:(n_inscr_polygons)){
      
      # init_angle <- mean(c(coerceTo2Pi(angles[1]), coerceTo2Pi(angles[2])))
      # # if(i == 2){print(init_angle)}
      # r <- cos(abs(angles[1] - init_angle)) * r
      # angles <- seq(init_angle, init_angle + 2*pi, length.out = n_vert + 1)
      # verts <- t(sapply(angles, function(a) polar2cart(a, r)))
      # verts <- t(sapply(1:nrow(verts), function(v) rot(verts[v,1], verts[v,2], polygon_speeds[i] * f1 * ifelse(i %% 2 == 0, -2, 1))))
      # angles <- sapply(1:n_vert, function(j) coerceTo2Pi(cart2polar(x = verts[j,1], y = verts[j,2])[1]))
      # angle_matrix[i,] <- angles
      
      n_vert <- n_vert
      side_angle <- (pi - 2*pi/n_vert) / 2
      third_angle <- pi - side_angle - pi/n_vert
      radius_reduction_ratio <- sin(side_angle) / sin(third_angle)
      # init_angle <- mean(c(coerceTo2Pi(angles[1]), coerceTo2Pi(angles[2])))
      init_angle <- starting_angle + pi/n_vert * (i-1)
      init_angle <- coerceTo2Pi(init_angle + polygon_speeds[i] * f * ifelse(i %% 2 == 0, -1, 1))
      # if(i == 2){print(init_angle)}
      # r <- cos(abs(angles[1] - init_angle)) * r
      r <- init_r * radius_reduction_ratio^(i-1)
      angles <- coerceTo2Pi(seq(init_angle, init_angle + 2*pi, length.out = n_vert + 1))
      
      # init_angle <- mean(c(angles[1], angles[2])) + (polygon_speeds[i]) * f * ifelse(i %% 2 == 0, -2, 1) #this line is wrong -- base if off initial angle state, not current angles?
      # angles <- seq(init_angle, init_angle + 2*pi, length.out = n_vert + 1)
      angle_matrix[i,] <- angles[-length(angles)]
      
    }
    
    #only have to do for one vertex, bc of symmetry, unless using multiple types of shapes or just one vertex at a time
    
    conduit_cols_array[,,f] <- conduit_cols_array[,,f-1]
    if(f %% recharge_rate == 0){vertex <- sample(1:n_vert, 1); conduit_cols_array[n_inscr_polygons,vertex,f] <- conduit_cols_array[n_inscr_polygons,vertex,f] + 1}
    for(i in (n_inscr_polygons-1):1){
      closest_angles <- sapply(1:n_vert, function(v) closestAngle(angle_matrix[i,v], angle_matrix[i+1,]))[1:2,]
      closest_angles <- rbind(closest_angles, donation = as.numeric(closest_angles[2,] < transfer_threshold)) #donation
      conduit_cols_array[i,,f] <- conduit_cols_array[i,,f] + unlist(closest_angles[3,]) * as.numeric(conduit_cols_array[i+1,unlist(closest_angles[1,]),f-1] > 1)
      conduit_cols_array[i+1,unlist(closest_angles[1,]),f] <- conduit_cols_array[i+1,unlist(closest_angles[1,]),f] - 
        unlist(closest_angles[3,]) * as.numeric(conduit_cols_array[i+1,unlist(closest_angles[1,]),f-1] > 1)
    }
    
    #let worms eat light
    worm_angles_matrix <- matrix(data = NA, nrow = n_chains, ncol = n_circles) #need to pre-roll in order to parallelize
    angle_circle_init <- angle_circle <- coerceTo2Pi(pi + worm_speed *  f)
    
    #find useful values
    side_angle <- (pi - 2*pi/n_vert) / 2
    weight_angle_thresh <- prop_angle_skip_edge * 2*pi/n_vert
    
    for(j in 1:n_chains){
      
      #get later initial angles
      if(j > 1){angle_circle <- angle_circle_init <- coerceTo2Pi((angle_circle_init + 2*pi/n_chains) %% (2*pi))}
      
      for(i in 1:n_circles){
        
        #initialize circle values
        r_circle <- r_circles[i]
        
        if(i > 1){
          sum_rads <- sum(r_circles[c(i, i-1)])
          # touch_r_next <- sqrt(touch_r^2 + sum_rads^2 - 2 * touch_r * sum_rads * cos(pi - side_angle - (max(angles[adj_verts]) - angle_circle)))
          touch_r_next <- sqrt(touch_r^2 + sum_rads^2 - 2 * touch_r * sum_rads * cos(side_angle + closestAngle(angle_circle, angles)$distance))
          # angle_circle <- angle_circle + asin(sin(pi - side_angle - (max(angles[adj_verts]) - angle_circle)) / touch_r_next * sum_rads)
          angle_circle <- angle_circle + asin(sin(side_angle + closestAngle(angle_circle, angles)$distance) / touch_r_next * sum_rads)
        }
        
        worm_angles_matrix[j,i] <- coerceTo2Pi(angle_circle)
        
      }
    }
    
    wormfood_array[,,f] <- wormfood_array[,,f-1]
    for(i in 1:n_chains){
      closest_angles <- sapply(1:n_circles, function(v) closestAngle(worm_angles_matrix[i,v], angle_matrix[1,]))[1:2,]
      closest_angles <- rbind(closest_angles, donation = as.numeric(closest_angles[2,] < transfer_threshold)) #donation
      wormfood_array[i,,f] <- wormfood_array[i,,f] + unlist(closest_angles[3,]) * as.numeric(conduit_cols_array[1,unlist(closest_angles[1,]),f-1] > 1) 
      angles_to_change <- unique(unlist(closest_angles[1,]))
      unaggregated_changes <- unlist(closest_angles[3,]) * as.numeric(conduit_cols_array[1,unlist(closest_angles[1,]),f-1] > 1) 
      amounts_to_change_by <- sapply(angles_to_change, function(a2c) sum(unaggregated_changes[unlist(closest_angles[1,]) == a2c]))
      new_conduit_cols_diff <- rep(0, length(conduit_cols_array[1,,f-1]))
      new_conduit_cols_diff[angles_to_change] <- new_conduit_cols_diff[angles_to_change] - amounts_to_change_by
      new_conduit_cols <- conduit_cols_array[1,,f] + new_conduit_cols_diff
      conduit_cols_array[1,,f] <- new_conduit_cols
    }
    
    
    
    #maybe dump a level for each circle that passes over? and drag a trail behind the worm
    # conduit_cols_array[,,f] <- matrix(data = sample(1:6, n_inscr_polygons * n_vert, replace = T, c(0.9, rep(0.1/n_extra_cols, n_extra_cols))), nrow = n_inscr_polygons)
    
  }


conduit_cols_array[conduit_cols_array > (n_extra_cols+1)] <- n_extra_cols + 1 

for(f in 1:nf){
  conduit_cols_array[,,f] <- sapply(1:n_vert, function(v) conduit_cols[as.numeric(conduit_cols_array[,v,f])])
}
plot(1:length(conduit_cols), 1:length(conduit_cols), col = conduit_cols, pch = 19, cex = 3)

}

#specify worm eating effects
wormfood_array_diffs <- array(data = 0, dim = c(n_chains, n_circles, nf)) #need to pre-roll in order to parallelize
for(f in 2:nf){wormfood_array_diffs[,,f] <- wormfood_array[,,f] - wormfood_array[,,f-1]}
n_worm_background_circles_to_draw <- array(data = 0, dim = c(n_chains, n_circles, nf)) #need to pre-roll in order to parallelize
length_propagation <- 1*nfps
ramp_up_down_ratio <- 0.1
#ramp up and then down AFTER eating
for(f in 2:nf){
  for(p in 1:length_propagation){
    if((f+p)>nf){
      next()
    }
    n_worm_background_circles_to_draw[,,f+p-1] <- n_worm_background_circles_to_draw[,,f+p-1] + wormfood_array_diffs[,,f] * (length_propagation - p + 1) 
  }
  for(p in 1:(length_propagation*ramp_up_down_ratio)){
    if((f-p)<1){
      next()
    }
    n_worm_background_circles_to_draw[,,f-p] <- n_worm_background_circles_to_draw[,,f-p] + wormfood_array_diffs[,,f] * (length_propagation * ramp_up_down_ratio - p + 1) 
  }
}
max_nbgc <- max(n_worm_background_circles_to_draw)
max_nbgc_brightness <- 0.05
max_bgc_radius_prop <- 0.35
max_bgc_radius_prop_main <- 0.15

# ~~~~~~~~~~~~~~~~~~~ #
# START THE ANIMATION #
# ~~~~~~~~~~~~~~~~~~~ #

foreach(f1=seq(1, nf, framethin), .packages = c("png")) %dopar% {
# for(f1 in seq(1, nf, framethin)){
# for(f1 in 721){

  
cols <- cols_matrix[f1,]
r <- init_r
  
if((f1*100) %% nf == 0){cat(paste0(f1*100/nf, " "))}

png(filename = paste0("~/Documents/geometric_animation/frames/frame_", 
                      paste0(rep(0, 5-nchar(((f1 - 1) / framethin) + 1)), collapse = ""), ((f1 - 1) / framethin) + 1,".png"), 
    width = pwidth, height = pheight)
par(mfrow = c(1, 1), mar = c(0,0,0,0), xpd = NA)

# work out initial polygon
init_vert <- polar2cart(starting_angle + polygon_speeds[1] * f1, init_r)
angles <- seq(cart2polar(init_vert[1], init_vert[2])[1], cart2polar(init_vert[1], init_vert[2])[1] + 2*pi, 
              length.out = n_vert + 1)
verts <- t(sapply(angles, function(a) polar2cart(a, r)))
# verts <- t(sapply(1:nrow(verts), function(v) rot(verts[v,1], verts[v,2], pi/16)))

# start plotting
plot(1,1,col = "white", xlim = c(-2,2), ylim = c(-2,2), xaxt = "n", yaxt = "n", frame.plot = F, xlab = "", ylab = "")

#draw background with stars
rect(-1E6, -1E6, 1E6, 1E6, col = "#000a1a")
star_pos <- t(sapply(1:nrow(star_pos_init), function(s) rot(star_pos_init[s,1], star_pos_init[s,2], star_speed * f1)))
star_cols_alpha <- sapply(1:n_stars, function(s) adjustcolor(star_cols[s], alpha.f = star_alphas[s,f1]))
points(star_pos[,1], star_pos[,2], pch = 19, col = star_cols_alpha, cex = star_size, xpd = NA)

#draw initial polygon
polygon(verts[,1], verts[,2], col = adjustcolor(cols[1], alpha.f = polygon_opacities[1]*max_opacity_val[f1]))
for(j in 1:(nrow(verts)-1)){
  conduit_color <- conduit_cols_array[1,j,f1]
  segments(x0 = 0, y0 = 0, x1 = verts[j,1] * 0.995, y1 = verts[j,2] * 0.995, col = conduit_color, lwd = 3)
  if(conduit_color != conduit_cols[1]){
    for(fadeout in 1:4){
      segments(x0 = 0, y0 = 0, x1 = verts[j,1] * 0.99, y1 = verts[j,2] * 0.99, col = adjustcolor(conduit_color, alpha.f = 1 / (fadeout + 5)), lwd = 3 + fadeout * 2.5)
    }
  }
}
# for(j in 1:n_vert){text(labels = paste0(i, ": ", round(coerceTo2Pi(cart2polar(x = verts[j,1], y = verts[j,2])[1]), 3)), x = verts[j,1], y = verts[j,2], col = "white")}

verts_init <- verts

#put circles / worm on border
angle_circle_init <- angle_circle <- coerceTo2Pi(pi + worm_speed *  f1)

#find useful values
angles_std <- coerceTo2Pi(angles)
side_angle <- (pi - 2*pi/n_vert) / 2
weight_angle_thresh <- prop_angle_skip_edge * 2*pi/n_vert

for(j in 1:n_chains){
  
  #get later initial angles
  if(j > 1){angle_circle <- angle_circle_init <- coerceTo2Pi((angle_circle_init + 2*pi/n_chains) %% (2*pi))}
  
  for(i in 1:n_circles){
    
    #initialize circle values
    r_circle <- r_circles[i]
    
    if(i > 1){
      sum_rads <- sum(r_circles[c(i, i-1)])
      # touch_r_next <- sqrt(touch_r^2 + sum_rads^2 - 2 * touch_r * sum_rads * cos(pi - side_angle - (max(angles[adj_verts]) - angle_circle)))
      touch_r_next <- sqrt(touch_r^2 + sum_rads^2 - 2 * touch_r * sum_rads * cos(side_angle + closestAngle(angle_circle, angles)$distance))
      # angle_circle <- angle_circle + asin(sin(pi - side_angle - (max(angles[adj_verts]) - angle_circle)) / touch_r_next * sum_rads)
      angle_circle <- angle_circle + asin(sin(side_angle + closestAngle(angle_circle, angles)$distance) / touch_r_next * sum_rads)
    }
    
    angle_circle <- coerceTo2Pi(angle_circle)
    
    adj_verts <- c(max(which(angles <= angle_circle)), 
                   min(which(angles >= angle_circle)))
    if(any(adj_verts == Inf | adj_verts == -Inf)){
      adj_verts <- c(max(which(angles <= (2*pi + angle_circle))), 
                     min(which(angles >= (2*pi + angle_circle))))
    }
    if(any(adj_verts == Inf | adj_verts == -Inf)){
      adj_verts <- c(max(which((2*pi + angles) <= (angle_circle))), 
                     min(which((2*pi + angles) >= (angle_circle))))
    }
    
    
    # angle_circle
    # angles[adj_verts]
    # angles_std[adj_verts]
    
    touch_r <- init_r * sin(side_angle) / abs(sin(pi - side_angle - (max(angles[adj_verts]) - angle_circle)))
    touch_pt <- polar2cart(angle_circle, touch_r)
    
    # #this gives angle in between vertices in closest face to current face
    # closest_vertex <- closestAngle(focal = angle_circle, thetas = angles_std)$index
    # adj_angle <- angles_std[closest_vertex]
    # adj_angle <- adj_angle + closestAngle(focal = angle_circle, thetas = angles_std)$direction * pi/n_vert
    # adj_angle <- adj_angle %% (2 * pi)
    # 
    # # prop_adj_angle <- max(weight_angle_thresh - abs(mean(c(adj_angle, mean(angles_std[adj_verts]) %% (2*pi))) - angle_circle), 0) / weight_angle_thresh / 2
    # prop_adj_angle <- (weight_angle_thresh - closestAngle(focal = angle_circle, thetas = angles_std)$distance) / (2*weight_angle_thresh)
    # prop_adj_angle <- smooth0to0.5(max(prop_adj_angle,0))
    
    #this gives angle in between vertices in closest face to current face
    closest_vertex <- closestAngle(focal = angle_circle, thetas = angles_std)$index
    adj_angle <- angles_std[closest_vertex]
    adj_angle <- adj_angle + closestAngle(focal = angle_circle, thetas = angles_std)$direction * pi/n_vert
    adj_angle <- coerceTo2Pi(adj_angle)
    
    # prop_adj_angle <- max(weight_angle_thresh - abs(mean(c(adj_angle, mean(angles_std[adj_verts]) %% (2*pi))) - angle_circle), 0) / weight_angle_thresh / 2
    prop_adj_angle <- (weight_angle_thresh - closestAngle(angle_circle, angles_std[closest_vertex])$distance) / (2*weight_angle_thresh)
    prop_adj_angle <- smooth0to0.5(max(prop_adj_angle,0))
    
    
    # angle_out <- mean(angles_std[adj_verts]) * (1 - prop_adj_angle) + adj_angle * prop_adj_angle
    
    curr_face_angle <- angles_std[closest_vertex]
    curr_face_angle <- curr_face_angle - closestAngle(focal = angle_circle, thetas = angles_std)$direction * pi/n_vert
    curr_face_angle <- coerceTo2Pi(curr_face_angle)
    if(abs(adj_angle - curr_face_angle) > pi){
      if(curr_face_angle > adj_angle){
        curr_face_angle <- curr_face_angle - 2*pi
      } else {
        adj_angle <- adj_angle - 2*pi
      }
    }
    
    angle_out <- curr_face_angle * (1 - prop_adj_angle) + adj_angle * prop_adj_angle
    
    # if(sign(diff(angles_std[adj_verts])) == 1){
    #   # print(paste0(diff(angles_std[adj_verts]), collapse = " "))
    #   angle_out <- mean(c(angles_std[adj_verts[1]]-2*pi, angles_std[adj_verts[2]])) * (1 - prop_adj_angle) + adj_angle * prop_adj_angle
    # } else {
    #   angle_out <- mean(angles_std[adj_verts]) * (1 - prop_adj_angle) + adj_angle * prop_adj_angle  
    # }
    
    
    
    # cent_circ <- polar2cart(angle_out - ifelse(sign(diff(angles_std[adj_verts])) == -1, pi, 0), r_circle) + touch_pt
    cent_circ <- polar2cart(angle_out, r_circle) + touch_pt
    
    r_circle_modified <- r_circle
    
    if(n_worm_background_circles_to_draw[j,i,f1] > 0){
      nbgc <- n_worm_background_circles_to_draw[j,i,f1]
      r_circle_modified <- r_circle + nbgc / max_nbgc * r_circle * max_bgc_radius_prop_main
      for(bgc in 1:nbgc){
        # print(paste0(j, ", ", i, ", ", round(r_circle + r_circle * bgc / max_nbgc * max_bgc_radius_prop, 2), ", ", round(nbgc / max_nbgc * (nbgc - bgc) / nbgc * max_nbgc_brightness, 3)))
        draw_circle(cent_circ[1], cent_circ[2], r_circle_modified + r_circle_modified * bgc / max_nbgc * max_bgc_radius_prop, 
                    col = adjustcolor(cols_circle[1], nbgc / max_nbgc * (nbgc - bgc) / nbgc * max_nbgc_brightness), lwd = 0.1, border = NA)
      }
    }
    
    
    
    
    draw_circle(cent_circ[1], cent_circ[2], r_circle_modified, col = cols_circle[i], lwd = 1, border = cols[1])
    # points(touch_pt[1], touch_pt[2])
    points(cent_circ[1], cent_circ[2], pch = 19, col = cols[1], cex = 2)
    if(i > 1){
      segments(x0 = cent_circ[1], y0 = cent_circ[2], x1 = prev_cent_circ[1], y1 = prev_cent_circ[2], col = cols[1], lwd = 3)
    }
    
    prev_cent_circ <- cent_circ
    
    if(debug_animation){
      text(x = cent_circ[1], y = cent_circ[2], labels = j, col = "white", cex = 2)
    }
    # print(paste0(i, "-- angle_circ:", angle_circle, ", ang_out:", angle_out, ", touch_r:", touch_r,
    #              ", prop_adj_ang:", prop_adj_angle, ", adj_ang:", adj_angle, ", closest_vert_ang:", angles[closest_vertex]))
    
  }
}



# #inscribe smaller polygon inside
for(i in 2:(n_inscr_polygons)){
  n_vert <- n_vert
  side_angle <- (pi - 2*pi/n_vert) / 2
  third_angle <- pi - side_angle - pi/n_vert
  radius_reduction_ratio <- sin(side_angle) / sin(third_angle)
  # init_angle <- mean(c(coerceTo2Pi(angles[1]), coerceTo2Pi(angles[2])))
  init_angle <- starting_angle + pi/n_vert * (i-1)
  init_angle <- coerceTo2Pi(init_angle + polygon_speeds[i] * f1 * ifelse(i %% 2 == 0, -1, 1))
  # if(i == 2){print(init_angle)}
  # r <- cos(abs(angles[1] - init_angle)) * r
  r <- init_r * radius_reduction_ratio^(i-1)
  angles <- seq(init_angle, init_angle + 2*pi, length.out = n_vert + 1)
  verts <- t(sapply(angles, function(a) polar2cart(a, r)))
  # verts <- t(sapply(1:nrow(verts), function(v) rot(verts[v,1], verts[v,2], polygon_speeds[i] * f1 * ifelse(i %% 2 == 0, -2, 1))))
  # polygon(verts[,1], verts[,2], col = cols[i], border = cols[n_inscr_polygons - i], lwd = 0.1)
  polygon(verts[,1], verts[,2], col = adjustcolor(cols[i], alpha.f = polygon_opacities[i]*max_opacity_val[f1]), border = ifelse(i != n_inscr_polygons, NA, cols[1]), lwd = 3)
  
  
  #find layer down's angles
  if(i < n_inscr_polygons){
    init_angle_next <- starting_angle + pi/n_vert * (i)
    init_angle_next <- coerceTo2Pi(init_angle_next + polygon_speeds[i+1] * f1 * ifelse((i+1) %% 2 == 0, -1, 1))
    r_next <- init_r * radius_reduction_ratio^(i)
    angles_next <- seq(init_angle_next, init_angle_next + 2*pi, length.out = n_vert + 1)
    verts_next <- t(sapply(angles_next, function(a) polar2cart(a, r_next)))
    # text(labels = i+1, x = verts_next[1,1], y = verts_next[1,2], col = "white")
  }
  
  #draw lines connecting vertices to center
  if(i != n_inscr_polygons){
    for(j in 1:(nrow(verts)-1)){
      conduit_color <- conduit_cols_array[i,j,f1]
      # text(labels = j, x = verts[j,1], y = verts[j,2], col = "white")
      if(debug_animation){
        if(T){text(labels = paste0(i, ", ", j,  ": ", round(coerceTo2Pi(cart2polar(x = verts[j,1], y = verts[j,2])[1]), 3)), x = verts[j,1], y = verts[j,2], col = "white")}
      }
      
      #find starting location
      # starting_loc <- c(0,0)
      layer_down_closest_angle <- closestAngle(angles[j], angles_next)
      layer_down_closest_vertices <- c(layer_down_closest_angle$index, (layer_down_closest_angle$index - layer_down_closest_angle$direction) %% n_vert)
      layer_down_closest_vertices[layer_down_closest_vertices == 0] <- n_vert
      starting_loc <- intersecting_point(xy0a = verts_next[layer_down_closest_vertices[1], ], xy0b = verts_next[layer_down_closest_vertices[2], ], 
                                         xy1a = c(0,0), xy1b = verts[j,])
      
      # if(j == 1){
      #   text(labels = j, x = verts[j,1], y = verts[j,2], col = "white", cex = 3)
      #   points(verts_next[layer_down_closest_vertices[1], 1],verts_next[layer_down_closest_vertices[1], 2], cex = 3, col = "red")
      #   points(verts_next[layer_down_closest_vertices[2], 1], verts_next[layer_down_closest_vertices[2], 2], cex = 3, col = "red")
      # }
      # points(starting_loc[1], starting_loc[2])
      
      
      segments(x0 = starting_loc[1], y0 = starting_loc[2], x1 = verts[j,1] * 0.995, y1 = verts[j,2] * 0.995, col = conduit_color, lwd = 3)
      if(conduit_color != conduit_cols[1]){
        for(fadeout in 1:4){
          segments(x0 = starting_loc[1], y0 = starting_loc[2], x1 = verts[j,1] * 0.99, y1 = verts[j,2] * 0.99, lwd = 3 + fadeout * 2.5, 
                   col = adjustcolor(conduit_color, alpha.f = 1 / (fadeout + 5) + (1 - polygon_opacities[i]*max_opacity_val[f1])))
        }
      }
    }
  }
}

#put signature in 
text("@markov.chain.monte.carlo", x = 1.675, y = -2.1, col = cols[2], xpd = NA, family = "Helvetica", cex = 1.5)

dev.off()
}

n_runs <- n_runs + 1
if(n_runs > n_runs_reset){
  stopCluster(cl = cl)
  rm(cl)
}

#render video w/ filters
base_filename <- paste0("worms_", pwidth, "x" , pheight)
raw_filename <- paste0("raw_", base_filename, ".mp4")
blur_filename <- paste0("blur_", base_filename, ".mp4")
final_filename <- paste0(base_filename, ".mp4")


if(render_video){
  if(file.exists(paste0("~/Documents/geometric_animation/", raw_filename))){file.remove(paste0("~/Documents/geometric_animation/", raw_filename))}
  if(file.exists(paste0("~/Documents/geometric_animation/", blur_filename))){file.remove(paste0("~/Documents/geometric_animation/", blur_filename))}
  if(file.exists(paste0("~/Documents/geometric_animation/", final_filename))){file.remove(paste0("~/Documents/geometric_animation/", final_filename))}
  
  system(paste0("cd Documents/geometric_animation; ffmpeg -r ", nfps / framethin," -f image2 -s ", pwidth, "x" , pheight," -i frames/frame_%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ", raw_filename))
  # system(paste0("cd Documents/geometric_animation; ffmpeg -i ", raw_filename," -vf gblur=sigma=2:steps=6 -pix_fmt yuv420p ", blur_filename))
  # system(paste0("cd Documents/geometric_animation; ffmpeg -i ", raw_filename," -i ", blur_filename, " -filter_complex \"blend=lighten\" ", final_filename))
}
