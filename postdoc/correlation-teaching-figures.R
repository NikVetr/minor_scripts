library(rgl)

n <- 1000
x1 <- seq(-1,1,length.out=n/2)
y1 <- sqrt(1-x1^2)
x2 <- rev(x1)[-1]
y2 <- -sqrt(1-x2^2)
x <- c(x1,x2)
y <- c(y1,y2)
par(pty="s")

#unit circle, with points added
plot(x,y, type = "l")
points(0,0,pch=19,cex=1)
abline(h=0, lty=3, col = adjustcolor(1,0.4))
abline(v=0, lty=3, col = adjustcolor(1,0.4))
points(sqrt(1/2), -sqrt(1/2), pch = 19, col = 2, cex = 2)
text(sqrt(1/2) + 0.15, -sqrt(1/2) - 0.15, 
     labels = expression(bgroup("{", sqrt(1/2) ~ "," ~ sqrt(-1/2), "}")),
     cex = 0.55)
points(-sqrt(1/2), sqrt(1/2), pch = 19, col = 2, cex = 2)
text(-sqrt(1/2) - 0.15, sqrt(1/2) + 0.15, 
     labels = expression(bgroup("{", sqrt(1/2) ~ "," ~ sqrt(-1/2), "}")),
     cex = 0.55)
segments(0,0, sqrt(1/2), -sqrt(1/2),col=2,lty=2,lwd=2)
segments(0,0, -sqrt(1/2), sqrt(1/2),col=2,lty=2,lwd=2)
points(0,0,pch=19,cex=1)
arc <- cbind(x,y)[x>-sqrt(1/2) & y>-sqrt(1/2),] / 4
lines(arc, lty = 2)
text(0.3, 0.2, labels = expression(theta == pi))
text(0.5, -0.5, labels = expression(theta == 0), pos = 4)

plot(c(1,-1), y = c(-1,1), pch = 19, col = 2, cex = 2, 
     xlab = "variable 1", ylab = "variable 2")
abline(a = 0, b = -1, lty = 2, col = 2, lwd = 2)
plot(c(-1,1), y = c(-1,1), pch = 19, col = 2, cex = 2, 
     xlab = "variable 1", ylab = "variable 2")
abline(a = 0, b = 1, lty = 2, col = 2, lwd = 2)

# Function to generate a unit sphere mesh
plot_sphere <- function(n = 50) {
  u <- seq(0, pi, length.out = n)
  v <- seq(0, 2*pi, length.out = n)
  
  x <- outer(sin(u), cos(v))
  y <- outer(sin(u), sin(v))
  z <- outer(cos(u), rep(1, length(v)))
  
  persp3d(x, y, z, col = "lightblue", alpha = 0.3, add = TRUE)
}

# Function to plot the plane x + y + z = 0
plot_plane <- function() {
  x <- seq(-1, 1, length.out = 100)
  y <- seq(-1, 1, length.out = 100)
  xy_grid <- expand.grid(x = x, y = y)
  z <- matrix(-(xy_grid$x + xy_grid$y), nrow = length(x), ncol = length(y))
  x_mat <- matrix(xy_grid$x, nrow = length(x), ncol = length(y))
  y_mat <- matrix(xy_grid$y, nrow = length(x), ncol = length(y))
  bad_z <- which(!(z > -1 & z < 1), arr.ind = T)
  z[bad_z] <- NA
  surface3d(x_mat, y_mat, z, color = "orange", alpha = 0.3)
}

# Function to plot the great circle (intersection of the sphere and plane)
plot_great_circle <- function() {
  x1 <- seq(-sqrt(2/3), sqrt(2/3), length.out = 100)
  y1 <- 1/2*(sqrt(2-3*x1^2) - x1)
  x2 <- rev(x1)[-1]
  y2 <- 1/2*(-sqrt(2-3*x2^2) - x2)
  x <- c(x1, x2)
  y <- c(y1, y2)
  z <- -(x + y)
  lines3d(x, y, z, col = "red", lwd = 5)
}

# Function to retrieve coordinates for all points on the unit sphere that are
# theta radians away from a vector w
find_angle_circle_coords <- function(w, theta){
  coefs <- c(A=cos(theta)/w[3],
             B=w[1]/w[3],
             C=w[2]/w[3])
  yf <- function(x, coefs, p = T){
    A <- coefs["A"]
    B <- coefs["B"]
    C <- coefs["C"]
    a <- 1 + C^2
    b <- 2*B*C*x - 2*A*C
    c <- (1 + B^2) * x^2 - 2*A*B*x + A^2 - 1
    y <- (-b + ifelse(p, 1, -1) * sqrt(b^2 - 4*a*c)) / (2 * a)
    # y <- ((A*C - B*C*x) + ifelse(p, 1, -1) * sqrt(-A^2 + 2*A*B*x - B^2*x^2 - C^2*x^2 + C^2 - x^2 + 1)) /
    (C^2 + 1)
    names(y) <- NULL
    return(y)
  }
  # xroots2 <- function(coefs, eps = 1E-6) {
  #   A <- coefs["A"]
  #   B <- coefs["B"]
  #   C <- coefs["C"]
  #   
  #   # Compute coefficients for the quadratic discriminant
  #   a <- -4 * ((1 + C^2) * (1 + B^2) - B^2 * C^2)
  #   b <- 8 * ((1 + C^2) * A * B - A * B * C)
  #   c <- 4 * A^2 * C^2 - 4 * (1 + C^2) * (A^2 - 1)
  #   
  #   # Compute discriminant
  #   discriminant <- b^2 - 4 * a * c
  #   
  #   x1 <- (-b - sqrt(discriminant)) / (2 * a)
  #   x2 <- (-b + sqrt(discriminant)) / (2 * a)
  #   xb <- c(x1,x2)
  #   return(xb - sign(xb) * eps)
  # }
  
  xroots <- function(coefs, eps = 1E-6){ #uses the quadratic formula
    A <- coefs["A"]
    B <- coefs["B"]
    C <- coefs["C"]
    
    a <- -(B^2 + C^2 + 1)
    b <- 2 * A * B
    c <- C^2 - A^2 + 1
    
    delta <- b^2 - 4 * a * c
    
    x1 <- (-b + sqrt(delta)) / (2 * a)
    x2 <- (-b - sqrt(delta)) / (2 * a)
    xb <- c(x1,x2)
    
    
    return(sort(xb) + c(1,-1) * eps)
  }
  
  xb <- xroots(coefs)
  x1 <- seq(xb[1], xb[2], length.out = 100)
  x2 <- rev(x1)[-1]
  y1 <- yf(x1, coefs)
  y2 <- yf(x2, coefs, F)
  x <- c(x1, x2)
  y <- c(y1, y2)
  z <- coefs["A"] - coefs["B"] * x - coefs["C"] * y  
  return(data.frame(x=x, y=y, z=z))
}


# Initialize the 3D plot
par3d(windowRect = c(20, 30, 800, 800))
clear3d()
plot_sphere()
points3d(0,0,0, size = 5)
axes3d()
par3d(
  userMatrix = matrix(c(
    0.01495839,  0.7961579,  0.6049041, 0,
    0.58204389,  0.4849756, -0.6527042, 0,
    -0.81301934,  0.3618442, -0.4561449, 0,
    0.00000000,  0.0000000,  0.0000000, 1
  ), nrow = 4, byrow = TRUE),
  zoom = 0.8246962,
  FOV = 83.84414
)

#save view
index <- 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

#plot the plane, circle, and axes
plot_plane()

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

plot_great_circle()

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

#plot a line w on the unit circle (aka a unit vector)
set.seed(1)
w <- rnorm(3)
w <- (w - mean(w)) / norm(w - mean(w), "2")
lines3d(rbind(0,w), col = "green", lwd = 5)

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

#find cone at angle around that circle
theta <- pi/8
ic <- find_angle_circle_coords(w, theta)
ic <- rbind(ic, tail(ic, 1))

# draw the cone extending past by some amount
tri <- list(x = rbind(0, ic$x, c(ic$x[-1], ic$x[1])), 
            y = rbind(0, ic$y, c(ic$y[-1], ic$y[1])), 
            z = rbind(0, ic$z, c(ic$z[-1], ic$z[1])))
extend_by <- 1.0
triangles3d(tri$x * extend_by, 
            tri$y * extend_by, 
            tri$z * extend_by, 
            col = "blue", alpha = 0.3)

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

#draw the intersection
lines3d(ic, col = "blue", lwd = 5)

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

#plot intersection of cone and earlier circle
#can get algebraic solution later, for now just
poss_pts <- ic[order(abs(rowSums(ic)))[1:3],]
dupe_inds <- c(arrayInd(order(as.matrix(dist(poss_pts)))[4], .dim = c(3,3)))
ipts <- poss_pts[c(setdiff(1:3, dupe_inds), min(dupe_inds)),]
points3d(ipts, size = 20, col = "purple")

#draw lines to these points
lines3d(rbind(0, ipts[1,]), col = "hotpink", lwd = 10)
lines3d(rbind(0, ipts[2,]), col = "hotpink", lwd = 10)

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

#plot another vector on the circle
v <- rnorm(3)
v <- (v - mean(v)) / norm(v - mean(v), "2")
lines3d(rbind(0,v), col = "magenta", lwd = 5)

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

#draw arc between lines to represent angles
arcs <- list(t(v - t(t(t(0:50/50)) %*% t(t(v - ipts[1,])))),
             t(v - t(t(t(0:50/50)) %*% t(t(v - ipts[2,])))))
arcs[[1]] <- arcs[[1]] / apply(arcs[[1]], 1, norm, "2") * 0.5
arcs[[2]] <- arcs[[2]] / apply(arcs[[2]], 1, norm, "2") * 0.75
lines3d(arcs[[1]], col = "black", lwd = 4)
text3d(arcs[[1]][floor(nrow(arcs[[1]])/2.5),] * 1.2, 
       text = "θ₁", col = "black", cex = 1.25)

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))

lines3d(arcs[[2]], col = "black", lwd = 4)
text3d(arcs[[2]][floor(nrow(arcs[[1]])/4),] * 1.15, 
       text = "θ₂", col = "black", cex = 1.25)

#save view
index <- index + 1
rgl.snapshot(paste0("~/Documents/Correlation_Images/3D_", index, ".png"))
