library(rgl)

options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)

triangulate_convex_fan <- function(P) {
  # Triangle fan: (1, i, i+1), i = 2..n-1
  # Assumes vertices are ordered around the polygon.
  P <- as.matrix(P)
  n <- nrow(P)
  if (n < 3) stop("polygon needs at least 3 points")
  idx <- cbind(1, 2:(n-1), 3:n)
  # reorder into a 3*(n-2) x 3 coordinate matrix for triangles3d
  V1 <- P[idx[,1], , drop = FALSE]
  V2 <- P[idx[,2], , drop = FALSE]
  V3 <- P[idx[,3], , drop = FALSE]
  rbind(V1, V2, V3)
}


is_coplanar <- function(P, tol = 1e-6) {
  # P: n x 3 matrix/data.frame
  P <- as.matrix(P)
  n <- nrow(P)
  if (n < 3) return(TRUE)
  # find first set of non-collinear points
  i <- 1; j <- 2; k <- 3
  v1 <- P[j, ] - P[i, ]
  v2 <- P[k, ] - P[i, ]
  tries <- 0
  while (tries < n && sqrt(sum((v1 / (sqrt(sum(v1^2)) + 1e-12))^2)) > 0 &&
         sqrt(sum((v2 / (sqrt(sum(v2^2)) + 1e-12))^2)) > 0 &&
         sqrt(sum(crossprod(crossprod(matrix(c(v1, v2), ncol = 3, byrow = TRUE))))) < tol) {
    # extremely defensive, but we'll just break out below; keep simple:
    break
  }
  # normal from first 3 points
  nrm <- pracma_cross(v1, v2)
  if (sum(abs(nrm)) < tol) return(TRUE) # degenerate -> treat as coplanar
  nrm <- nrm / sqrt(sum(nrm^2))
  d <- -sum(nrm * P[i, ])
  dist <- abs(P %*% nrm + d)
  all(dist < tol)
}

# small cross product without extra deps
pracma_cross <- function(a, b) {
  c(a[2]*b[3]-a[3]*b[2],
    a[3]*b[1]-a[1]*b[3],
    a[1]*b[2]-a[2]*b[1])
}

triangulate_convex_fan(square_xy)

plot_polygons_rgl <- function(polys, colors = NULL, alpha = 0.7, wire = TRUE) {
  if (!is.list(polys)) stop("'polys' must be a list of n x 3 matrices/data.frames")
  if (is.null(colors)) colors <- rep_len(rainbow(length(polys)), length(polys))
  if (length(colors) != length(polys)) stop("length(colors) must match length(polys)")
  
  open3d()                      # opens a headless rgl device
  bg3d("white")
  
  for (i in seq_along(polys)) {
    P <- as.matrix(polys[[i]])
    if (ncol(P) != 3L || nrow(P) < 3L) next
    
    # simple triangle-fan fill assuming convex & ordered
    if (!is_coplanar(P)) {
      cat(sprintf("[debug] polygon %d not coplanar within tolerance; drawing outline only\n", i))
      if (wire) lines3d(P_closed, col = colors[i], lwd = 2)
      next
    }
    if (nrow(P) == 3) {
      # triangle
      tris <- P
      triangles3d(tris, col = colors[i], alpha = alpha)
    } else if (nrow(P) == 4) {
      # triangle
      quads <- P
      quads3d(quads, col = colors[i], alpha = alpha)
    } else {
      # triangulate as a convex fan
      tris <- triangulate_convex_fan(P)
      triangles3d(tris, col = colors[i], alpha = alpha)
      cat(sprintf("[debug] polygon %d triangulated into %d triangles\n", i, nrow(tris)/3))
    }
    
    if (wire) lines3d(rbind(P, P[1,]), col = "black", lwd = 1.5)
  }
  
  axes3d(); aspect3d(1,1,1)
  rglwidget()                   # return the widget (critical for your build)
}

# demo polygons
square_xy <- rbind(c(0,0,0), c(1,0,0), c(1,1,0), c(0,1,0)) + 0.2
tri_tilt  <- rbind(c(0,0,0), c(1,0.2,0.8), c(0.2,1.0,0.5))
pent      <- rbind(c(0,0,0), c(1.2,0.1,0.3), c(1.4,1.1,0.6), 
                   c(0.3,1.4,0.4), c(-0.2,0.7,0.2))

# plot_polygons_rgl(list(square_xy, tri_tilt, pent),
#                   colors = c("#1f77b4", "#2ca02c", "#d62728"), alpha = 0.6, wire = TRUE)
plot_polygons_rgl(list(square_xy, tri_tilt),
                  colors = c("#1f77b4", "#2ca02c"), alpha = 0.6, wire = T)
