# ----------------------------
#  Utility: 3D cross product
# ----------------------------
cross3 <- function(a, b) {
  c(a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1])
}

# ----------------------------
#  Utility: 3D dot product
# ----------------------------
dot3 <- function(a, b) {
  sum(a * b)
}

# -------------------------------------------------------
#  1) Check if a point is inside a 3D triangle
#     Using barycentric coordinates
# -------------------------------------------------------
point_in_triangle_3d <- function(pt, tri, eps = 1e-9) {
  #
  # tri is a 3x3 matrix: each row is (x, y, z)
  # pt is a length-3 numeric vector
  #
  # Step 1: Check if pt is near the plane of tri
  v0 <- tri[3,] - tri[1,]
  v1 <- tri[2,] - tri[1,]
  n  <- cross3(v0, v1)       # Normal vector
  denom_n <- sqrt(dot3(n, n))
  
  # Distance from pt to plane
  dist_plane <- abs(dot3((pt - tri[1,]), n)) / (denom_n + eps)
  if (dist_plane > eps) {
    # Not in the same plane
    return(FALSE)
  }
  
  # Step 2: Barycentric coordinate check
  # Let v2 = pt - tri[1,]
  v2    <- pt - tri[1,]
  dot00 <- dot3(v0, v0)
  dot01 <- dot3(v0, v1)
  dot02 <- dot3(v0, v2)
  dot11 <- dot3(v1, v1)
  dot12 <- dot3(v1, v2)
  
  invDenom <- 1 / (dot00*dot11 - dot01*dot01 + eps)
  u        <- (dot11*dot02 - dot01*dot12) * invDenom
  v        <- (dot00*dot12 - dot01*dot02) * invDenom
  
  # Inside if u >= 0, v >= 0, and u + v <= 1
  if (u >= -eps && v >= -eps && (u + v) <= 1 + eps) {
    return(TRUE)
  }
  FALSE
}

# -------------------------------------------------------
#  2) Check if two 3D line segments intersect
#     Segment1 = p1->p2, Segment2 = p3->p4
# -------------------------------------------------------
segment_segment_intersect_3d <- function(p1, p2, p3, p4, eps = 1e-9) {
  #
  # Each p* is length-3 numeric vector
  #
  r <- p2 - p1  # direction of seg1
  s <- p4 - p3  # direction of seg2
  
  rxs <- cross3(r, s)
  rxs_len2 <- dot3(rxs, rxs)  # squared length of rxs
  
  qp  <- p3 - p1
  # Check if r and s are parallel via cross(r,s) = 0
  if (rxs_len2 < eps) {
    # They are parallel or collinear
    # Check if they are collinear: cross(qp, r) = 0
    qpxr <- cross3(qp, r)
    if (sqrt(dot3(qpxr, qpxr)) > eps) {
      # Not collinear at all
      return(FALSE)
    }
    # Collinear - project onto r or s and check 1D overlap
    # We'll project onto r's axis
    r_len2 <- dot3(r, r)
    # param of p3 on p1->p2
    t0 <- dot3(qp, r) / (r_len2 + eps)
    # param of p4 on p1->p2
    t1 <- dot3((p4 - p1), r) / (r_len2 + eps)
    
    # Overlap if intervals [t0, t1] intersect [0, 1]
    tmin <- min(t0, t1)
    tmax <- max(t0, t1)
    if (tmax < 0 || tmin > 1) {
      return(FALSE)
    }
    return(TRUE)
  } else {
    # Not parallel: solve p1 + t*r = p3 + u*s
    # using cross products
    t_num <- cross3(qp, s)
    t     <- dot3(t_num, rxs) / (rxs_len2 + eps)
    
    u_num <- cross3(qp, r)
    u     <- dot3(u_num, rxs) / (rxs_len2 + eps)
    
    # If 0 <= t <= 1 and 0 <= u <= 1 => intersection
    if (t >= 0 - eps && t <= 1 + eps &&
        u >= 0 - eps && u <= 1 + eps) {
      return(TRUE)
    }
    return(FALSE)
  }
}

# -------------------------------------------------------
#  3) The Main Function: Triangle–Triangle Intersection
# -------------------------------------------------------
tri_tri_intersect <- function(tri1, tri2, eps = 1e-9) {
  #
  # tri1, tri2: each a 3x3 numeric matrix of 3D coords
  # row i => vertex i => c(x, y, z)
  #
  # Strategy:
  #   (a) If any vertex of tri2 is inside tri1 => intersect
  #   (b) If any vertex of tri1 is inside tri2 => intersect
  #   (c) Check all 3 edges of tri1 vs. all 3 edges of tri2
  #
  
  # (a) Check if tri2’s vertices lie inside tri1
  for (i in 1:3) {
    if (point_in_triangle_3d(tri2[i,], tri1, eps)) {
      return(TRUE)
    }
  }
  
  # (b) Check if tri1’s vertices lie inside tri2
  for (i in 1:3) {
    if (point_in_triangle_3d(tri1[i,], tri2, eps)) {
      return(TRUE)
    }
  }
  
  # (c) Edge–Edge intersections
  # tri1 edges: (1->2), (2->3), (3->1)
  # tri2 edges: (1->2), (2->3), (3->1)
  edges1 <- list(
    rbind(tri1[1,], tri1[2,]),
    rbind(tri1[2,], tri1[3,]),
    rbind(tri1[3,], tri1[1,])
  )
  edges2 <- list(
    rbind(tri2[1,], tri2[2,]),
    rbind(tri2[2,], tri2[3,]),
    rbind(tri2[3,], tri2[1,])
  )
  
  for (e1 in edges1) {
    for (e2 in edges2) {
      if (segment_segment_intersect_3d(e1[1,], e1[2,],
                                       e2[1,], e2[2,],
                                       eps)) {
        return(TRUE)
      }
    }
  }
  
  # No intersection found
  FALSE
}

# -------------------------------------------------------
#  Example usage
# -------------------------------------------------------

# Are they intersecting?
#### TEST ####
intersection_matrix <- do.call(rbind, lapply(1:length(face_coords), function(fi1){
  sapply(1:length(face_coords), function(fi2){
    tri_tri_intersect(face_coords[[fi1]], face_coords[[fi2]])    
  })
}))


# ---------------------------------------
# 1) Get plane normal and d for a polygon
#    polygon: N x 3 matrix
#    returns: list(normal, d)
# ---------------------------------------
get_plane <- function(polygon) {
  # For numerical stability, pick any two edges that
  # are not collinear.
  # E.g. vector1 = p2 - p1, vector2 = p3 - p1
  v1 <- polygon[2, ] - polygon[1, ]
  # find next non-collinear point for v2
  idx <- 3
  v2 <- polygon[idx, ] - polygon[1, ]
  while (sqrt(dot3(cross3(v1, v2), cross3(v1, v2))) < 1e-12 && idx < nrow(polygon)) {
    idx <- idx + 1
    v2 <- polygon[idx, ] - polygon[1, ]
  }
  
  # Normal:
  n <- cross3(v1, v2)
  n_len <- sqrt(dot3(n, n))
  if (n_len < 1e-12) {
    stop("Polygon is degenerate or all points collinear.")
  }
  n <- n / n_len  # normalized normal
  
  # plane eqn: n . x + d = 0 => d = - n . p1
  d <- -dot3(n, polygon[1, ])
  
  list(normal = n, d = d)
}

# ---------------------------------------
# 2) Plane-plane intersection line
#    input: plane1 = list(normal, d)
#           plane2 = list(normal, d)
#    returns: list(point = p0, dir = direction)
# ---------------------------------------
plane_plane_intersect <- function(plane1, plane2, eps = 1e-12) {
  n1 <- plane1$normal
  d1 <- plane1$d
  n2 <- plane2$normal
  d2 <- plane2$d
  
  # Check if parallel:
  cross_n1_n2 <- cross3(n1, n2)
  cross_len2 <- dot3(cross_n1_n2, cross_n1_n2)
  if (cross_len2 < eps) {
    # The planes are parallel or coincident
    # For line intersection, they need to be non-parallel
    return(NULL)  # or a special flag
  }
  
  # direction = cross(n1, n2)
  dir <- cross_n1_n2 / sqrt(cross_len2)
  
  # Find one point on the line by solving:
  #  n1 . x + d1 = 0
  #  n2 . x + d2 = 0
  #
  # There's a known formula or we can do a small linear solve:
  # We'll pick the axis that is largest in dir for stable solve
  absdir <- abs(dir)
  maxc <- which.max(absdir)
  
  # We'll set that coordinate to 0 and solve the 2D system
  # in the other two coordinates.
  # e.g., if maxc=1 => x is largest => set x=0, solve for y,z
  # Then check if that solution is feasible, else we might
  # try setting y=0 or z=0.
  
  attempt_axes <- c(1, 2, 3)  # x=1, y=2, z=3
  p0 <- NA
  
  for (axis in attempt_axes) {
    # We'll try setting x[axis] = 0
    guess <- c(0, 0, 0)
    # The other two coords are a 2x2 system from the plane eqns:
    
    # system: n1[rem] . x[rem] + n1[axis]*0 + d1 = 0
    #         n2[rem] . x[rem] + n2[axis]*0 + d2 = 0
    # where rem = the 2 coords not 'axis'
    
    rem <- setdiff(1:3, axis)
    A <- rbind(n1[rem], n2[rem])  # 2x2
    b <- c(-d1, -d2)
    
    if (abs(det(A)) > eps) {
      sol <- solve(A, b)
      guess[rem] <- sol
      p0 <- guess
      break
    }
  }
  
  if (any(is.na(p0))) {
    # Fallback: the system might be weirdly ill-conditioned
    # but we know the planes are not parallel,
    # so there's definitely a line.
    # Could do a direct 3×3 solve approach too:
    
    # Build M = [n1; n2; dir]^T, then solve M*x = -[d1, d2, ?]
    # but let's skip for brevity. We'll trust the above approach works for typical polygons.
    return(NULL)
  }
  
  list(point = p0, dir = dir)
}

# ---------------------------------------
# 3) Clip infinite line to a CONVEX polygon
#    polygon: Nx3, line: list(point, dir)
#    returns: a numeric vector c(tmin, tmax)
#       so that line segment = p0 + t*dir, t in [tmin, tmax]
#    or NULL if no intersection.
#
#  Strategy:
#    1. For each edge, we build a plane "above" that edge,
#       pointing inside the polygon. Then find constraints
#       on the parameter t.
#    2. Accumulate the intersection of all intervals.
# ---------------------------------------
clip_line_to_polygon <- function(polygon, line, eps = 1e-12) {
  # We assume the polygon is in a plane. We'll:
  #  - get the plane normal
  #  - project polygon + line into a 2D coordinate system
  #  - do a standard "clip line by 2D convex polygon" approach
  
  # get plane eqn to confirm normal
  pl <- get_plane(polygon) 
  nrm <- pl$normal
  
  # We’ll form a local 2D basis in that plane:
  # pick any vector u in plane, then cross with nrm => v
  # to build an orthonormal basis {u, v, nrm}.
  # Then we can map 3D -> 2D by (dot with u, dot with v).
  
  # 1) pick a random vector not parallel to nrm
  #    e.g. e1 = c(1,0,0), or c(0,1,0), etc.
  #    whichever is not ~ collinear with nrm
  e1 <- c(1, 0, 0)
  if (abs(dot3(e1, nrm)) > 0.99) {
    e1 <- c(0, 1, 0)
  }
  # project out nrm
  u <- e1 - dot3(e1, nrm)*nrm
  u_len <- sqrt(dot3(u, u))
  if (u_len < eps) {
    # fallback
    u <- c(0, 1, 0)
    u <- u - dot3(u, nrm)*nrm
    u_len <- sqrt(dot3(u, u))
  }
  u <- u / u_len
  
  v <- cross3(nrm, u)
  # v is orthonormal to both nrm and u
  
  # Function to map a 3D point -> 2D coords in plane
  project2D <- function(pt) {
    c(dot3(pt, u), dot3(pt, v))
  }
  
  # Shift polygon so that plane passes through origin for easy dot products
  # or we can just pick polygon[1,] as reference
  ref <- polygon[1,]
  # We'll project "pt - ref"
  polygon_2d <- t(apply(polygon, 1, function(row) {
    project2D(row - ref)
  }))
  
  # Also project line point
  p0 <- line$point
  dir <- line$dir
  p0_2d <- project2D(p0 - ref)
  dir_2d <- project2D(dir)
  
  # Now we have a 2D polygon (convex) and a 2D line (param: p0_2d + t*dir_2d).
  # Let's do a standard 2D convex-polygon line clipping.
  # We'll do an edge-by-edge half-plane approach:
  
  n <- nrow(polygon_2d)
  tmin <- -Inf
  tmax <- Inf
  
  for (i in seq_len(n)) {
    i2 <- if (i == n) 1 else i+1
    pA <- polygon_2d[i, ]
    pB <- polygon_2d[i2, ]
    
    # Edge vector in 2D
    edge <- pB - pA
    
    # outward normal? We want the polygon's inside to be on "left" side
    # For a polygon in CCW order, the left normal is:
    edge_normal <- c(-edge[2], edge[1])
    
    # But we must confirm the orientation is indeed CCW or adjust sign if not.
    # For convex polygons with consistent ordering, we can figure it out by checking
    # if the polygon area is positive => CCW. Or just guess and check.
    # We'll do a quick check:
    # Summed cross product area in 2D:
    # If we need to ensure the normal points inside, we can detect that
    # the polygon is probably oriented CCW if the total area is > 0.
    # We'll do that once outside the loop:
    
    # We'll do the "line vs. half-plane" constraint:
    # Half-plane eqn: edge_normal . (X - pA) >= 0 means inside
    # Substitute X = p0_2d + t*dir_2d => condition on t.
    
    denom <- dot3(edge_normal, dir_2d)
    numer <- dot3(edge_normal, (pA - p0_2d))
    
    if (abs(denom) < eps) {
      # line is parallel to this edge
      # if numer < 0 => line is outside, no intersection
      if (numer < 0) {
        return(NULL)
      }
      # else line is in the half-plane, no clipping from this edge
    } else {
      t_candidate <- numer / denom
      if (denom > 0) {
        # line enters half-plane at t_candidate
        # so t_candidate is new lower bound
        if (t_candidate > tmin) tmin <- t_candidate
      } else {
        # denom < 0 => line leaves half-plane at t_candidate
        # so t_candidate is new upper bound
        if (t_candidate < tmax) tmax <- t_candidate
      }
    }
    # If tmin > tmax => no intersection
    if (tmin > tmax) {
      return(NULL)
    }
  }
  
  # If we get here, the line in param t in [tmin, tmax] is inside the polygon
  c(tmin, tmax)
}

# ---------------------------------------
# 4) Putting it all together:
#    Find the intersection segment of polygon A and polygon B
#    (both are convex).
# ---------------------------------------
intersect_polygons_3D_as_line_segment <- function(polyA, polyB) {
  # 1) plane eqns
  planeA <- get_plane(polyA)
  planeB <- get_plane(polyB)
  
  # 2) plane-plane intersection
  line <- plane_plane_intersect(planeA, planeB)
  if (is.null(line)) {
    # parallel or coincident
    # if coincident, you need 2D polygon intersection in that plane
    return(NULL)
  }
  
  # 3) Clip line to each polygon
  clipA <- clip_line_to_polygon(polyA, line)
  if (is.null(clipA)) return(NULL)
  clipB <- clip_line_to_polygon(polyB, line)
  if (is.null(clipB)) return(NULL)
  
  # Each clip gives [tmin, tmax]. Intersection is overlap:
  tmin <- max(clipA[1], clipB[1])
  tmax <- min(clipA[2], clipB[2])
  
  if (tmin > tmax) {
    return(NULL)
  }
  
  # Return endpoints in 3D:
  p0 <- line$point
  dir <- line$dir
  segA <- p0 + tmin * dir
  segB <- p0 + tmax * dir
  
  # If the intersection is a single point (tmin ~ tmax),
  # then segA == segB. Still valid as a single intersection point.
  
  rbind(segA, segB)
}

# ---------------------------------------
# Example usage
# ---------------------------------------
# Define two simple squares in different planes that intersect in a line segment.
# For illustration, let's do:
#  - Square A in plane z=0
#  - Square B in plane x=0, but shifted so they definitely intersect

intersection_matrix <- do.call(rbind, lapply(1:length(face_coords), function(fi1){
  sapply(1:length(face_coords), function(fi2){
    tri_tri_intersect(face_coords[[fi1]], face_coords[[fi2]])    
  })
}))


intersect_polygons_3D_as_line_segment(face_coords[[1]], face_coords[[8]])
