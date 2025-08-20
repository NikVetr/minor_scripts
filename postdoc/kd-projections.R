# install.packages("pracma")  # if not yet installed
library(pracma)  # for nullspace()
library(expm)
library(doParallel)
library(foreach)
library(rgl)
library(grid)
library(png)

#source hypercube optimization script
source("~/scripts/minor_scripts/postdoc/collapse_hcube_efficient.R")
source("~/scripts/minor_scripts/postdoc/true_text_dim_functions.R")

# # Dimensions
# k <- 3  # ambient dimension
# p <- 2  # dimension of the subspace we want
# 
# # 1) Generate a random k x k matrix
# M <- matrix(rnorm(k * k), nrow = k)
# 
# # 2) Orthonormalize M via QR, then take the first p columns
# #    -> columns of B are p orthonormal vectors in R^3
# Q <- qr.Q(qr(M))       # Q is k x k orthonormal
# B <- Q[, 1:p, drop=FALSE]  # B is k x p = 3 x 2
# 
# cat("Basis matrix B (each column is an orthonormal vector in R^3):\n")
# print(B)
# cat("Check columns are orthonormal (B^T B):\n")
# print(t(B) %*% B)
# 
# # 3) The column space of B is our p-dimensional subspace in R^3
# #    The null space of B^T is its orthogonal complement
# N <- nullspace(t(B))  # B^T is p x k = 2 x 3, so the null space is 1-dimensional in R^3
# 
# cat("\nNull-space basis vector (orthogonal complement):\n")
# print(N)
# 
# # 4) Demonstrate orthogonality
# #    a) Pick a random combination of the basis columns -> a point in the subspace
# alpha <- rnorm(p)                   # random coefficients
# x_subspace <- B %*% alpha           # 3x1 vector in the subspace
# 
# #    b) Pick a random multiple of the null-space vector -> a point in orth. complement
# beta <- rnorm(1)
# x_null <- N %*% beta                # 3x1 vector in the orthogonal complement
# 
# #    c) Check their dot product
# dot_val <- sum(x_subspace * x_null) # Should be ~ 0 if truly orthogonal
# cat(sprintf("\nDot product of subspace vector and null-space vector: %g\n", dot_val))

#### functions ####
generate_hypercube_faces <- function(k, d) {
  stopifnot(d >= 1 && d <= k-1)
  
  vertices <- as.matrix(expand.grid(replicate(k, c(0, 1), simplify = FALSE))) - 0.5
  
  # Identify which dimensions to fix (k-d)
  fixed_dims_combinations <- combn(k, k - d, simplify = FALSE)
  
  face_list <- list()
  for (fixed_dims in fixed_dims_combinations) {
    # Each combination of fixed dimensions has 2^(k-d) ways to set values ±0.5
    fixed_vals_list <- expand.grid(replicate(k - d, c(-0.5, 0.5), simplify=FALSE))
    
    for (fixed_vals_idx in seq_len(nrow(fixed_vals_list))) {
      fixed_vals <- fixed_vals_list[fixed_vals_idx, ]
      
      # Select vertices matching fixed values exactly
      condition <- rep(TRUE, nrow(vertices))
      for (j in seq_along(fixed_dims)) {
        condition <- condition & (vertices[, fixed_dims[j]] == fixed_vals[[j]])
      }
      
      face_verts <- which(condition)
      face_list[[length(face_list) + 1]] <- sort(face_verts)
    }
  }
  
  face_matrix <- do.call(rbind, face_list)
  return(face_matrix)
}



generate_hypercube <- function(k, face_dim = NA) {
  
  if(is.na(face_dim)){
    face_dim <- k-1
  }
  
  # Generate all binary combos: 2^k x k
  binary_matrix <- as.matrix(
    expand.grid(replicate(k, c(0, 1), simplify = FALSE))
  )
  # Shift from [0,1] to [-0.5, 0.5]
  vertices <- binary_matrix - 0.5
  edge_matrix <- generate_hypercube_faces(k, 1)
  face_matrix <- generate_hypercube_faces(k, face_dim)
  
  return(list(vertices = vertices, 
              edges = edge_matrix,
              faces = face_matrix))
}


extract_subrotations <- function(rotmat) {
  k <- nrow(rotmat)
  result <- list()
  for (d in 2:(k - 1)) {
    V <- rotmat[, 1:d]
    # Project into the first d dimensions
    V_proj <- V[1:d, , drop = FALSE]
    # QR decomposition to orthonormalize projected vectors
    qr_decomp <- qr(V_proj)
    Q <- qr.Q(qr_decomp)
    
    for (j in 1:d) {
      if (sum(Q[,j] * V_proj[,j]) < 0) {
        Q[,j] <- -Q[,j]
      }
    }
    # Ensure right-handedness: check determinant
    if (det(Q) < 0) {
      Q[, d] <- -Q[, d]
    }
    R <- diag(k)
    R[1:d, 1:d] <- Q
    result[[paste0("R_", d)]] <- R
  }
  return(result)
}

pad_with_identity <- function(mat, k){
  out <- diag(k)
  out[1:nrow(mat), 1:ncol(mat)] <- mat
  out
}

d2rot <- function(theta) {matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)}

len <- function(x) {sqrt(sum(x^2))}

dkrot <- function(k, thetas) {
  planes <- t(combn(k, 2))  # same ordering as in `dkrot()`
  planes <- planes[order(apply(planes, 1, max)),]
  angles <- data.frame(cbind(planes, theta = thetas))
  rotmats <- lapply(angles$theta, d2rot)
  big_rotmats <- lapply(1:choose(k, 2), function(i){
    bigrm <- diag(k)
    bigrm[angles$V1[i], angles$V1[i]] <- rotmats[[i]][1,1]
    bigrm[angles$V1[i], angles$V2[i]] <- rotmats[[i]][1,2]
    bigrm[angles$V2[i], angles$V1[i]] <- rotmats[[i]][2,1]
    bigrm[angles$V2[i], angles$V2[i]] <- rotmats[[i]][2,2]
    bigrm
  })
  rotmat <- Reduce("%*%", big_rotmats)
  rotmat
}

rotation_matrix <- function(axis, angle_rad) {
  axis <- axis / sqrt(sum(axis^2))  # normalize
  K <- matrix(c(
    0, -axis[3], axis[2],
    axis[3], 0, -axis[1],
    -axis[2], axis[1], 0), nrow = 3, byrow = TRUE)
  R <- diag(3) + sin(angle_rad)*K + (1 - cos(angle_rad))*(K %*% K)
  return(R)
}

mcprint <- function(...){system(sprintf('printf "%s"', paste0(..., collapse="")))}

#### define unit hypercube ####
k <- 9
hcube <- generate_hypercube(k, face_dim = 2)
axis_bounds <- sqrt(k*0.5^2)
nrots <- 1440
nrots_per_full_rot <- 120

#now rotate it

#specify rotation
thetas <- runif(choose(k, 2), 3*pi/4, 3*pi/2)
rotmat <- dkrot(k, thetas)
poss_pt <- c(axis_bounds, rep(0, k-1))
pp_rot <- c(rotmat %*% poss_pt)
angle <- acos(sum(poss_pt * pp_rot) / axis_bounds^2)
prop_full_rot <- 1 / (2 * pi / angle)

#identify fractional rotation matrix
rm_log  <- logm(rotmat)
rm_factor <- nrots_per_full_rot * prop_full_rot
rm_frac <- expm(1/rm_factor * rm_log)
rm_fracs <- c(list(diag(k)), extract_subrotations(rm_frac), list(rm_frac))

#set view direction (ie, plane to project to)
M <- matrix(rnorm(k * k), nrow = k)
Q <- qr.Q(qr(M))
plane_vectors <- Q[, 1:2, drop=FALSE] 
plane_vectors <- cbind(c(1, rep(0,k-1)),
                       c(0, 1, rep(0,k-2)))
view_dir <- qr.Q(qr(plane_vectors), complete=TRUE)[, 3, drop = F]

#file paths
vid_dir <- "~/shape_rotation/"
frames_dir <- paste0(vid_dir, "frames/")
if(!dir.exists(frames_dir)){
  dir.create(frames_dir, recursive = T)  
} else {
  all(file.remove(list.files(frames_dir, full.names = T)))
}

#colors
face_alpha <- 1 / 1.25^k
# face_alpha <- 1
face_col <- "lightskyblue"
t0_face_col <- "lightskyblue"
face_cols <- adjustcolor(colorRampPalette(c(t0_face_col, face_col))(nrots), face_alpha)
face_cols <- c(face_cols, rev(face_cols))
# face_col <- "lightgrey"
polyborder_col <- adjustcolor(1, alpha.f = 1)
polyborder_col <- NA
seg_col <- "deepskyblue4"
t0_seg_col <- "deepskyblue4"
seg_cols <- colorRampPalette(c(t0_seg_col, seg_col))(nrots)
seg_cols <- c(seg_cols, rev(seg_cols))

#specify gradual expansion of the dimensions of rotation
dim_break_pts <- head(floor(seq(1, nrots, length.out = k+2)), k) #linear unfolding
dim_break_pts <- head(ceiling(nrots - seq(sqrt(nrots), 1, length.out = k+2)^2 + 1), k) #quadratic unfolding
pow_breaks <- 2
dim_break_pts <- head(ceiling(nrots - seq(nrots^(1/pow_breaks), 1, length.out = k+2)^pow_breaks + 1), k) #quadratic unfolding
dim_break_pts <- c(dim_break_pts, nrots+1)
dim_break_pts[1] <- 1
frac_1d2d_shift <- 0.75
frac_2d1d_shift <- 1 - frac_1d2d_shift
weight_ord <- (sin(-50:50/100*pi) + 1) / 2
n_dim_per_tim <- rep(1:k, diff(dim_break_pts))

#### compute forward rotation ####
rot_hcube <- hcube
rot_mats <- list()
cumrot_mats <- list()
full_hcubes <- list()
hcubes_proj <- list()
depths_proj <- list()
for(i in 1:nrots){
  
  #find which dims to rotate in
  ndim2r <- true_ndim <- sum(i >= dim_break_pts)
  if(ndim2r == 1){
    ndim2r <- 2
  }
  bordering_breaks <- dim_break_pts[sum(i >= dim_break_pts) + c(0,1)]
  trans_prop <- (i - bordering_breaks[1]) / diff(bordering_breaks)
  trans_ind <- floor(trans_prop * 100) + 1
  wnd <- weight_ord[trans_ind]
  
  #construct rotation matrix
  rm_to_use <- rm_fracs[[ndim2r]]
  
  #adjust for smooth transitions, blending current and previous rotation matrices
  # rm_to_use[ndim2r, ] <- rm_frac[ndim2r, ] * wnd + diag(k)[ndim2r, ] * (1-wnd)
  # rm_to_use[ndim2r, ndim2r] <- rm_frac[ndim2r, ndim2r]
  # rm_to_use[, ndim2r] <- rm_frac[, ndim2r] * wnd + diag(k)[, ndim2r] * (1-wnd)
  
  #compute rotation and project to the viewing plane
  rot_hcube$vertices <- t(rm_to_use %*% t(rot_hcube$vertices))
  full_hcubes[[i]] <- rot_hcube$vertices
  hcubes_proj[[i]] <- full_hcubes[[i]] %*% plane_vectors
  depths_proj[[i]] <- full_hcubes[[i]] %*% view_dir
  
  #record the rotation matrix that got us here
  rot_mats[[i]] <- rm_to_use
  if(i==1){
    cumrot_mats[[i]] <- rm_to_use
  } else {
    cumrot_mats[[i]] <- rm_to_use %*% cumrot_mats[[i-1]]  
  }
  
  #get a rotating line in there too
  if(true_ndim == 1){
    # square_verts <- cbind(hcubes_proj[[i]][duplicated(hcubes_proj[[i]]),], 0)
    square_verts <- cbind(hcubes_proj[[i]][tail(rot_hcube$faces, 1),], 0)
    v1 <- square_verts[1,]
    v2 <- square_verts[4,]
    rotation_axis <- v2 - v1
    rotation_axis <- rotation_axis / sqrt(sum(rotation_axis^2)) # normalize
    if(trans_prop > frac_1d2d_shift){
      sq_ra <- pi / 2  - 
        (trans_prop - frac_1d2d_shift) / (1 - frac_1d2d_shift) * pi / 2
    } else {
      sq_ra <- pi / 2  
    }
    sq_rm <- rotation_matrix(rotation_axis, sq_ra)
    rot_sq <- t(sq_rm %*% t(square_verts))[,1:2]
    hcubes_proj[[i]] <- matrix(c(t(rot_sq)), 
                               nrow = nrow(hcubes_proj[[i]]), 
                               ncol = ncol(hcubes_proj[[i]]), 
                               byrow = T)
  }
}

#### compute backward rotation ####
for(di in 1:(k-1)){
  
  #initialize target settings
  if(di == 1){
    
    nrots_curr <- nrots
    hcube_init <- tail(full_hcubes, 1)[[1]]
    R_init <- tail(rot_mats, 1)[[1]]
    n_steps <- tail(diff(dim_break_pts), 1)
    
  } else {
    
    nrots_curr <- length(full_hcubes)
    d2r <- (k-di+2):k
    hcube_init_full <- tail(full_hcubes, 1)[[1]][,-d2r]
    sub_verts <- !duplicated(round(tail(hcubes_proj, 1)[[1]], 10))
    hcube_init <- hcube_init_full[sub_verts,]
    subrots <- extract_subrotations(tail(rot_mats, 1)[[1]][1:(k - di + 2),
                                                           1:(k - di + 2)])
    R_init <- tail(subrots, 1)[[1]][-((k - di + 2)),-((k - di + 2))]
    n_steps <- tail(diff(dim_break_pts), di-1)[1]
    n_dim_per_tim <- c(n_dim_per_tim, rep(k - di + 1, n_steps))
    
  }
  
  #perform optimization and extract rotation matrices
  if((k-di) > 1){
    collapsed_hcube <- collapse_hcube(hcube_init, R_init, 
                                      n_steps = n_steps, shape_blend = 0.8,
                                      n_p = 2)
    rm2use <- collapsed_hcube$R_mats
  } else {
    rm2use <- lapply(1:n_steps, function(foo) R_init)
  }
  
  #pad these to have identity on outside
  if(di != 1){
    rm2use <- lapply(rm2use, function(x){
      new_rm <- diag(k)
      new_rm[1:dim(x)[1], 1:dim(x)[1]] <- x
      new_rm
    })
  }
  
  #compute projected coordinates
  for(i in (nrots_curr+1):(nrots_curr+n_steps)){
    #construct rotation matrix
    rm_to_use <- rm2use[[i-nrots_curr]]
    
    #compute rotation and project to the viewing plane
    rot_hcube$vertices <- t(rm_to_use %*% t(rot_hcube$vertices))
    full_hcubes[[i]] <- rot_hcube$vertices
    hcubes_proj[[i]] <- full_hcubes[[i]] %*% plane_vectors
    depths_proj[[i]] <- full_hcubes[[i]] %*% view_dir
    
    #record the rotation matrix that got us here
    rot_mats[[i]] <- rm_to_use
    cumrot_mats[[i]] <- rm_to_use %*% cumrot_mats[[i-1]]  
  }
}

#collapse square to line
nrots_curr <- length(full_hcubes)
n_steps <- diff(dim_break_pts[1:2])
n_dim_per_tim <- c(n_dim_per_tim, rep(2, n_steps * frac_2d1d_shift))
n_dim_per_tim <- c(n_dim_per_tim, rep(1, n_steps * c(1-frac_2d1d_shift)))
for(i in (nrots_curr+1):(nrots_curr+n_steps)){
  
  #update rot_hcube with constant rotation from last step
  rot_hcube$vertices <- t(rm_to_use %*% t(rot_hcube$vertices))
  
  #recover rotation axis
  hd_hcube <- rot_hcube$vertices[,1:2]
  square_verts <- cbind(hd_hcube[!duplicated(round(hd_hcube, 5)),], 0)
  v1 <- square_verts[1,]
  v2 <- square_verts[4,]
  rotation_axis <- v2 - v1
  rotation_axis <- rotation_axis / sqrt(sum(rotation_axis^2)) # normalize
  
  #identify position in line rotation
  trans_prop <- (i - nrots_curr) / n_steps
  
  #construct rotation angle
  if(trans_prop > frac_2d1d_shift){
    sq_ra <- pi / 2  
  } else {
    sq_ra <- pi / 2  - 
      (trans_prop - frac_1d2d_shift) / (1 - frac_1d2d_shift) * pi / 2
  }
  
  #rotate square into line
  sq_rm <- rotation_matrix(rotation_axis, sq_ra)
  rot_sq <- t(sq_rm %*% t(square_verts))[,1:2]
  
  #record variables
  hcubes_proj[[i]] <- matrix(c(t(rot_sq)), 
                             nrow = nrow(hcubes_proj[[i-1]]), 
                             ncol = ncol(hcubes_proj[[i-1]]), 
                             byrow = T)
  full_hcubes[[i]] <- rot_hcube$vertices
  depths_proj[[i]] <- full_hcubes[[i]] %*% view_dir
  rot_mats[[i]] <- rm_to_use
  cumrot_mats[[i]] <- rm_to_use %*% cumrot_mats[[i-1]]
}

#keep rotating line until we have returned to start

#compute rotation matrices
curr_line <- rot_sq[c(which.min(rot_sq[,1]), which.max(rot_sq[,1])),]
starting_rot_sq <- hcubes_proj[[1]][!duplicated(round(hcubes_proj[[1]], 5)),]
starting_line <- starting_rot_sq[c(which.min(starting_rot_sq[,1]), which.max(starting_rot_sq[,1])),]
curr_rot_mat <- rm_to_use[1:2,1:2]
curr_rot_angle <- -acos(curr_rot_mat[1,1])
starting_rot_mat <- rot_mats[[1]][1:2, 1:2]
starting_rot_angle <- -acos(starting_rot_mat[1,1])
rotation_needed <- find_min_directed_rotation(curr_line, starting_line, curr_rot_angle)
rotation_steps <- generate_smooth_rotation_steps(curr_rot_angle, starting_rot_angle, rotation_needed)

#### return to start ####
#rotate line to start
nrots_curr <- length(full_hcubes)
n_steps <- length(rotation_steps)
n_dim_per_tim <- c(n_dim_per_tim, rep(1, n_steps + 100))
for(i in (nrots_curr+1):(nrots_curr+n_steps)){
  
  #update rot_hcube with constant rotation from last step
  rotation_angle <- rotation_steps[i-nrots_curr]
  rm_to_use <- pad_with_identity(d2rot(rotation_angle), k)
  rot_hcube$vertices <- t(rm_to_use %*% t(rot_hcube$vertices))
  
  #recover rotation axis
  hd_hcube <- rot_hcube$vertices[,1:2]
  square_verts <- cbind(hd_hcube[!duplicated(round(hd_hcube, 5)),], 0)
  v1 <- square_verts[1,]
  v2 <- square_verts[4,]
  rotation_axis <- v2 - v1
  rotation_axis <- rotation_axis / sqrt(sum(rotation_axis^2)) # normalize
  
  #identify position in line rotation
  sq_ra <- pi / 2  
  
  #rotate square into line
  sq_rm <- rotation_matrix(rotation_axis, sq_ra)
  rot_sq <- t(sq_rm %*% t(square_verts))[,1:2]
  
  #record variables
  hcubes_proj[[i]] <- matrix(c(t(rot_sq)), 
                             nrow = nrow(hcubes_proj[[i-1]]), 
                             ncol = ncol(hcubes_proj[[i-1]]), 
                             byrow = T)
  full_hcubes[[i]] <- rot_hcube$vertices
  depths_proj[[i]] <- full_hcubes[[i]] %*% view_dir
  rot_mats[[i]] <- rm_to_use
  cumrot_mats[[i]] <- rm_to_use %*% cumrot_mats[[i-1]]
}


#set up parallel environment
use_multicore <- TRUE
if(use_multicore){
  if (!exists("cl")) {
    cl <- makeCluster(12)  # or number of cores
    registerDoParallel(cl)
  }
}

#find plotting ranges
max_val <- cummax(sapply(lapply(hcubes_proj, abs), max))

#get opacity variable
opa <- rep(1, length(hcubes_proj))

#generate PNGs for label with fancy R
hypercube_names <- c(
  "point (ℝ⁰)",
  "line segment (ℝ¹)",
  "square (ℝ²)",
  "cube (ℝ³)",
  "tesseract (ℝ⁴)",
  "penteract (ℝ⁵)",
  "hexeract (ℝ⁶)",
  "hepteract (ℝ⁷)",
  "octeract (ℝ⁸)",
  "ennearact (ℝ⁹)"
)
hypercube_names <- setNames(hypercube_names, sapply(strsplit(hypercube_names, " "), head, 1))
labels_dir <- paste0(vid_dir, "labels/")
if (!dir.exists(labels_dir)) dir.create(labels_dir, recursive = TRUE)
label_paths <- paste0(labels_dir, names(hypercube_names), ".png")

get_text_dims_px <- function(label, fontsize = 48, dpi = 96, fontfamily = "sans", padding_factor = 0.25) {
  tmp_file <- tempfile(fileext = ".png")
  png(tmp_file, width = 50, height = 20, units = "in", res = dpi, type = "cairo")
  on.exit({
    dev.off()
    unlink(tmp_file)
  })
  
  grid.newpage(recording = FALSE)
  gt <- grid::textGrob(label, gp = gpar(fontsize = fontsize, fontfamily = fontfamily))
  width_inches <- grid::convertWidth(grid::grobWidth(gt), "inches", valueOnly = TRUE)
  height_inches <- grid::convertHeight(grid::grobHeight(gt), "inches", valueOnly = TRUE)
  
  padded_width_inches <- width_inches * (1 + 2 * padding_factor)
  padded_height_inches <- height_inches * (1 + 2 * padding_factor) + (height_inches * 0.1) # Extra 10% height padding
  
  width_px <- ceiling(padded_width_inches * dpi)
  height_px <- ceiling(padded_height_inches * dpi)
  width_px <- max(width_px, 1)
  height_px <- max(height_px, 1)
  
  return(list(width = width_px, height = height_px))
}

render_label_png <- function(label, filename, fontsize = 48, dpi = 96, fontfamily = "sans", padding_factor = 0.25) {
  dims <- get_text_dims_px(label, fontsize = fontsize, dpi = dpi, fontfamily = fontfamily, padding_factor = padding_factor)
  png(filename, width = dims$width, height = dims$height, units = "px",
      bg = "transparent", type = "cairo", res = dpi)
  on.exit(dev.off())
  grid.newpage(recording = FALSE)
  grid.text(label, x = 0.5, y = 0.50,
            just = c("center", "center"),
            gp = gpar(fontsize = fontsize, col = "black", fontfamily = fontfamily))
  return(dims)
}

#render png labels and recover dimension info
label_dims <- cbind(w = rep(0, length(hypercube_names)),
                    h = rep(0, length(hypercube_names)))
for (i in seq_along(hypercube_names)) {
  name <- names(hypercube_names)[i]
  label <- hypercube_names[[name]]
  label_dims[i,] <- unlist(render_label_png(label, label_paths[i], fontsize = 48, dpi = 96, padding_factor = 0.25))
}

# Read the PNG, get dimension info
img_data <- readPNG(png_file)
img_height <- dim(img_data)[1]
img_width  <- dim(img_data)[2]
aspect     <- img_height / img_width  # or height/width

# Let's draw the quad so it’s about 2 units wide in X
quad_width  <- 1
quad_height <- quad_width * aspect

#### draw static pictures ####
use_rgl <- F
if(use_rgl){
  rgl_r_version <- function(){}
  # foo <- foreach(i=1:length(hcubes_proj), .packages=c("rgl")) %dopar% {
  # foo <- foreach(i=1:10, .packages=c("rgl")) %dopar% {
  for(i in 1:length(hcubes_proj)){
  # for(i in 1:200){
    
    #report on progress
    mcprint(paste0(i, " "))
    
    # Recover hcube vertices for this frame
    hcube_proj <- hcubes_proj[[i]]
    depth_proj <- depths_proj[[i]]
    
    #rgl setup
    if (i == 1) {
      open3d()
      bg3d("white")
      par3d(windowRect = c(0, 0, 800, 800))
      
      # Dummy points to fix bounding box
      dummy_pts <- expand.grid(
        x = c(-1, 1) * axis_bounds,
        y = c(-1, 1) * axis_bounds,
        z = c(-1, 1) * axis_bounds
      )
      points3d(dummy_pts$x, dummy_pts$y, dummy_pts$z, alpha = 0)
      
      # Remove default light and add overhead light
      clear3d("lights")  
      light3d(x = 2, y = 2, z = 2, 
              diffuse = "white",
              specular = "white",
              ambient = "white",
              viewpoint.rel = T
      )
      
      # set a view and lock
      view3d(theta = 0, phi = 0, fov = 0, zoom = 0.75)
      fixed_matrix <- par3d("userMatrix")
      
    }
    
    #set material properties
    material3d(
        color = face_col,
        alpha = opa[i],
        ambient = "black",
        specular = "white",
        emission = "black",
        shininess = 50,
        lit = TRUE,
        smooth = TRUE
      )
    points3d(dummy_pts$x, dummy_pts$y, dummy_pts$z, alpha = 0)
    
    #recover face coords in 3D
    verts <- hcube_proj[!duplicated(round(hcube_proj, 9)),]
    nverts <- nrow(verts)
    if(nverts >= 4){
      
      face_coords <- lapply(1:nrow(rot_hcube$faces), function(fi){
        f_inds <- rot_hcube$faces[fi, chull(hcube_proj[rot_hcube$faces[fi, ], ])]
        if(length(f_inds) == 4){
          return(cbind(hcube_proj[f_inds, ], depth_proj[f_inds]))  
        } else {
          return(NA)
        }
      })
      face_coords <- face_coords[!is.na(face_coords)]  
      
      #recover order of faces
      face_depths <- sapply(face_coords, function(x) mean(x[,3]))
      face_order <- order(face_depths, decreasing=T)  # farthest to nearest
      
    } else {
      #here we have a line, so we can just plot it once?
      face_coords <- list(verts)
      face_order <- 1
    }
    
    # Filter out any NAs or improperly formed faces
    valid_face_indices <- !sapply(face_coords, function(x) any(is.na(x)))
    face_coords <- face_coords[valid_face_indices]
    
    #recover only the unique face_coords
    face_coords <- unique(face_coords)
    
    #and also just the ones that are actual polygons and not eg points
    face_coords <- face_coords[sapply(face_coords, function(fcs){
      mean(apply(fcs, 2, var)) > 1E-4
    })]
    
    # plot faces
    
    #one face at a time
    # for (face_matrix in face_coords) {
    #   quads3d(face_matrix, col = current_face_col, alpha = 1)
    # }
    
    #all faces together
    if (length(face_coords) > 0) {
      # Combine all vertices into one matrix for quads3d
      all_face_vertices <- do.call(rbind, face_coords)
      
      # Determine face color for this frame
      current_face_col <- face_cols[min(i, length(face_cols))]
      current_seg_col <- seg_cols[min(i, length(seg_cols))]
      
      # Plot opaque quads. Color is recycled for all vertices.
      if(nrow(all_face_vertices) < 4){
        all_face_vertices <- cbind(all_face_vertices, 0)
        lines3d(all_face_vertices, col = current_seg_col, lwd = 4)  
      } else {
        quads3d(all_face_vertices)#, col = current_face_col, alpha = 1)
        #outline the faces
        num_quads <- nrow(all_face_vertices) / 4
        indices <- t(matrix(1:(num_quads * 4), nrow = 4)) # Rows: quad vertex indices (1,2,3,4), (5,6,7,8), ...
        seg_indices <- c(t(cbind(indices[,1], indices[,2], indices[,2], indices[,3],
                                 indices[,3], indices[,4], indices[,4], indices[,1])))
        segments3d(all_face_vertices[seg_indices, ], lwd = 4)#, col = current_seg_col) # Adjust color/lwd
      }
    }
    
    #write in text for shape
    
    #first recover text properties and location
    text_pos <- c(-axis_bounds, axis_bounds, axis_bounds)
    num_dim <- n_dim_per_tim[i]
    label_dim = label_dims[num_dim + 1,]
    label_ar <- label_dim["w"] / label_dim["h"]
    label_h_prop_bounds <- 0.075
    label_h <- label_h_prop_bounds * axis_bounds * 2
    label_bounds <- list(
      x = -axis_bounds + c(0, label_h, label_h, 0) * label_ar,
      y = axis_bounds - c(0, 0, label_h, label_h)
    )

    # Now add the textured quad
    quads3d(
      x = label_bounds$x,
      y = label_bounds$y,
      z = c(0,0,0,0),
      texcoords = list(
        x = c(0,1,1,0),
        y = c(1,1,0,0)
      ),
      texture = label_paths[num_dim+1],
      lit     = T,
      alpha   = 1,
      blend   = c("src_alpha","one_minus_src_alpha"),
      textype = "rgba",             # Use RGBA (alpha included)
      texmode = "replace",          # Don’t multiply texture by object color
      color   = "white"
    )
    
    #reset view angle
    view3d(theta = 0, phi = 0, fov = 0, zoom = 0.75, userMatrix = fixed_matrix)
    
    #plot to file
    output_filename <- paste0(frames_dir, paste0(rep(0, 5 - nchar(i)), collapse = ""), i, ".png")
    rgl.snapshot(filename = output_filename, fmt = "png")
    
    #now clear the plotting pane to prepare for the next one
    clear3d("shapes")
    
    if(i == length(hcubes_proj)){
     close3d() 
    }
    
  }
} else {
  base_r_version <- function(){}
  foo <- foreach(i=1:length(hcubes_proj), .packages=c("png")) %dopar% {
    
    #report on progress
    mcprint(paste0(i, " "))
    
    #recover hcube vertices
    hcube_proj <- hcubes_proj[[i]]
    depth_proj <- depths_proj[[i]]
    
    #write to file
    png(filename = paste0(frames_dir, paste0(rep(0, 5 - nchar(i)), collapse = ""), i, ".png"), 
        width = 800, height = 800, type = "cairo", pointsize = 35)
    
    plot.new()
    par(mar = c(0,0,0,0))
    plot.window(xlim = c(-axis_bounds, axis_bounds),
                ylim = c(-axis_bounds, axis_bounds), asp = 1)
    # plot.window(xlim = c(-max_val[i], max_val[i]),
    #             ylim = c(-max_val[i], max_val[i]), asp = 1)
    
    
    #recover face coords in 3D
    verts <- hcube_proj[!duplicated(round(hcube_proj, 9)),]
    nverts <- nrow(verts)
    if(nverts >= 4){
      
      face_coords <- lapply(1:nrow(rot_hcube$faces), function(fi){
        f_inds <- rot_hcube$faces[fi, chull(hcube_proj[rot_hcube$faces[fi, ], ])]
        if(length(f_inds) == 4){
          return(cbind(hcube_proj[f_inds, ], depth_proj[f_inds]))  
        } else {
          return(NA)
        }
      })
      face_coords <- face_coords[!is.na(face_coords)]  
      
      #recover order of faces
      face_depths <- sapply(face_coords, function(x) mean(x[,3]))
      face_order <- order(face_depths, decreasing=T)  # farthest to nearest
      
    } else {
      #here we have a line, so we can just plot it once?
      face_coords <- list(verts)
      face_order <- 1
    }
    
    
    #plot faces of the polygon
    for(f in face_order) {
      polygon(face_coords[[f]][,1:2], col = face_cols[min(i, length(face_cols))], 
              border = seg_cols[min(i, length(seg_cols))], 
              lwd = 2, xpd = NA)
    }
    
    #plot edges of the polygon
    # segments(x0 = hcube_proj[rot_hcube$edges[,1],1], 
    #          y0 = hcube_proj[rot_hcube$edges[,1],2], 
    #          x1 = hcube_proj[rot_hcube$edges[,2],1], 
    #          y1 = hcube_proj[rot_hcube$edges[,2],2], col = seg_cols[i], lwd = 2)
    
    dev.off()
  }
  
}

#check a single frame gen
# for(fii in 1:length(face_coords)){
#   
#   plot(face_coords[[1]][,1:2], asp = 1, xlim = range(do.call(rbind, face_coords)[,1]),
#        ylim = range(do.call(rbind, face_coords)[,2]))
#   for(fi in 1:length(face_coords)){
#     lines(rbind(face_coords[[fi]][,1:2], face_coords[[fi]][1,1:2]), asp = 1)
#   }
#   
#   polygon(face_coords[[fii]][,1:2], col = face_cols[min(i, length(face_cols))], 
#           border = seg_cols[min(i, length(seg_cols))], 
#           lwd = 2, xpd = NA)
#   
# }
# 


#### animate ####
# stitch together animation
fps <- 120
thin <- 1

output_temp <- "rotation_tmp.mp4"
output_final <- "rotation.mp4"

# Generate animation video (no reverse, no concat)
system(paste0(
  "cd ", vid_dir, "; ",
  "ffmpeg -y -r ", fps / thin, " -f image2 -s 800x800 ",
  "-i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ", output_temp
))

# Atomic rename into place
file.rename(
  file.path(vid_dir, output_temp),
  file.path(vid_dir, output_final)
)

# Optional: play the video
cat("mpv ", file.path(vid_dir, output_final), "\n")

