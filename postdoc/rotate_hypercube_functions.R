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