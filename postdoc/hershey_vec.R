source("~/scripts/minor_scripts/postdoc/true_text_dim_functions.R")
source("~/repos/polylines/R/functions.R")
source("~/scripts/minor_scripts/postdoc/interpolate_signature_frequency_series.R")

# Read the input text file
font_name <- "HersheyScript"
# string <- paste0(letters[5:15], collapse = "")
# string <- paste0(letters, collapse = "")
# string <- "a pacharinsak"
# string <- "a gates"
# string <- "a morrill" 
string_chars <- strsplit(string, "")[[1]]

hershey_vec <- get_hershey_coords(string = string, font_name = font_name)
hershey_vecs <- lapply(string_chars, function(string_char){
  get_hershey_coords(string = string_char, font_name = font_name)
})
hershey_ncoords <- sapply(hershey_vecs, nrow)
hershey_inds <- rep(1:length(hershey_ncoords), hershey_ncoords)
split_vec <- split(hershey_vec, hershey_inds)

#find individual strokes
tol <- 1E-3
split_vec <- lapply(split_vec, function(sv){
  diffs <- apply(sv, 2, diff)
  sv$diff <- 0
  sv$diff[-1] <- sqrt(rowSums(diffs^2))
  sv$ss <- sv$zd <- abs(sv$diff) < tol #ss is for "same segment", zd for "zero diff"
  sv$zd[1] <- F
  sv$ss[which(sv$ss) + 1] <- T
  
  #when we have a run of F, we expect the next F after the first to be the continuation
  rless <- rle(sv$ss)
  rless_starts <- c(0, cumsum(rless$lengths))[1:length(rless$values)] + 1
  first_F <- rless_starts[!rless$values]
  sv$ss[first_F + 1] <- T
  
  #identify segments and split
  sv$seg <- cumsum(!sv$ss) + 1
  sv <- sv[!sv$zd,]
  return(split(sv, sv$seg))
})
split_vec <- unlist(split_vec, recursive = F)
nc <- length(split_vec)

#connect ends
first_last <- lapply(split_vec, function(sv){
  sv[c(1,nrow(sv)), c("x", "y")]
})
fl_cat <- do.call(rbind, first_last)
fl_dists <- fields::rdist(fl_cat[1:nc*2,],
                          fl_cat[1:nc*2-1,])
continues <- which(fl_dists < 1E-6, arr.ind = T)
continues <- continues[continues[,1] != continues[,2],]
comp <- igraph::components(igraph::graph_from_data_frame(continues, directed = FALSE))$membership
cont_sets <- split(as.numeric(names(comp)), comp)
indep_components <- setdiff(1:nc, unlist(cont_sets))
all_sets <- c(indep_components, cont_sets)
all_sets <- all_sets[order(sapply(all_sets, head, 1))]
split_vec <- lapply(all_sets, function(si){
  vecs <- split_vec[si]
  num_comps <- length(si)
  if(num_comps > 1){
    vecs[2:num_comps] <- lapply(2:num_comps, function(ci){
      vecs[[ci]][-1,]
    })
    return(do.call(rbind, vecs))
  } else {
    return(vecs[[1]])
  }
})

#smooth output
desired_num_points <- 200
smoothing_variance_threshold <- 0.9999
smooth_split_vec <- lapply(seq_along(split_vec), function(i){
  
  #interpolate coords
  coord_matrix <- split_vec[[i]]
  
  #jitter points first
  jitter_sd <- sqrt((var(coord_matrix$x) + var(coord_matrix$y)) / 1E9)
  coord_matrix[,c("x", "y")] <- coord_matrix[,c("x", "y")] + 
    rnorm(n = nrow(coord_matrix) * 2, sd = jitter_sd)
  
  cmeans <- colMeans(coord_matrix[,1:2])
  coord_matrix$x <- coord_matrix$x - cmeans["x"]
  coord_matrix$y <- coord_matrix$y - cmeans["y"]
  interpolated_matrix_with_i <- interpolate_to_length_kd(
    x = coord_matrix,
    n = desired_num_points,
    equal_increments = TRUE
  )
  
  wf <- c("haar", "la8", "d4", "d6", "mb4")[5]
  #smooth coords
  interpolated_xy <- interpolated_matrix_with_i[, 1:2, drop = FALSE]
  x_interpolated <- interpolated_xy[, 1]
  x_smooth_reflected <- wavelet_smooth(c(x_interpolated, rev(x_interpolated)),
                                       var_thresh = smoothing_variance_threshold, 
                                       wavelet_filter = wf)
  # x_smooth_reflected <- fourier_smooth(c(x_interpolated, rev(x_interpolated)))
  x_final <- unlist(x_smooth_reflected)[1:desired_num_points]
  
  y_interpolated <- interpolated_xy[, 2]
  y_smooth_reflected <- wavelet_smooth(c(y_interpolated, rev(y_interpolated)),
                                       var_thresh = smoothing_variance_threshold, 
                                       wavelet_filter = wf)
  # y_smooth_reflected <- fourier_smooth(c(y_interpolated, rev(y_interpolated)))
  y_final <- unlist(y_smooth_reflected)[1:desired_num_points]
  
  #return result
  out <- data.frame(x = x_final + cmeans["x"], y = y_final + cmeans["y"])
  return(out)
  
})

#plot output
plot_test <- F
if(plot_test){
  plot.new()
  plot.window(xlim = range(hershey_vec$x), ylim = range(hershey_vec$y), asp = 1)
  
  for(i in 1:length(split_vec)){
    # lines(x = split_vec[[i]]$x, y = split_vec[[i]]$y, col = i, lwd = 4)
    lines(x = smooth_split_vec[[i]]$x, y = smooth_split_vec[[i]]$y, col = i, lwd = 2)
    # polylines(x = smooth_split_vec[[i]]$x, y = smooth_split_vec[[i]]$y, lwd = 3)
  }
}
#record for later
orig_smooth_vec <- smooth_split_vec
