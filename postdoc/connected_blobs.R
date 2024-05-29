pixel_coords <- matrix(c(2, 4,
                         3, 4,
                         1, 3,
                         2, 3,
                         3, 3,
                         4, 3,
                         1, 2,
                         2, 2,
                         3, 2,
                         4, 2,
                         2, 1,
                         3, 1,
                         7, 3,
                         8, 3,
                         9, 3,
                         6, 2,
                         7, 2,
                         8, 2), 
                       ncol = 2, byrow = TRUE)
colnames(pixel_coords) <- c("x","y")

#retrieve adjacency matrix
adj <- (as.matrix(dist(pixel_coords, method = "manhattan")) == 1) * 1
#or for diagonal connectedness, can do eg
#adj <- (as.matrix(dist(pixel_coords, method = "euclidean")) <= (sqrt(2) + 1E-6)) * 1
adj[adj==0] <- NA

#find shortest paths w/ floyd-warshall
shortest_path <- Rfast::floyd(adj)
shortest_path[shortest_path > 0] <- 1
diag(shortest_path) <- 1

#find connected blobs and isolate largest one(s)
n_connected <- rowSums(shortest_path, na.rm = T)
largest_blob_size <- max(n_connected)
largest_blob_indices <- which(n_connected == largest_blob_size)
n_largest_blobs <- length(largest_blob_indices) / largest_blob_size
split_largest_blob_indices <- lapply(1:n_largest_blobs, function(i) integer(0))
for(i in 1:n_largest_blobs){
  split_largest_blob_indices[[i]] <- which(shortest_path[largest_blob_indices[1],] == 1)
  largest_blob_indices <- setdiff(largest_blob_indices, split_largest_blob_indices[[i]])
}

#plot the result
plot(pixel_coords)
points(pixel_coords[split_largest_blob_indices[[1]],], col = 2, pch= 19) #plot first largest blob


#alternatively, with igraph
library(igraph)

# Assuming 'adj' is your adjacency matrix
adj <- (as.matrix(dist(pixel_coords, method = "manhattan")) == 1)
#or for diagonal connectedness, can do eg
#adj <- (as.matrix(dist(pixel_coords, method = "euclidean")) <= (sqrt(2) + 1E-6)) * 1
g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
blobs <- components(g)
blob_sizes <- components$csize
largest_blob_size <- max(blob_sizes)
largest_blob_indices <- which(component_sizes == largest_blob_size)
split_largest_blob_indices <- lapply(largest_blob_indices, function(x) {
  which(blobs$membership == x)
})


#let's have a larger size?
n_blobs <- 1E2
n_pix_per_blob <- sample(x = 10:1E2, n_blobs, replace = T)
blob_seeds <- cbind(x = sample(1:1E7, n_blobs, replace = F),
                   y = sample(1:1E7, n_blobs, replace = F))
n_expansions <- ceiling(sqrt(max(n_pix_per_blob)))
displacement_mats <- lapply(1:n_expansions, function(i){
  expand.grid(-ceiling(i/2):ceiling(i/2), -ceiling(i/2):ceiling(i/2))
})

blob_expansion_n <- choose(1:10, 2)
pixel_coords <- do.call(rbind, lapply(1:n_blobs, function(i){
  t(blob_seeds[i,] + t(displacement_mats[[ceiling(sqrt(n_pix_per_blob[i]))]]))
}))
pixel_coords <- pixel_coords[sample(1:nrow(pixel_coords)),]
adj <- (as.matrix(dist(pixel_coords, method = "manhattan")) == 1)

floyd_version <- function(adj){
  adj[adj==0] <- NA
  
  #find shortest paths w/ floyd-warshall
  shortest_path <- Rfast::floyd(adj)
  shortest_path[shortest_path > 0] <- 1
  diag(shortest_path) <- 1
  
  #find connected blobs and isolate largest one(s)
  n_connected <- rowSums(shortest_path, na.rm = T)
  largest_blob_size <- max(n_connected)
  largest_blob_indices <- which(n_connected == largest_blob_size)
  n_largest_blobs <- length(largest_blob_indices) / largest_blob_size
  split_largest_blob_indices <- lapply(1:n_largest_blobs, function(i) integer(0))
  for(i in 1:n_largest_blobs){
    split_largest_blob_indices[[i]] <- which(shortest_path[largest_blob_indices[1],] == 1)
    largest_blob_indices <- setdiff(largest_blob_indices, split_largest_blob_indices[[i]])
  }
  return(split_largest_blob_indices)
}

igraph_version <- function(adj){
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  blobs <- components(g)
  largest_blob_size <- max(blobs$csize)
  largest_blobs <- which(blobs$csize == largest_blob_size)
  largest_blob_indices <- which(blobs$membership %in% largest_blobs)
  split_largest_blob_indices <- split(largest_blob_indices, blobs$membership[largest_blob_indices])
  return(split_largest_blob_indices)
}

microbenchmark::microbenchmark(igraph_version(adj),
                               floyd_version(adj), 
                               times = 1)
