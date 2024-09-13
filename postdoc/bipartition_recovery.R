library(ape)
library(ComplexHeatmap)
library(phangorn)

# Function to compute the bipartition of a tree
bipart_exists <- function(r1, r2, p1, p2 = 2, use_tree = T){
  R1 <- diag(p1) + r1 - diag(p1) * r1
  R2 <- diag(p1+p2)
  R2[1:p1, 1:p1] <- R1
  R2[1,p1+1] <- R2[p1+1,1] <- r2
  
  d <- 1-abs(R2)
  colnames(d) <- rownames(d) <- paste0("v", 1:(p1+p2))
  tree <- nj(d)
  net <- neighborNet(d)
  plot(net)
  
  if(use_tree){
    biparts <- prop.part(tree)  
  } else {
    biparts <- prop.part(net)
  }
  
  bipart.exists <- any(vapply(biparts, function(x) length(x) == p && all(x == 1:p), 
                              FUN.VALUE = T))
  return(bipart.exists)
}

rvals <- 0:100/100
p <- 5
bp.exists <- do.call(rbind, parallel::mclapply(rvals, function(r1) 
  vapply(rvals, function(r2)
    bipart_exists(r1, r2, p), FUN.VALUE = T), 
  mc.cores = 12))
line_of_transition <- apply(bp.exists, 1, function(x) min(which(!x)))
plot(rvals, rvals[line_of_transition], xlim = c(0,1), ylim = c(0,1), type = "l")
plot(which(bp.exists, arr.ind = T), xlim = c(0,101), ylim = c(0,101))
# Heatmap(bp.exists * 1)
