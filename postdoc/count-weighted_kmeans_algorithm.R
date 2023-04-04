wkmeans <- function(x, w, k, kmeans.iter.max = NA, n = NA, eps = 1E-6){
  
  #hyperparameters
  if(is.na(n)){
    n <- k * 10
  }
  
  if(is.na(kmeans.iter.max)){
    kmeans.iter.max <- k * 5  
  }
  
  kmeans_out <- lapply(1:n, function(kmeans_run){
    
    #initialize centers w/ kmeans++
    #TODO need to weigh probs by count!
    kmeans_iter <- 1
    converged <- F
    cents <- matrix(NA, k, 4)
    cents[1,] <- x[sample(1:nrow(x), 1),]
    for(ki in 2:k){
      sqdists <- sqdists_to_centers(x, cents[1:(ki-1),])
      cls <- apply(sqdists, 1, which.min)
      smallest_sqdists <- sqdists[cbind(1:length(cls), cls)]
      cents[ki,] <- x[sample(1:length(cls), 1, prob = smallest_sqdists),]
    }
    
    #evaluate initial cluster assignments
    sqdists <- sqdists_to_centers(x, cents)
    cls <- apply(sqdists, 1, which.min)
    
    #run kmeans iteratively
    while(!converged & kmeans_iter < kmeans.iter.max){
      
      #increment iter
      kmeans_iter <- kmeans_iter + 1
      
      #check that all clusters present and reassign cents if not
      prev_cents <- cents
      if(!all(1:k %in% cls)){
        empty_cls <- (1:k)[!(1:k %in% cls)]
        for(ki in empty_cls){
          sqdists <- sqdists_to_centers(x, cents)
          cls <- apply(sqdists, 1, which.min)
          smallest_sqdists <- sqdists[cbind(1:length(cls), cls)]
          cents[ki,] <- x[sample(1:length(cls), 1, prob = smallest_sqdists),]
        }
        sqdists <- sqdists_to_centers(x, cents)
        cls <- apply(sqdists, 1, which.min)
      }
      
      #find new centers
      cents <- do.call(rbind, lapply(1:k, function(ki){
        clusmems <- which(cls == ki) #members of the cluster
        clusvals <- x[clusmems,] 
        cluscounts <- w[clusmems]
        apply(t(clusvals) * rep(cluscounts, each = 4), 1, sum) / sum(cluscounts)
      }))
      
      #check for convergence
      converged <- all(abs(cents - prev_cents) < eps)
      
      #evaluate new cluster assignments
      if(!converged){
        sqdists <- sqdists_to_centers(x, cents)
        cls <- apply(sqdists, 1, which.min)
      }
      
    }
    
    #calculate total variance / sum of squares
    total_ss <- sum(sqdists[cbind(1:length(cls), cls)] * w)
    
    #return value
    list(converged = converged, k = k, kmeans_iter = kmeans_iter, cents = cents, cls = cls, total_ss = total_ss, obs_counts = w)
    
  })
  
  #pick best run
  total_sss <- sapply(kmeans_out, function(kmeans_indiv) kmeans_indiv$total_ss)
  optimal_cluster_index <- which.min(total_sss)
  optimal_cluster <- kmeans_out[[optimal_cluster_index]]
  return(optimal_cluster)
  
}

wkmeans(x, w, k)
