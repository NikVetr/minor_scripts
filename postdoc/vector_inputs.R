my_sum <- function(x, y, z) {
    arg_vals <- as.list(environment())
  nargs <- length(arg_vals)
  arg_lens <- sapply(arg_vals, length)
  if(any(arg_lens > 1)){
    out <- numeric(max(arg_lens))
    for(i in 1:max(arg_lens)){
      arg_inds <- (i - 1) %% arg_lens + 1
      arg_vals_i <- lapply(setNames(1:nargs, names(arg_vals)), function(j) arg_vals[[j]][arg_inds[j]])
      out[i] <- do.call(my_sum, arg_vals_i)
    }
    return(out)
  } else {
    return(x + y + z)  
  }
}

my_sum(1,2,3:5)

my_sum2 <- function(x, y, z, ...) {
  arg_vals <- as.list(environment())
  arg_vals <- c(arg_vals, list(...))
  nargs <- length(arg_vals)
  arg_lens <- sapply(arg_vals, length)
  if(any(arg_lens > 1)){
    out <- numeric(max(arg_lens))
    for(i in 1:max(arg_lens)){
      arg_inds <- (i - 1) %% arg_lens + 1
      arg_vals_i <- lapply(setNames(1:nargs, names(arg_vals)), function(j) arg_vals[[j]][arg_inds[j]])
      out[i] <- do.call(my_sum2, arg_vals_i)
    }
    return(out)
  } else {
    return(sum(x, y, z, ...))  
  }
}

my_sum2(1, 2, 3, 10:100)
