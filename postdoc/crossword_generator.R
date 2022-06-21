n_words <- 20
n_letters <- table(sample(5:10, n_words, replace = T))
possible_words <- words::words
words <- unlist(sapply(as.integer(names(n_letters)), function(nl) sample(possible_words$word[possible_words$word_length == nl], 
                                               n_letters[as.character(nl)], replace = F)))
total_nchar <- nchar(paste0(words, collapse = ""))

xword <- matrix("", ncol = max(sqrt(total_nchar * 5), nchar(words)), nrow = max(sqrt(total_nchar * 5), nchar(words)))
used <- rep(F, n_words)

init_xword <- function(xword, words, used){
  tp <- sample(which(!used), size = 1)
  nr <- nrow(xword)
  nc <- ncol(xword)
  chars <- strsplit(words[tp], "")[[1]]
  start <- (floor(nc/2) - floor(length(chars) / 2))
  end <- start + length(chars) - 1
  start <- start - max(0, end - nc)
  end <- end - max(0, end - nc)
  xword[floor(nr/2), start:end] <- chars
  used[tp] <- T
  return(list(xword, used))
}

initialized_xword <- init_xword(xword, words, used)
xword <- initialized_xword[[1]]
used <- initialized_xword[[2]]
  
add_xword <- function(xword, words, used){
  tp <- sample(which(!used), size = 1)
  nr <- nrow(xword)
  nc <- ncol(xword)
  chars <- strsplit(words[tp], "")[[1]]
  matches <- sapply(setNames(chars, chars), function(let) which(let == xword, arr.ind = T))
  matches <- cbind(do.call(rbind, matches), char_i = rep(1:length(chars), sapply(matches, function(x) nrow(x))))
  
  if(nrow(matches) == 0){
    return(list(xword, used))
  }
  
  horiz_scores <- sapply(1:nrow(matches), function(ri) check_compatibility_horiz(chars, xword, matches[ri,]))
  vert_scores <- sapply(1:nrow(matches), function(ri) check_compatibility_vert(chars, xword, matches[ri,]))
  
  if(all(c(horiz_scores, vert_scores) < 1)){
    return(list(xword, used))
  }
  
  add_dir <- c("horiz", "vert")[which.max(c(max(horiz_scores), max(vert_scores)))]
  
  if(add_dir == "vert"){
    loc <- matches[which.max(vert_scores),]
    xword[(loc[1] - loc[3] + 1) : (loc[1] - loc[3] + length(chars)), loc[2]] <- chars
  } else {
    loc <- matches[which.max(horiz_scores),]
    xword[loc[1], (loc[2] - loc[3] + 1) : (loc[2] - loc[3] + length(chars))] <- chars
  }
  
  used[tp] <- T
  
  return(list(xword, used))
  
}


check_compatibility_vert <- function(chars, xword, loc){
  xlocs <- (loc[1] - loc[3] + 1) : (loc[1] - loc[3] + length(chars))
  if(any(xlocs > ncol(xword)) | any(xlocs < 1)){
    return(-1)
  }  
  xword_word = xword[xlocs, loc[2]]
  if(any(!(xword_word == chars | xword_word == ""))){
    return(-1)
  } else {
    return(sum(xword_word == chars))
  }
}

check_compatibility_horiz <- function(chars, xword, loc){
  ylocs <- (loc[2] - loc[3] + 1) : (loc[2] - loc[3] + length(chars))
  if(any(ylocs > nrow(xword)) | any(ylocs < 1)){
    return(-1)
  }
  xword_word = xword[loc[1], ylocs]
  if(any(!(xword_word == chars | xword_word == ""))){
    return(-1)
  } else {
    return(sum(xword_word == chars))
  }
}


add_xword(xword, words, used)

generate_xword <- function(words, tolerance = 5){
  xword <- matrix("", ncol = max(sqrt(total_nchar * 10), nchar(words)), nrow = max(sqrt(total_nchar * 10), nchar(words)))
  used <- rep(F, n_words)  
  
  initialized_xword <- init_xword(xword, words, used)
  xword <- initialized_xword[[1]]
  used <- initialized_xword[[2]]
  
  dead_ends <- 0
  while(dead_ends < tolerance & !all(used)){
    curr_used <- used
    new_xword <- add_xword(xword, words, used)
    xword <- new_xword[[1]]
    used <- new_xword[[2]]
    if(all(curr_used == used)){dead_ends <- dead_ends + 1}
  }
  
  return(xword)
}

trim_matrix <- function(xword){
  empty_rows <- sapply(1:nrow(xword), function(ri) all(xword[ri,] == ""))
  empty_cols <- sapply(1:ncol(xword), function(ci) all(xword[,ci] == ""))
  xword <- xword[!empty_rows, !empty_cols]
  return(xword)
}

lapply(1:20, function(x) trim_matrix(generate_xword(words, tolerance = 50)))
