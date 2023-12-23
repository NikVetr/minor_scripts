library(hash)

#randomly generate xword entries
n_words <- 50
n_letters <- table(sample(5:10, n_words, replace = T))
possible_words <- words::words
words <- unlist(sapply(as.integer(names(n_letters)), function(nl) sample(possible_words$word[possible_words$word_length == nl], 
                                               n_letters[as.character(nl)], replace = F)))

#or read in your own
words_and_clues <- readLines("~/2023-Dec_Xmas-Card-Xword.txt", warn = F)
words_and_clues <- readLines("~/2023-Dec_Xmas-Card-Xword_Russian.txt", warn = F, encoding = "UTF-8")

words_and_clues <- words_and_clues[grepl(":", words_and_clues)]
words_and_clues <- data.frame(trimws(do.call(rbind, strsplit(words_and_clues, split = ":"))))
colnames(words_and_clues) <- c("clues", "words")
words <- tolower(gsub(" ", "", words_and_clues$words))
n_words <- length(words)
clues <- setNames(words_and_clues$clues, words)
full_words <- setNames(words_and_clues$words, words)

total_nchar <- nchar(paste0(words, collapse = ""))
n_words <- length(words)
xword <- matrix("", ncol = ceiling(max(sqrt(total_nchar * 5), nchar(words))) * 2, 
                nrow = ceiling(max(sqrt(total_nchar * 5), nchar(words))) * 2)
used <- setNames(rep(F, n_words), words)
word_locs <- lapply(setNames(words, words), function(x) return(NA))

init_xword <- function(xword, words, used, word_locs){
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
  word_locs[[tp]] <- cbind(floor(nr/2), start:end)
  return(list(xword = xword, used = used, word_locs = word_locs))
}

initialized_xword <- init_xword(xword, words, used, word_locs)
xword <- initialized_xword[["xword"]]
used <- initialized_xword[["used"]]
word_locs <- initialized_xword[["word_locs"]]


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

find_words <- function(xword){
  
  horiz_words <- unlist(apply(xword, 1, function(x){
    rles <- rle(x != "")
    rlebs <- c(0, cumsum(rles$lengths))
    rlebs <- cbind(rlebs[-length(rlebs)]+1, rlebs[-1])
    word_inds <- matrix(rlebs[rles$values & rles$lengths != 1,], ncol = 2)
    if(nrow(word_inds) != 0){
      out <- sapply(1:nrow(word_inds), function(wi) paste0(x[word_inds[wi,1]:word_inds[wi,2]], collapse = ""))
    } else {
      out <- character(0)
    }
    return(out)
  }))
  
  vert_words <- unlist(apply(xword, 2, function(x){
    rles <- rle(x != "")
    rlebs <- c(0, cumsum(rles$lengths))
    rlebs <- cbind(rlebs[-length(rlebs)]+1, rlebs[-1])
    word_inds <- matrix(rlebs[rles$values & rles$lengths != 1,], ncol = 2)
    if(nrow(word_inds) != 0){
      out <- sapply(1:nrow(word_inds), function(wi) paste0(x[word_inds[wi,1]:word_inds[wi,2]], collapse = ""))
    } else {
      out <- character(0)
    }
    return(out)
  }))
  
  return(c(horiz_words, vert_words))
}

find_n_extra <- function(xword, words){
  found_words <- find_words(xword)
  bad_words <- setdiff(found_words, words)
  missing_words <- setdiff(words, found_words)
  
  n_unmatched_bad_words <- 0
  if(length(bad_words) == 0 & length(missing_words) == 0){
    n_extra <- 0
  } else {
    
    if(length(missing_words) > 0){
      n_extra <- unlist(sapply(missing_words, function(word){
        subword_found <- grepl(word, bad_words)
        extra_letters <- nchar(bad_words[subword_found]) - nchar(word)
        if(length(extra_letters) == 0) return(0) else return(extra_letters)
      }))
    } else {
      n_extra <- rep(0, max(c(length(missing_words)), 1))
    }
    
    if(length(missing_words) > 0){
      n_unmatched_bad_words <- sapply(bad_words, function(bad_word){
          if(any(sapply(missing_words, function(word) grepl(word, bad_word)))){
            return(0)
          } else {
            return(nchar(bad_word))
          }
        })
    } else {
      n_unmatched_bad_words <- nchar(bad_words)
    }
  }
  
  if(length(n_unmatched_bad_words) == 0){
    n_unmatched_bad_words <- 0
  }
  
  return(c(n_extra, unlist(n_unmatched_bad_words)))
}

quick_find_n_extra <- function(xword, loc, chars, add_dir){
  xword[loc[1], loc[2]]
  if(add_dir == "vert"){
    charlocs <- ((loc[1] - loc[3] + 1) : (loc[1] - loc[3] + length(chars)))[-loc[3]]
    adj_letters <- c(xword[charlocs, loc[2]-1], 
                     xword[charlocs, loc[2]+1], 
                     xword[c(charlocs[1]-1, tail(charlocs, 1) + 1), loc[2]]) != ""
  } else {
    charlocs <- (loc[2] - loc[3] + 1) : (loc[2] - loc[3] + length(chars))[-loc[3]]
    adj_letters <- c(xword[loc[1]-1, charlocs],
                     xword[loc[1]+1, charlocs],
                     xword[loc[1], c(charlocs[1]-1, tail(charlocs, 1) + 1)]) != ""
  }
  return(sum(adj_letters))
}

horiz_vert <- function(word_locs){
  sapply(word_locs[!sapply(word_locs, function(wlx) all(is.na(wlx)))], function(wli) 
    c("h","v")[which(apply(wli, 2, function(wlinds) wlinds[1] == wlinds[2]))])
}

score_xword <- function(xword, words, word_locs, pic_locs = NA){
  
  n_words_filled <- sum(!sapply(word_locs, function(x) all(is.na(x))))
  if(n_words_filled == 1){
    return(0)
  }
  
  trimmed_xword_info <- trim_matrix(xword, word_locs)
  trimmed_xword <- trimmed_xword_info$xword
  trimmed_word_locs <- trimmed_xword_info$word_locs
  
  #score xword
  penalties <- list(
    
    #penalty #1 is the total area of the xword
    area = prod(dim(trimmed_xword)),
    
    #penalty #2 is the deviation from the desired dimensional ratio
    area_ratio = {
      desired_area_ratio <- 1
      ratio_ratio <- desired_area_ratio / (ncol(trimmed_xword) / nrow(trimmed_xword))
      max(c(ratio_ratio, 1/ratio_ratio))},
    
    #3 is the deviation from the desired vertical vs horiz word ratio
    hv_ratio = {
      desired_hv_ratio <- 1
      obs_hv_count <- table(horiz_vert(word_locs))
      hv_ratio_ratio <- desired_hv_ratio / (obs_hv_count["h"] / obs_hv_count["v"])
      if(all(is.na(hv_ratio_ratio))){
        0
      } else {
        max(c(hv_ratio_ratio, 1/hv_ratio_ratio))  
      }
    },
    
    #4 is the connectivity of the matrix
    n_connections = {
      filled_words <- !sapply(word_locs, function(wlx) all(is.na(wlx)))
      ifelse(sum(filled_words > 4), 
        apply(xword_adjacency_matrix(word_locs[filled_words]), 2, sum),
        0)
      },
    
    #5 is the number of 4 & 6 length cycles in the xword
    n_cycles = {
      filled_words <- !sapply(word_locs, function(wlx) all(is.na(wlx)))
      if(sum(filled_words > 4)){
        adjmat <- xword_adjacency_matrix(word_locs[filled_words])
        g <- igraph::graph_from_adjacency_matrix(adjmat, mode = "undirected")
        sapply(c(4,6), function(cycle_size){
          n_cycles <- igraph::graph.motifs(g, size = cycle_size)
          sum(n_cycles, na.rm = T)
        })
      } else {
        0
      }
    },
    
    #6 is how circular the overall shape is
    #find convex hull, then find variance of distances from center
    #the lower the variance, the more circular the shape
    circularity = {
      letter_locs <- do.call(rbind, trimmed_word_locs)
      letter_locs <- letter_locs[!is.na(letter_locs[,1]),]
      letter_locs <- letter_locs[!duplicated(letter_locs),]
      # convex_hull <- letter_locs[chull(letter_locs),]
      concave_hull <- concaveman::concaveman(letter_locs, concavity = 1)
      center <- apply(concave_hull, 2, mean)
      letter_dists <- sqrt(apply((t(t(concave_hull) - center))^2, 1, sum))
      var(letter_dists)
    },
    
    #7 compares the xword to a user-supplied picture
    pic_similarity = {
      if(all(is.na(pic_locs))){
        0
      } else {
        char_locs <- fill_adj_elements(xword = trimmed_xword, word_locs = trimmed_word_locs)$char_locs
        pic_comparison_dists <- evaluate_picture_concordance(char_locs, pic_locs)
        unlist(pic_comparison_dists)
      }
    }
      
  )
  
  #positive penalties are bad
  #negative are good
  penalty <- sqrt(penalties[["area"]]) + # + makes xword smaller / more compact
    (penalties[["area_ratio"]])^2 + # + makes xword minimum bounds more square-like
    (penalties[["hv_ratio"]])^3 + # + makes xword have more similar #s of horiz vs vert words
    -sum(penalties[["n_connections"]]^2) + # + makes xword have fewer word intersections
    -sum(penalties[["n_cycles"]]^1) + # + makes xword have fewer cycles of degree 4 and 6
    sum(penalties[["pic_similarity"]]^4) + # + makes xword more like the picture supplied
    penalties[["circularity"]]^2 # + makes xword more circular
  return(penalty)
  
}

softmax <- function(x) exp(x) / sum(exp(x))

add_xword <- function(xword, words, used, word_locs){
  
  unused_inds <- which(!used)
  if(length(unused_inds) == 1){
    tp <- unused_inds
  } else {
    tp <- sample(unused_inds, size = 1)  
  }
  
  nr <- nrow(xword)
  nc <- ncol(xword)
  chars <- strsplit(words[tp], "")[[1]]
  matches <- lapply(setNames(chars, chars), function(let) which(let == xword, arr.ind = T))
  matches <- data.frame(cbind(do.call(rbind, matches), 
                              char_i = rep(1:length(chars), sapply(matches, function(x) nrow(x)))))
  
  if(nrow(matches) == 0){
    return(list(xword = xword, used = used, word_locs = word_locs, success = F))
  }
  
  #evaluate compatibility with existing xword
  horiz_scores <- sapply(1:nrow(matches), function(ri){
    check_compatibility_horiz(chars, xword, unlist(matches[ri,]))
  })

  vert_scores <- sapply(1:nrow(matches), function(ci){
    check_compatibility_vert(chars, xword, unlist(matches[ci,]))
  })

  matches$easy_check <- !(horiz_scores == -1 & vert_scores == -1)
  
  if(all(!matches$easy_check)){
    return(list(xword = xword, used = used, word_locs = word_locs, success = F))
  }
  
  matches$direction <- NA
  matches$direction[matches$easy_check] <- c("horiz", "vert")[apply(cbind(horiz_scores[matches$easy_check], vert_scores[matches$easy_check]), 
                                                                    1, 
                                                                    function(x) which(x != -1)[1])]
  
  #now do harder check for passing words
  matches$hard_check <- matches$easy_check
  matches$score <- Inf
  poss_scores <- data.frame(do.call(rbind, lapply(which(matches$easy_check), function(mi){
    add_dir <- matches$direction[mi]
    loc <- unlist(matches[mi,1:3])
    temp_xword <- xword
    temp_word_locs <- word_locs
    if(add_dir == "vert"){
      temp_xword[(loc[1] - loc[3] + 1) : (loc[1] - loc[3] + length(chars)), loc[2]] <- chars
      temp_word_locs[[tp]] <- cbind((loc[1] - loc[3] + 1) : (loc[1] - loc[3] + length(chars)), loc[2])
    } else {
      temp_xword[loc[1], (loc[2] - loc[3] + 1) : (loc[2] - loc[3] + length(chars))] <- chars
      temp_word_locs[[tp]] <- cbind(loc[1], (loc[2] - loc[3] + 1) : (loc[2] - loc[3] + length(chars)))
    }
    c(n_extra = sum(find_n_extra(temp_xword, words)), score = score_xword(xword = temp_xword, words, word_locs = temp_word_locs))
    # c(n_extra = quick_find_n_extra(xword, loc, chars, add_dir), score = score_xword(temp_xword, words, temp_word_locs))
  })))
  matches$hard_check[matches$easy_check] <- poss_scores$n_extra == 0
  matches$score[matches$easy_check] <- poss_scores$score
  
  if(all(!matches$hard_check)){
    return(list(xword = xword, used = used, word_locs = word_locs, success = F))
  }
  
  #now add word in that passes both checks
  good_matches <- matches[matches$hard_check,] 
  # winning_match <- which.min(good_matches$score)
  winning_match <- sample(nrow(good_matches), 1, prob = softmax(-(good_matches$score - min(good_matches$score))))
  
  add_dir <- good_matches$direction[winning_match]
  
  if(add_dir == "vert"){
    loc <- unlist(good_matches[winning_match,1:3])
    xword[(loc[1] - loc[3] + 1) : (loc[1] - loc[3] + length(chars)), loc[2]] <- chars
    word_locs[[tp]] <- cbind((loc[1] - loc[3] + 1) : (loc[1] - loc[3] + length(chars)), loc[2])
  } else {
    loc <- unlist(good_matches[winning_match,1:3])
    xword[loc[1], (loc[2] - loc[3] + 1) : (loc[2] - loc[3] + length(chars))] <- chars
    word_locs[[tp]] <- cbind(loc[1], (loc[2] - loc[3] + 1) : (loc[2] - loc[3] + length(chars)))
  }
  
  used[tp] <- T
  
  return(list(xword = xword, used = used, word_locs = word_locs, success = T))
  
}

# Function to compare two words and count character overlaps (vectorized)
compare_words <- function(word1, word2) {
  
  # Convert words to lowercase and split into characters
  chars1 <- strsplit(tolower(word1), "")[[1]]
  chars2 <- strsplit(tolower(word2), "")[[1]]
  
  # Count occurrences of each character
  table1 <- table(chars1)
  table2 <- table(chars2)
  
  # Find common characters and multiply their counts
  common_chars <- intersect(names(table1), names(table2))
  overlaps <- sum(table1[common_chars] * table2[common_chars])
  
  return(overlaps)
}


# Function to create a matrix of pairwise character overlaps (efficiently)
create_overlap_matrix <- function(words) {
  # Vectorized computation of overlaps for all combinations of words
  overlap_matrix <- outer(words, words, Vectorize(compare_words))
  diag(overlap_matrix) <- 0
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- words
  return(overlap_matrix)
}


#initial generation function
generate_xword <- function(words, max_attempts_word_placement = 20, max_attempts_generate_xword = 5){
  
  total_nchar <- nchar(paste0(words, collapse = ""))
  n_words <- length(words)
  
  # for(debug_seed in 7:100){
  # set.seed(debug_seed)  
  
  generated_xword <- F
  xword_attempt_num <- 1
  while(!generated_xword & xword_attempt_num <= max_attempts_generate_xword){
    
    print(paste0("attempt #", xword_attempt_num, " to generate xword"))
    
    #initalize empty xword objects
    xword <- matrix("", ncol = ceiling(max(sqrt(total_nchar * 10), nchar(words))), 
                    nrow = ceiling(max(sqrt(total_nchar * 10), nchar(words))))
    used <- rep(F, n_words)
    word_locs <- lapply(setNames(words, words), function(x) return(NA))
    
    #place first word
    initialized_xword <- init_xword(xword, words, used, word_locs)
    xword <- initialized_xword[["xword"]]
    used <- initialized_xword[["used"]]
    word_locs <- initialized_xword[["word_locs"]]
    
    #place all the remaining words
    finished_xword <- finish_xword(xword, words, used, word_locs, max_attempts_word_placement)
    xword <- finished_xword[["xword"]]
    used <- finished_xword[["used"]]
    word_locs <- finished_xword[["word_locs"]]
    
    if(all(used)){
      generated_xword <- T
    } else {
      xword_attempt_num <- xword_attempt_num + 1
      print(paste0("FAILED to generate xword"))
    }
    
  }
  
  # }
  
  if(!generated_xword){
    print(paste0("reached max # of attempts (", max_attempts_generate_xword, "), raise <max_attempts_generate_xword> or use different words"))
    char_overlap_mat <- create_overlap_matrix(words)
    match_difficulty <- setNames(apply(diag(1/nchar(words)) %*% char_overlap_mat %*% diag(1/nchar(words)), 1, sum), words)
    match_difficulty <- sort(match_difficulty)
    print(paste0("we suspect the hardest words to match are ", paste0(head(names(match_difficulty), n = 3), collapse = " + ")))
    return(NA)
  }
  
  return(list(xword = xword, word_locs = word_locs))
}

finish_xword <- function(xword, words, used, word_locs, max_attempts_word_placement = 20){
  
  word_attempt_num <- 1
  while(word_attempt_num <= max_attempts_word_placement & !all(used)){
    curr_used <- used
    new_xword <- add_xword(xword, words, used, word_locs)
    xword <- new_xword[["xword"]]
    used <- new_xword[["used"]]
    word_locs <- new_xword[["word_locs"]]
    if(all(curr_used == used)){
      word_attempt_num <- word_attempt_num + 1
    } else {
      cat(paste0("(", words[curr_used != used], ") "))
    }
  }
  
  return(list(xword = xword, 
              words = words, 
              used = used, 
              word_locs = word_locs, 
              success = (word_attempt_num <= max_attempts_word_placement)))
  
}

trim_matrix <- function(xword, word_locs){
  empty_rows <- sapply(1:nrow(xword), function(ri) all(xword[ri,] == ""))
  empty_cols <- sapply(1:ncol(xword), function(ci) all(xword[,ci] == ""))
  new_xword <- xword[!empty_rows, !empty_cols]
  empty_row_inds <- which(empty_rows)
  empty_col_inds <- which(empty_cols)
  new_word_locs <- lapply(word_locs, function(winds){
    if(all(is.na(winds))){
      return(NA)
    } else {
      new_winds <- winds
      new_winds[,1] <- winds[,1] - sum(min(winds[,1]) > empty_row_inds)
      new_winds[,2] <- winds[,2] - sum(min(winds[,2]) > empty_col_inds)
      return(new_winds  )
    }
  })
  return(list(xword = new_xword, word_locs = new_word_locs))
}

#find adjacency matrix of words
xword_adjacency_matrix <- function(word_locs){
  dupmat <- data.frame(word = rep(names(word_locs), sapply(word_locs, nrow)),
                       loc = do.call(rbind, word_locs))
  dupmat$duplicated = duplicated(dupmat[,c("loc.1", "loc.2")]) | 
    duplicated(dupmat[,c("loc.1", "loc.2")], fromLast = T)
  dupmat <- dupmat[dupmat$duplicated,]
  word_intersects <- split(dupmat$word, 
                           apply(dupmat[,c("loc.1", "loc.2")], 1, paste0, collapse = "-"))
  word_intersects <- do.call(rbind, word_intersects)
  adj_mat <- matrix(0, nrow = length(word_locs), ncol = length(word_locs))
  rownames(adj_mat) <- colnames(adj_mat) <- names(word_locs)
  adj_mat[word_intersects] <- adj_mat[word_intersects[,2:1]] <- 1
  return(adj_mat)
}

#find words that can be removed (only one intersection or whose removal leaves alt paths)
removable_words <- function(word_locs){
  word_names <- names(word_locs)
  absent_words <- sapply(word_locs, function(wlx) all(is.na(wlx)))
  
  adj_mat <- xword_adjacency_matrix(word_locs[!absent_words])
  cut_vertices <- igraph::articulation_points(
    igraph::graph_from_adjacency_matrix(
      adj_mat, mode = "undirected", diag = FALSE))
  out <- setdiff(word_names[!absent_words], rownames(adj_mat)[as.numeric(cut_vertices)])
  out <- match(out, word_names)
  return(out)
}

#check word loc validity -- make sure first character of word_locs is in the right spot in xword
valid_word_locs <- function(xword, word_locs){
  all(sapply(1:length(word_locs), function(i){
    all(strsplit(names(word_locs)[i], split = "")[[1]] == xword[word_locs[[i]]])
  }))
}


#swap word location for another
swap_word <- function(xword, words, word_locs){
  
  n_words_to_remove <- sample(length(words):1, 1, F, softmax(1:length(words)/5))
  # n_words_to_remove <- 1

  prop_xword <- xword
  used <- rep(T, length(words))
  for(i in 1:n_words_to_remove){
    word_to_prop <- sample(removable_words(word_locs), 1) #maybe change to prioritize poorly fitting words?
    old_locs <- word_locs[[word_to_prop]]
    other_locs <- do.call(rbind, word_locs[-word_to_prop])
    remove_locs <- !(paste0(old_locs[,1], "-", old_locs[,2]) %in% 
                       paste0(other_locs[,1], "-", other_locs[,2]))
    prop_xword[old_locs[remove_locs,]] <- ""
    prop_word_locs <- word_locs
    prop_word_locs[word_to_prop] <- NA
    used[word_to_prop] <- F
  }
  
  # prop_xword <- add_xword(xword = prop_xword, words = words, used = used, word_locs = prop_word_locs)
  prop_xword <- finish_xword(xword = prop_xword, words = words, used = used, word_locs = prop_word_locs)
  #debugging code
  # valid_word_locs(xword, word_locs)
  # valid_word_locs(prop_xword[["xword"]], prop_xword[["word_locs"]])
  # if(!valid_word_locs(prop_xword[["xword"]], prop_xword[["word_locs"]])){
  #   save(list = c("xword", "words", "word_locs"), file = "~/xword_debug.R", envir = .GlobalEnv)
  #   load("~/xword_debug.R")
  #   break
  # }
  
  
  if(prop_xword$success){
    return(list(xword = prop_xword[["xword"]], word_locs = prop_xword[["word_locs"]]))
  } else {
    return(list(xword = xword, word_locs = word_locs))
  }
}

number_words <- function(word_locs){
  
  #snag info first first letters of each word
  first_cells <- data.frame(do.call(rbind, lapply(word_locs, head, n=1)))
  first_cells$cell_id <- apply(first_cells, 1, paste0, collapse = "-")
  first_cells$wind <- 1:length(word_locs)
  first_cells$word <- names(word_locs)
  distrank_from_origin <- order(first_cells[,1]^2 + first_cells[,2]^2)
  first_cells <- first_cells[distrank_from_origin,]
  first_cells$num <- 1:nrow(first_cells)
  
  #split and order
  split_cells <- split(first_cells, first_cells$cell_id)
  first_cells <- do.call(rbind, lapply(split_cells, function(x){x$num <- min(x$num); x}))
  first_cells$num <- sapply(first_cells$num, function(ni) ni - sum(duplicated(first_cells$num[(ni > first_cells$num)])))
  first_cells$hv <- horiz_vert(word_locs)[first_cells$wind]
  first_cells_nodupes <- first_cells[!duplicated(first_cells$num),]
  
  word_order <- first_cells[order(first_cells$num),]
  word_order <- split(word_order[,c("num", "word")], word_order$hv)
  
  return(list(first_cells = first_cells, word_order = word_order))
}

usr2dev <- function() {
  # Get current par settings
  usr <- par("usr")
  plt <- par("plt")
  dev_size <- dev.size("in")
  pin <- par("pin")
  
  # Calculate scaling factors
  scale_x <- (usr[2] - usr[1]) / (plt[2] - plt[1])
  scale_y <- (usr[4] - usr[3]) / (plt[4] - plt[3])
  
  # Calculate device limits in user coordinates
  x_left <- usr[1] - plt[1] * scale_x
  x_right <- usr[2] + (1 - plt[2]) * scale_x
  y_bottom <- usr[3] - plt[3] * scale_y
  y_top <- usr[4] + (1 - plt[4]) * scale_y
  
  # Return the device coordinates in user units
  return(c(x_left, x_right, y_bottom, y_top))
}



draw_xword <- function(xword, word_locs, write_words = T, 
                       cell_color = "grey99", background_color = "white",
                       print2screen = F, border_lwd = 2, family = "Arial Unicode MS",
                       num_cex = 0.75, border_text_color = "black", cell_scale = 0.9,
                       drop_shadow = T, shadow_length = 0.1){
  
  #initialize plot
  plot(x = NA, y = NA, xaxt = "n", yaxt = "n", 
       xlim = c(0, ncol(xword)+1) * xyrat(), ylim = c(-nrow(xword)-1, 0),
       frame = F, xlab = "", ylab = "", family = family)
  dev_coords <- usr2dev()
  rect(xleft = dev_coords[1],
       xright = dev_coords[2],
       ybottom = dev_coords[3],
       ytop = dev_coords[4],
       col = background_color, xpd = T)
  
  #draw boxes
  occ_cells <- which(xword != "", arr.ind = T)
  if(drop_shadow){
    for(i in 1:nrow(occ_cells)){
      rect(xleft = (occ_cells[i,2] - 1/2 * cell_scale + shadow_length * cell_scale) * xyrat(),
           xright = (occ_cells[i,2] + 1/2 * cell_scale + shadow_length * cell_scale) * xyrat(),
           ybottom = -occ_cells[i,1] - 1/2 * cell_scale - shadow_length * cell_scale,
           ytop = -occ_cells[i,1] + 1/2 * cell_scale - shadow_length * cell_scale,
           col = border_text_color, lwd = border_lwd, border = border_text_color)
    }
  }
  for(i in 1:nrow(occ_cells)){
    rect(xleft = (occ_cells[i,2] - 1/2 * cell_scale) * xyrat(),
         xright = (occ_cells[i,2] + 1/2 * cell_scale) * xyrat(),
         ybottom = -occ_cells[i,1] - 1/2 * cell_scale,
         ytop = -occ_cells[i,1] + 1/2 * cell_scale,
         col = cell_color, lwd = border_lwd, border = border_text_color)
    if(write_words){
      text(x = occ_cells[i,2] * xyrat(), y = -occ_cells[i,1], 
           labels = xword[occ_cells[i,1],occ_cells[i,2]], font = 2, family = family,
           col = border_text_color)
    }
  }
  
  #draw numbers in boxes
  word_numbers <- number_words(word_locs)
  
  for(i in 1:nrow(word_numbers$first_cells)){
    num <- word_numbers$first_cells$num[i]
    num_w <- strwidth(num, cex = num_cex, family = family)
    num_h <- strheight(num, cex = num_cex, family = family)
    text(x = (word_numbers$first_cells[i,2] - 7/16 * cell_scale) * xyrat() + num_w / 2, 
         y = -word_numbers$first_cells[i,1] + 7/16 * cell_scale - num_h / 2, 
         labels = num, cex = num_cex, family = family,
         col = border_text_color)
  }
  
  #print out word order
  if(print2screen){
    cat(paste0("\nhorizontal: ", paste0(word_numbers$word_order$h$num, ". ", word_numbers$word_order$h$word, collapse = ", ")))
    cat(paste0("\n\nvertical: ", paste0(word_numbers$word_order$v$num, ". ", word_numbers$word_order$v$word, collapse = ", ")))
  }
  
}

xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (exp(x) + 1)

fill_adj_elements <- function(xword, word_locs) {
 
  #find indices of all unique filled cells
  char_locs <- do.call(rbind, word_locs)
  char_locs <- char_locs[!is.na(char_locs[,1]),]
  char_locs <- char_locs[!duplicated(char_locs),] + 1
  
  #shift these locations in the four cardinal directions
  xwords_shifted <- lapply(setNames(c("up", "down", "left", "right"), 
                                    c("up", "down", "left", "right")), 
                           function(i) matrix(F, ncol = ncol(xword) + 2, nrow = nrow(xword) + 2))
  xwords_shifted$up[t(t(char_locs) + c(0,-1))] <- T
  xwords_shifted$down[t(t(char_locs) + c(0,1))] <- T
  xwords_shifted$left[t(t(char_locs) + c(-1,0))] <- T
  xwords_shifted$right[t(t(char_locs) + c(1,0))] <- T
  
  #sum these shifted xwords to find cells with >=2 cardinal neighbors
  n_neighbors <- Reduce("+", xwords_shifted)
  spots_to_fill <- n_neighbors >= 2
  spots_to_fill <- spots_to_fill[1:(nrow(spots_to_fill)-1), 1:(ncol(spots_to_fill)-1)]
  new_char_locs <- rbind(char_locs, which(spots_to_fill, arr.ind = T))
  new_char_locs <- new_char_locs[!duplicated(new_char_locs),] - 1
  new_filled_xword <- matrix(F, ncol = ncol(xword), nrow = nrow(xword))
  new_filled_xword[new_char_locs] <- T
  
  #check our work  
  # plot(new_char_locs, pch = 15, col = 2)
  # points(char_locs - 1, pch = 15)
  
  return(list(new_filled_xword = new_filled_xword, char_locs = new_char_locs))
  
}

# Create hash table from matrix
create_hash_table <- function(matrix = NA, coords = NA) {
  h <- hash()
  if(all(is.na(coords))){
    coords <- which(matrix == 1, arr.ind = TRUE)
  }
  for (i in 1:nrow(coords)) {
    key <- paste(coords[i, 1], coords[i, 2], sep = ",")
    .set(h, key, TRUE)
  }
  h
}

# Find nearest distances
find_nearest_manhattan_distances <- function(hash_A, hash_B, max_distance) {
  distances <- numeric(length = length(keys(hash_A)))
  i <- 1

  for (key in keys(hash_A)) {
    coord <- as.numeric(strsplit(key, ",")[[1]])
    found <- FALSE
    for (dist in 0:max_distance) {
      if (found) break
      for (dx in -dist:dist) {
        for (dy in -dist:dist) {
          if (abs(dx) + abs(dy) == dist) { # Only consider perimeter at each distance
            check_key <- paste(coord[1] + dx, coord[2] + dy, sep = ",")
            if (has.key(check_key, hash_B)) {
              distances[i] <- dist
              found <- TRUE
              break
            }
          }
        }
      }
    }
    if (!found) distances[i] <- NA
    i <- i + 1
  }
  distances
}

# Function to find the Euclidean distance
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Find nearest distances using Euclidean distance
find_nearest_euclidean_distances <- function(hash_A, hash_B, max_distance) {
  distances <- numeric(length = length(keys(hash_A)))
  i <- 1
  
  for (key in keys(hash_A)) {
    coord <- as.numeric(strsplit(key, ",")[[1]])
    found <- FALSE
    dist <- 0
    while (!found && dist <= max_distance) {
      for (dx in -dist:dist) {
        for (dy in -dist:dist) {
          if (euclidean_distance(0, 0, dx, dy) <= dist) { # Check if the point is within the circle
            check_key <- paste(coord[1] + dx, coord[2] + dy, sep = ",")
            if (has.key(check_key, hash_B)) {
              distances[i] <- euclidean_distance(coord[1], coord[2], coord[1] + dx, coord[2] + dy)
              found <- TRUE
              break
            }
          }
        }
        if (found) break
      }
      dist <- dist + 1
    }
    if (!found) distances[i] <- NA
    i <- i + 1
  }
  distances
}

evaluate_picture_concordance <- function(char_locs, pic_locs, use_manhattan_dist = F){
  
  #find useful info about these two tables
  dim_pic <- apply(pic_locs, 2, range)
  pic_locs <- t(t(pic_locs) - dim_pic[1,]) + 1
  dim_pic <- t(t(dim_pic) - dim_pic[1,]) + 1
  dim_char <- apply(char_locs, 2, range)
  
  # pic_density <- nrow(pic_locs) / prod(dim_pic[2,])
  # xword_density <- nrow(new_char_locs) / prod(dim_xword[2,])
  
  #square both xword and picture, ie expand dims so they both 
  #occupy a square shape corresponding to the larger dimension
  ldim_pic <- which.max(dim_pic[2,])
  pic_buff <- c(ifelse(ldim_pic == 1, 0, ceiling(abs(diff(dim_pic[2,])) / 2)),
                ifelse(ldim_pic == 1, ceiling(abs(diff(dim_pic[2,])) / 2), 0))
  
  ldim_char <- which.max(dim_char[2,])
  char_buff <- c(ifelse(ldim_char == 1, 0, ceiling(abs(diff(dim_char[2,])) / 2)),
                ifelse(ldim_char == 1, ceiling(abs(diff(dim_char[2,])) / 2), 0))
  
  char_locs <- t(t(char_locs) + char_buff)
  pic_locs <- t(t(pic_locs) + pic_buff)
  
  #now that they are both inside squares, resize the larger one to the dimension of the smaller one
  #almost always this'll be resizing the pic but idk should check anyway
  pic_larger <- dim_pic[2,ldim_pic] > dim_char[2,ldim_char]
  if(pic_larger){
    pic_locs <- pic_locs * dim_char[2,ldim_char] / dim_pic[2,ldim_pic]
    pic_locs <- ceiling(pic_locs)
    pic_locs <- pic_locs[!duplicated(pic_locs),]
  } else {
    char_locs <- char_locs * dim_pic[2,ldim_pic] / dim_char[2,ldim_char]
    char_locs <- ceiling(char_locs)
    char_locs <- char_locs[!duplicated(char_locs),]
  }
  
  #now that these are the same dim, we can compare them
  #creating minimum hash tables and finding the manhattan distance
  xlims <- range(c(pic_locs[,1], char_locs[,1]))
  ylims <- range(c(pic_locs[,2], char_locs[,2]))
  pic_not_in_char <- pic_locs[!duplicated(rbind(char_locs, pic_locs))[(nrow(char_locs)+1):(nrow(char_locs) + nrow(pic_locs))],]
  char_not_in_pic <- char_locs[!duplicated(rbind(pic_locs, char_locs))[(nrow(pic_locs)+1):(nrow(pic_locs) + nrow(char_locs))],]
  
  # plot(char_locs, pch = 15, xlim = xlims, ylim = ylims, col = adjustcolor(1, 0.2))
  # points(pic_locs, pch = 15, col = adjustcolor(2, 0.2))
  # points(pic_not_in_char, col = adjustcolor(3, 1))
  # points(char_not_in_pic, col = adjustcolor(4, 1))
  
  
  hash_char <- create_hash_table(coords = char_locs)
  hash_pic <- create_hash_table(coords = pic_locs)
  
  if(nrow(char_not_in_pic) > 0){
    hash_char_not_in_pic <- create_hash_table(coords = char_not_in_pic)  
    if(use_manhattan_dist){
      distances_char_2_pic <- find_nearest_manhattan_distances(hash_char_not_in_pic, hash_pic, 
                                                               max_distance = diff(xlims) + diff(ylims))
    } else {
      distances_char_2_pic <- find_nearest_euclidean_distances(hash_char_not_in_pic, hash_pic, 
                                                               max_distance =   sqrt(diff(xlims)^2 + diff(ylims)^2))      
    }
  } else {
    distances_char_2_pic <- 0
  }
  
  if(nrow(pic_not_in_char) > 0){
    hash_pic_not_in_char <- create_hash_table(coords = pic_not_in_char)
    if(use_manhattan_dist){
      distances_pic_2_char <- find_nearest_manhattan_distances(hash_pic_not_in_char, hash_char, 
                                                               max_distance = diff(xlims) + diff(ylims))  
    } else {
      distances_pic_2_char <- find_nearest_euclidean_distances(hash_pic_not_in_char, hash_char, 
                                                               max_distance =   sqrt(diff(xlims)^2 + diff(ylims)^2))
    }
  } else {
    distances_pic_2_char <- 0
  }
  
  return(list(distances_char_2_pic = distances_char_2_pic,
              distances_pic_2_char = distances_pic_2_char))
  
}

#### run MCMC ####


#load in a picture whose shape we want to target
pic <- png::readPNG("~/Downloads/black_heart.png")
has_transparency <- min((c(1,0) - mean(pic[,,4] > 0.5)) * c(1,-1)) > c(0.05)
if(has_transparency){
  pic_mat <- pic[,,4] > 0.01
  pic_locs <- which(pic_mat, arr.ind = T)
} else {
  pic_mat <- (pic[,,1] + pic[,,2] + pic[,,3]) <= 1
  pic_locs <- which(pic_mat, arr.ind = T)
}

#set MCMC params
n_iter <- 100

# generate and draw an initial xword
raw_xword <- generate_xword(words, max_attempts_word_placement = 20)
xword <- raw_xword[["xword"]]
word_locs <- raw_xword[["word_locs"]]
curr_score <- score_xword(xword, words, word_locs)
trimmed_xword_data <- trim_matrix(xword, word_locs)

ptsize_x <- 5
font_family_for_nums = "Home Christmas"
# font_family_for_nums = "Arial Unicode MS"
png("~/xword.png", family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
    height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
par(xpd = NA, mar = c(0,0,0,0))
draw_xword(trimmed_xword_data$xword, trimmed_xword_data$word_locs, border_lwd = 2 * ptsize_x,
           family = font_family_for_nums, num_cex = 1)
dev.off()


#run MCMC
for(i in 1:n_iter){
  cat(paste0("(i: ", i, ", "))
  prop_raw_xword <- swap_word(xword, words, word_locs)
  prop_score <- score_xword(prop_raw_xword[["xword"]], words, prop_raw_xword[["word_locs"]])
  # prop_score <- prop_score +
  #   ifelse(sum(find_n_extra(prop_raw_xword[["xword"]], words)) == 0,
  #      0,
  #      Inf)
  accept_prob <- invlogit(curr_score - prop_score)
  cat(paste0("\ndNC", ": ", sum(prop_raw_xword$xword != "") - sum(xword != ""),", "))
  if(rbinom(1, 1, prob = accept_prob) == 1){
    
    xword <- prop_raw_xword[["xword"]]
    word_locs <- prop_raw_xword[["word_locs"]]
    curr_score <- prop_score
    cat(paste0("dS", ": ", curr_score))
    
    #draw if desired
    trimmed_xword <- trim_matrix(xword, word_locs)$xword
    trimmed_word_locs <- trim_matrix(xword, word_locs)$word_locs
    png("~/xword.png", family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
        height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
    par(xpd = NA, mar = c(0,0,0,0))
    draw_xword(trimmed_xword, trimmed_word_locs, border_lwd = 2 * ptsize_x,
               family = "Arial Unicode MS", num_cex = 1)
    dev.off()
    
  } else {
    cat(paste0("dS", ": ", "^"))
  }
  cat(")\n")

}


#function for running MCMC
run_MCMC <- function(xword, word_locs, n_iter = 100, print_progress = F){

  words <- names(word_locs)
  curr_score <- score_xword(xword, words, word_locs)
  for(i in 1:n_iter){
    if(print_progress){cat(paste0("(i: ", i, ", "))}
    
    prop_raw_xword <- swap_word(xword, words, word_locs)
    prop_score <- score_xword(prop_raw_xword[["xword"]], words, prop_raw_xword[["word_locs"]])
    accept_prob <- 1 - invlogit(curr_score - prop_score)
    
    if(print_progress){cat(paste0("\ndNC", ": ", sum(prop_raw_xword$xword != "") - sum(xword != ""),", "))}
    
    #accept-reject
    if(rbinom(1, 1, prob = accept_prob) == 1){
      
      xword <- prop_raw_xword[["xword"]]
      word_locs <- prop_raw_xword[["word_locs"]]
      curr_score <- prop_score
      if(print_progress){cat(paste0("dS", ": ", curr_score))}
    
    } else {
      if(print_progress){cat(paste0("dS", ": ", "^"))}
    }
    if(print_progress){cat(")\n")}
  }
  
  return(list(xword = xword, word_locs = word_locs))
}

#generate random xwords
xword_dir <- "~/xwords/russian_runs/"
if(!dir.exists(xword_dir)) dir.create(xword_dir)
ptsize_x <- 5
xword_indices <- 101:110
for(i in xword_indices){
  print(i)
  
  raw_xword <- generate_xword(words, max_attempts_word_placement = 20)
  
  #save original and also refined version (MCMC output)
  for(j in 1:2){
    if(j == 2){
      refined_xword <- run_MCMC(raw_xword[["xword"]], raw_xword[["word_locs"]], n_iter = 50)
    } else {
      refined_xword <- raw_xword
    }
    
    save(refined_xword, file = paste0(xword_dir, "xword_", i, ifelse(j==1, "-orig", ""), ".RData"))
    trimmed_xword_data <- trim_matrix(refined_xword[["xword"]], refined_xword[["word_locs"]])
    
    #draw xword with words filled in
    png(paste0(xword_dir, "xword_", i, ifelse(j==1, "-orig", ""), ".png"), family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
        height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
    par(xpd = NA, mar = c(0,0,0,0))
    draw_xword(trimmed_xword_data$xword, trimmed_xword_data$word_locs, border_lwd = 2 * ptsize_x,
               family = "Arial Unicode MS", num_cex = 1)
    dev.off()
    
    #draw blank xword
    png(paste0(xword_dir, "xword_", i, ifelse(j==1, "-orig", ""), "_blank.png"), family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
        height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
    par(xpd = NA, mar = c(0,0,0,0))
    draw_xword(trimmed_xword_data$xword, trimmed_xword_data$word_locs, border_lwd = 2 * ptsize_x, write_words = F,
               family = font_family_for_nums, num_cex = 1)
    dev.off()
    
    #draw shape of xword
    png(paste0(xword_dir, "xword_", i, ifelse(j==1, "-orig", ""), "_shape.png"), family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
        height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
    par(xpd = NA, mar = c(0,0,0,0))
    draw_xword(trimmed_xword_data$xword, trimmed_xword_data$word_locs, border_lwd = 2 * ptsize_x, write_words = F, 
               cell_color = "black", 
               background_color = "white",
               family = font_family_for_nums, num_cex = 1, drop_shadow = F, cell_scale = 1)
    dev.off()
    
    #save words and clues
    word_numbers <- number_words(trimmed_xword_data$word_locs)
    
    sink(paste0(xword_dir, "xword_", i, ifelse(j==1, "-orig", ""), "_answers-clues.txt"))
    
    cat("ANSWERS:\n\n")
    cat(paste0("\nHorizontal: ", paste0(word_numbers$word_order$h$num, ". ", full_words[word_numbers$word_order$h$word], collapse = " ")))
    cat(paste0("\n\nVertical: ", paste0(word_numbers$word_order$v$num, ". ", full_words[word_numbers$word_order$v$word], collapse = " ")))
    
    cat("\n\nCLUES:\n\n")
    cat(paste0("\nHorizontal: ", paste0(word_numbers$word_order$h$num, ". ", clues[word_numbers$word_order$h$word], collapse = " ")))
    cat(paste0("\n\nVertical: ", paste0(word_numbers$word_order$v$num, ". ", clues[word_numbers$word_order$v$word], collapse = " ")))
    
    sink()
    
  }
  
}

#### further refinement ####
# iterate over favorites for further refinements
xword_dir <- "~/xwords/fave-iterations/"
if(!dir.exists(xword_dir)) dir.create(xword_dir)
xword_indices <- 1:50
fave_i <- 244
ptsize_x <- 5
for(i in xword_indices){
  
  load(paste0("~/xwords/runs/xword_", fave_i, ".RData"))
  
  if(i != 1){
    refined_xword <- run_MCMC(refined_xword[["xword"]], refined_xword[["word_locs"]], n_iter = 10)  
  }
  save(refined_xword, file = paste0(xword_dir, "xword_", fave_i, "-", i, ".RData"))
  trimmed_xword_data <- trim_matrix(refined_xword[["xword"]], refined_xword[["word_locs"]])
  
  #draw xword with words filled in
  png(paste0(xword_dir, "xword_", fave_i, "-", i, ".png"), family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
      height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
  par(xpd = NA, mar = c(0,0,0,0))
  draw_xword(trimmed_xword_data$xword, trimmed_xword_data$word_locs, border_lwd = 2 * ptsize_x,
             family = font_family_for_nums, num_cex = 1)
  dev.off()
  
  #draw blank xword
  png(paste0(xword_dir, "xword_", fave_i, "-", i, "_blank.png"), family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
      height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
  par(xpd = NA, mar = c(0,0,0,0))
  draw_xword(trimmed_xword_data$xword, trimmed_xword_data$word_locs, border_lwd = 2 * ptsize_x, write_words = F,
             family = font_family_for_nums, num_cex = 1)
  dev.off()
  
  #draw shape of xword
  png(paste0(xword_dir, "xword_", fave_i, "-", i, "_shape.png"), family = font_family_for_nums, width = ncol(trimmed_xword_data$xword) * 30 * ptsize_x, 
      height = nrow(trimmed_xword_data$xword) * 30 * ptsize_x, pointsize = 12 * ptsize_x)
  par(xpd = NA, mar = c(0,0,0,0))
  draw_xword(trimmed_xword_data$xword, trimmed_xword_data$word_locs, border_lwd = 2 * ptsize_x, write_words = F, 
             cell_color = "black", 
             background_color = "white",
             family = font_family_for_nums, num_cex = 1)
  dev.off()
  
  #save words and clues
  word_numbers <- number_words(trimmed_xword_data$word_locs)
  
  sink(paste0(xword_dir, "xword_", fave_i, "-", i, "_answers-clues.txt"))
  
  cat("ANSWERS:\n\n")
  cat(paste0("\nHorizontal: ", paste0(word_numbers$word_order$h$num, ". ", full_words[word_numbers$word_order$h$word], collapse = " ")))
  cat(paste0("\n\nVertical: ", paste0(word_numbers$word_order$v$num, ". ", full_words[word_numbers$word_order$v$word], collapse = " ")))
  
  cat("\n\nCLUES:\n\n")
  cat(paste0("\nHorizontal: ", paste0(word_numbers$word_order$h$num, ". ", clues[word_numbers$word_order$h$word], collapse = " ")))
  cat(paste0("\n\nVertical: ", paste0(word_numbers$word_order$v$num, ". ", clues[word_numbers$word_order$v$word], collapse = " ")))
  
  sink()
  
}