txt <- c("    file.rename(from = paste0(\"~/repos/ldsc/custom_genesets/cluster_\",", " cli, \".(chr)_\",", " cri, \".annot.gz\"), ")

extract_paste0 <- function(txt){
  
  #identify location of script keyword? 
  #identify location of paste0( string
  string <- paste0(txt, collapse = "")
  paste0_locs <- unlist(gregexpr("paste0", string))
  
  #check we have paste0s in here
  if(all(paste0_locs == -1)){
    return("")
  }
  
  #process last paste0, then swap it in w/ a placeholder and call this function again
  last_paste0 <- tail(paste0_locs, 1)
  open_paren_locs <- unlist(gregexpr("\\(", string))
  first_open_paren <- open_paren_locs[min(which(open_paren_locs > last_paste0))]
  
  #check we have an open paren in here
  if(all(is.na(open_paren_locs))){
    return("")
  }
  
  #now find the corresponding close paren
  #initialize counter to 1, 
  #count next open paren as +1
  #count next closing paren as -1
  #when counter reaches 0, extract contents between parens
  close_paren_locs <- unlist(gregexpr(")", string))
  
  if(all(close_paren_locs == -1) | all(close_paren_locs <= first_open_paren)){
    return("")
  }
  
  paren_pos <- rep(0, tail(close_paren_locs, 1) - first_open_paren)
  paren_pos[(open_paren_locs - first_open_paren + 1)[sign(open_paren_locs - first_open_paren + 1) == 1]] <- 1
  paren_pos[(close_paren_locs - first_open_paren + 1)[sign(close_paren_locs - first_open_paren + 1) == 1]] <- -1
  corresp_close_paren_disp <- min(which(cumsum(paren_pos) == 0))
  corresp_close_paren <- first_open_paren + corresp_close_paren_disp - 1
  
  if(all(is.na(corresp_close_paren))){
    return("")
  }
  
  #find commas and quotes
  #then find all commas not inside of quotes
  target_substr <- substr(string, first_open_paren, corresp_close_paren)
  quote_locs <- sort(c(unlist(gregexpr("\"", target_substr)), unlist(gregexpr("\'", target_substr))))
  quote_locs <- quote_locs[quote_locs > 0]
  quote_locs <- t(matrix(quote_locs, nrow = 2))
  quoted_content <- apply(quote_locs, 1, function(x) substr(target_substr, x[1] + 1, x[2] - 1))
  # comma_locs <- unlist(gregexpr(",", target_substr))
  
  #actually let's keep this simple
  #and just return quoted content separated by *
  last_paste0_contents <- paste0(quoted_content, collapse = "*")
  
  #return result or call function again
  if(length(paste0_locs) == 1){
    return(last_paste0_contents)
  } else {
    next_outer_string <- c(substr(string, 1, first_open_paren - 7), 
                           paste0("\"", last_paste0_contents, "\""),
                           substr(string, corresp_close_paren + 1, nchar(string)))
    return(extract_paste0(next_outer_string))
  }
  
}

extract_paste0(txt)
extract_paste0(c("abc", "xyz", "paste0(\"y_\",  foo, \"_bar_\",", "paste0(\"nested_\", i), \".txt\")", "asdf"))
