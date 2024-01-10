
remove_bottom <- function(x, replacement){
  nobot <- gsub("g|j|p|q|y|,|\\(|\\)|Q", replacement, x)
  nobot <- gsub("\\_s*\\{[^\\)]+\\}", replacement, nobot) #underscore in brackets
  nobot <- gsub("_[a-z|0-9|A-Z]{1}", replacement, nobot) #underscore w/ one letter following
  nobot
}

remove_top <- function(x, replacement){
  notop <- gsub("\\^s*\\{[^\\)]+\\}", replacement, x)
  notop <- gsub("\\^[a-z|0-9|A-Z]{1}", replacement, notop)
  notop
}

remove_tb <- function(x, replacement){
  remove_top(remove_bottom(x, replacement), replacement)
}

split_and_return_vec <- function(string, keyword){
  s <- paste0(string, " ")
  pattern <- sprintf("\\b%s\\b", keyword)
  segments <- unlist(strsplit(s, pattern, perl = TRUE))
  segments[length(segments)] <- substr(segments[length(segments)], 1, nchar(segments[length(segments)]) - 1)
  if(length(segments) > 1){
    segs <- rep("", 2 * length(segments) - 1)
    segs[2*(1:length(segments))-1] <- segments
    segs[2*(1:(length(segments)-1))] <- keyword
    return(segs)
  } else {
    return(segments)  
  }
  
}


wrap_in_dollarsigns <- function(string, exclude) {
  
  # Prepare regex pattern for each exclusion
  pattern_list <- lapply(exclude, function(keyword) {
    if (keyword %in% c("(", ")")) {
      # Directly escape parentheses for regex
      return(paste0("\\", keyword))
    } else {
      # Use word boundaries for other keywords
      return(paste0("\\b", keyword, "\\b"))
    }
  })
  
  # Create a combined regex pattern
  combined_pattern <- paste(pattern_list, collapse = "|")
  
  # Split the string on the pattern while retaining the separators
  segments <- unlist(strsplit(string, combined_pattern, perl = TRUE))
  separators <- regmatches(string, gregexpr(combined_pattern, string, perl = TRUE))[[1]]
  
  # #separate segments from excluded words
  # segments <- string
  # for(keyword in exclude){
  #   segments <- unlist(sapply(segments, function(s) split_and_return_vec(s, keyword)))
  # }
  
  if(length(segments) > 1){
    
    # Process segments and separators
    processed_segments <- vector("list", length = length(segments) + length(separators))
    processed_segments[c(TRUE, FALSE)] <- segments
    processed_segments[c(FALSE, TRUE)] <- separators  
    segments <- unlist(processed_segments)
    
    excluded_segments <- segments != "" & (segments %in% exclude)
    segments[excluded_segments] <- sapply(which(excluded_segments), function(ei){
      if(ei != 1 & ei != length(excluded_segments)){
        return(paste0(" ", segments[ei], " "))
      } else if (ei != 1){
        return(paste0(" ", segments[ei]))
      } else {
        return(paste0(segments[ei], " "))
      }
    })
    
    target_segments <- segments != "" & !excluded_segments
    segments[target_segments] <- paste0("$", gsub(" ", "~", trimws(segments[target_segments])), "$")
    return(paste0(segments, collapse = ""))
    
  } else {
    
    #swap in ~ for spaces
    segments <- gsub(" ", "~", trimws(segments))
    return(paste0("$", segments, "$", collapse = "")) 
  }
  
}

xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}



text3 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, col = 1, replacement = "a",
                  return_bb = F, return_verts = F, dline = F, vexp = T, hexp = F, srt = 0, ...){
  
  
  #convert text label to expression
  if(all(class(labels) == "character")){
    if (!grepl("^\\$.+\\$$", labels)) {
      
      # Enclose in dollar signs if there aren't any (hack to get commas etc. to render well)
      excluded_words <- c("in", "for", "'", "else", "(", ")", "anyone")
      if(any(labels %in% excluded_words)){
        
        word_expression <- labels
        
      } else {
        
        wrapped_labels <- wrap_in_dollarsigns(labels, excluded_words)
        word_expression <- correct_l2xp_vec(x = wrapped_labels)
        
      }
      
    }
  } else {
    word_expression <- labels
  }
  
  #use the plain label for the expression?
  if(!vexp){
    word_expression <- labels
  }
  
  #find general params
  strh <- strheight(word_expression, cex = cex, ...)
  base_strh <- strheight(correct_l2xp_vec("GIs"), cex = cex, ...)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strh / 2, 0)
  
  #adjust in case of ledding
  nobot <- remove_bottom(labels, replacement)
  ebot <- strheight(latex2exp::TeX(nobot), cex = cex) - strh
  
  notop <- remove_top(labels, replacement)
  etop <- strheight(latex2exp::TeX(notop), cex = cex) - strh
  
  nobottop <- remove_tb(labels, replacement)
  ebottop <- strheight(latex2exp::TeX(nobottop), cex = cex) - strh
  
  #ugh this was obnoxious to figure out
  ebt_delta <- ebottop - (ebot + etop)
  adj_ledding <- ifelse(abs(ebt_delta) > 1E-6, 
                        ebot / 2 - (ebottop - ebot) / 2, 
                        ebot / 2 - etop / 2)
  adj_ledding <- adj_ledding - ifelse(base_strh > strh, (base_strh - strh) / 2, 0)
  
  #use the plain label for the expression?
  if(!hexp){
    word_expression <- labels
  }
  
  #find width param (unaffected by leading issues)
  strw <- strwidth(word_expression, cex = cex, ...)
  
  #adjust base location
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strw / 2, 0)
  
  #find  bounding box in the unrotated case
  unrotated_bounding_box <- list(xleft = adj_x - strw / 2, 
                                 xright = adj_x + strw / 2, 
                                 ybottom = adj_y - strh / 2 + adj_ledding, 
                                 ytop = adj_y + strh / 2 + adj_ledding)
  
  
  #compute offsets for rotation
  theta <- 2 * pi * srt / 360
  if(abs(theta %% (2*pi)) < 1E-6){
    adj_x_rot <- 0
    adj_y_rot <- 0  
  }
  
  #find axis of rotation
  axis_rot <- box_midpoint <- c(x = (unrotated_bounding_box[["xleft"]] + unrotated_bounding_box[["xright"]]) / 2,
                y = (unrotated_bounding_box[["ybottom"]] + unrotated_bounding_box[["ytop"]]) / 2)
  
  if(any(pos %in% c(1,3))){
    if(any(pos == 1)){
      axis_rot["y"] <- unrotated_bounding_box[["ytop"]]
    }
    if(any(pos == 3)){
      axis_rot["y"] <- unrotated_bounding_box[["ybottom"]]
    }
  }
  if(any(pos %in% c(2,4))){
    if(any(pos == 2)){
      axis_rot["x"] <- unrotated_bounding_box[["xright"]]
    }
    if(any(pos == 4)){
      axis_rot["x"] <- unrotated_bounding_box[["xleft"]]
    }
  }
  
  rotated_adj <- rotate_points(pts = box_midpoint, 
                               axis_rot = axis_rot,
                               theta = theta) - box_midpoint
  adj_x_rot <- rotated_adj[,1]
  adj_y_rot <- rotated_adj[,2]
  verts <- rbind(tr = c(unrotated_bounding_box[["xright"]],unrotated_bounding_box[["ytop"]]),
                 tl = c(unrotated_bounding_box[["xleft"]],unrotated_bounding_box[["ytop"]]),
                 bl = c(unrotated_bounding_box[["xleft"]],unrotated_bounding_box[["ybottom"]]),
                 br = c(unrotated_bounding_box[["xright"]],unrotated_bounding_box[["ybottom"]]))
  verts <- rotate_points(pts = verts, 
                         axis_rot = axis_rot,
                         theta = theta)
  
  if(return_bb){ #just return the box that bounds the text without plotting
    return(bounding_box)
  } else if(return_verts){
    return(verts)
  } else {
    #print the text itself
    text(x = adj_x + adj_x_rot, 
         y = adj_y + adj_ledding + adj_y_rot, 
         labels = word_expression, pos = NULL, cex = cex, col = col, srt = srt, ...)
    
    #draw a box around it if desired
    if(drect){
      polygon(x = verts[,1], y = verts[,2], border = col)
    }
    if(dline){
      abline(h=y - strheight(latex2exp::TeX("GIs"), cex = cex) / 2, lwd = 0.5)  
    }
  }
  
}


rad2rot <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)

rotate_points <- function(pts, axis_rot, theta) {
  
  #convert to matrix format (for lin alg operations)
  if(!("matrix" %in% class(pts))){
    pts <- matrix(pts, ncol = 2)
  }
  
  #find useful variables
  rotmat <- rad2rot(theta)
  xyr <- xyrat()
  
  #center and rotate
  pts_eq_axis <- t(t(pts) - axis_rot)
  pts_eq_axis[,2] <- pts_eq_axis[,2] * xyr
  rot_centered_points <- t(rotmat %*% t(as.matrix(pts_eq_axis)))
  
  #transform back to orig space and uncenter
  rot_centered_points[,2] <- rot_centered_points[,2] / xyr
  rot_points <- t(t(rot_centered_points) + axis_rot)
  
  #return result
  colnames(rot_points) <- c("x", "y")
  return(rot_points)
}


plot(1,1,xlim = c(-2,2), ylim = c(-30,30), col = "white", xlab = "", ylab = "")
axrot <- c(x = -1/2, y = 15)
ptval <- c(x = -1, y = 5)
points(ptval["x"], ptval["y"], col = 1, pch = 19)
points(axrot["x"], axrot["y"], col = 1, pch = 8)
points(rotate_points(ptval, axrot, pi/3), col = 2, pch = 19)

plot(1,1,xlim = c(-2,2), ylim = c(-30,30), col = "white", xlab = "", ylab = "")
text3(1, 1, pos = c(3,4), labels = "testing", cex = 2, drect = T)
text3(1, 1, pos = c(3,4), labels = "testing", cex = 2, srt = 15, col = 2, drect = T)
text3(1, 1, pos = c(3,4), labels = "testing", cex = 2, srt = 135, col = 3, drect = T)

#####