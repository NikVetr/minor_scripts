# Install and load necessary packages
library(extrafont)

#### functions ####

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




text_cols <- function(string, cols, x, y, cex = 1, ...){
  for(char_i in 1:nchar(string)){
    txt_exploded <- c(substr(string, 1, char_i-1), substr(string, char_i, char_i), substr(string, char_i+1, nchar(string)))
    text2(x = x, y = y, #text3(x = x, y = y,  
          labels = bquote(phantom(.(txt_exploded[1])) * .(txt_exploded[2]) * phantom(.(txt_exploded[3]))), 
          col = cols[char_i], cex = cex, ...)
  }
}

find_optimal_cex_and_lines <- function(txt, rect_coords, rect_rescaling_ratio = 1, str_height_rescaling_ratio = 1, fixed_cex = NA,
                                       newlines = NA){
  
  strwidths <- strwidth(correct_l2xp_vec(txt))
  strheight <- (strheight("M\nM") - strheight("M")) * str_height_rescaling_ratio
  space_width_min <- strwidth(" ")
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectwidth <- abs(rect_coords$x1 - rect_coords$x0) * rect_rescaling_ratio
  rectheight <- abs(rect_coords$y1 - rect_coords$y0) * rect_rescaling_ratio
  
  if(all(is.na(newlines))){
    newlines <- rep(F, length(txt))
  }
  
  # ceiling(cumsum(strwidths) / rectwidth)
  data <- list(strwidths = strwidths, strheight = strheight, space_width_min = space_width_min, rectwidth = rectwidth, rectheight = rectheight,
               newlines = newlines)
  par <- log((rectwidth * rectheight) / ((sum(strwidths) + space_width_min * length(strwidths)) * strheight) * 0.5) #initialize cex
  while(compute_space_remaining(data, par) == Inf){
    par <- par + log(0.95)
  }
  # plot(1:120/100, sapply(log(1:120/100), function(cex) compute_space_remaining(data, cex)), type = "l")
  if(is.na(fixed_cex)){
    opt_cex <- suppressWarnings(optimx::optimx(par = par, fn = compute_space_remaining, data = data, hess = NULL, 
                                               method = c('nlminb', 'nlm', 'BFGS', 'Nelder-Mead')[4], hessian=FALSE, #can't compute hessian bc of sharp jumps when new line is formed? or maybe not?
                                               control = list(maxit = 1E4, trace = 0, kkt=FALSE)))
    
    #check if it worked
    compute_space_remaining(data = data, par = opt_cex$p1)
    compute_space_remaining(data = data, par = opt_cex$p1 + 0.0001)
    
    return(list(cex = exp(opt_cex$p1), 
                words_on_lines = put_words_on_lines(data = data, par = opt_cex$p1),
                space_width_min = space_width_min * exp(opt_cex$p1),
                vertical_space = strheight * exp(opt_cex$p1),
                vertical_space_noWS = strheight("M") * exp(opt_cex$p1))
    )
  } else {
    return(list(cex = fixed_cex, 
                words_on_lines = put_words_on_lines(data = data, par = fixed_cex),
                space_width_min = space_width_min * fixed_cex,
                vertical_space = strheight * fixed_cex,
                vertical_space_noWS = strheight("M") * fixed_cex)
    )
  }
  
  
}

compute_space_remaining <- function(data, par){
  put_words_on_lines(data, par, return_space_remaining = T)
}

put_words_on_lines <- function(data, par, return_space_remaining = F){
  
  #clarify par-cex relationship
  cex <- exp(par)
  
  #unwrap data
  strwidths_mod <- data$strwidths * cex
  strheight_mod <- data$strheight * cex
  space_width_min_mod <- data$space_width_min * cex
  rectwidth_mod <- data$rectwidth
  rectheight_mod <- data$rectheight
  newlines <- data$newlines
  
  #check that no words are wider than a line
  if(any(strwidths_mod > rectwidth_mod) & return_space_remaining){
    return(Inf)
  }
  
  txi <- 1
  linei <- 1
  current_width <- strwidths_mod[txi] + space_width_min_mod #if first line has a line break, need to clear current_width
  txt_lines <- list()
  txt_lines[[1]] <- txi
  txt_lines[[2]] <- integer(0)
  
  #check if first word has a line break
  linei <- linei + newlines[txi]
  current_width <- ifelse(newlines[txi], 0, current_width)
  
  while(txi < length(txt)){
    
    txi <- txi + 1
    
    prop_current_width <- current_width + strwidths_mod[txi]
    add_to_line <- (prop_current_width < rectwidth_mod) | (length(txt_lines[[linei]]) == 0)
    
    if(add_to_line){
      txt_lines[[linei]] <- c(txt_lines[[linei]], txi)
      current_width <- prop_current_width + space_width_min_mod
    } else {
      current_width <- strwidths_mod[txi] + space_width_min_mod
      linei <- linei + 1
      txt_lines[[linei]] <- txi
    }
    
    if(newlines[txi]){
      linei <- linei + 1
      current_width <- 0
      txt_lines[[linei]] <- integer(0)
    }
    
    #print for debugging
    # print(paste0("line: ", linei, ", txt_i: ", txi, ", txt: ", txt[txi], ", nl: ", newlines[txi], ", curr_w: ", round(current_width, 4)))
    
  }
  
  if(return_space_remaining){
    last_line_width_remaining <- rectwidth_mod - sum(strwidths_mod[txt_lines[[linei]]])
    current_height <- linei * strheight_mod
    vspace_remaining <- (rectheight_mod - current_height) * rectwidth_mod
    hspace_remaining <- (last_line_width_remaining * strheight_mod)
    
    if(vspace_remaining < 0 | hspace_remaining < 0){
      return(Inf)
    } else {
      total_space_remaining <- vspace_remaining + hspace_remaining
      return(total_space_remaining)
    }
  } else {
    return(txt_lines)
  }
  
}

text_wrapped_words <- function(txt, rect_coords, optimal_word_placement_inf, justified = F, str_height_lower_start_ratio = 0, 
                               str_width_lefter_start_ratio = 0, rect_rescaling_ratio = 1, col = "black", multicolor_words = F, cols_list,
                               vertically_justified = F,...){
  
  ws_height <- optimal_word_placement_inf$vertical_space - optimal_word_placement_inf$vertical_space_noWS
  cex <- optimal_word_placement_inf$cex
  curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
  curr_y <- rect_coords$y0 - optimal_word_placement_inf$vertical_space_noWS * (1 + str_height_lower_start_ratio) / 2
  nlines <- length(optimal_word_placement_inf$words_on_lines)
  strwidths_plotting <- strwidth(correct_l2xp_vec(txt), cex = cex)
  
  space_left_on_lines <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    abs(rect_coords$x0 - rect_coords$x1) * rect_rescaling_ratio - sum(strwidths_plotting[optimal_word_placement_inf$words_on_lines[[linei]]]))
  justified_space_between_words <- sapply(1:length(optimal_word_placement_inf$words_on_lines), function(linei)
    space_left_on_lines[linei] / (length(optimal_word_placement_inf$words_on_lines[[linei]]) - 1))
  
  if(vertically_justified){
    space_taken_vertically <- optimal_word_placement_inf$vertical_space_noWS * nlines + 
      (optimal_word_placement_inf$vertical_space - optimal_word_placement_inf$vertical_space_noWS) * (nlines - 1)
    space_left_vertically <- abs(rect_coords$y0 - rect_coords$y1) - space_taken_vertically
    extra_leading <- space_left_vertically / (nlines - 1)
  } else {
    extra_leading <- 0
  }
  
  words_written <- 0
  for(linei in 1:nlines){
    for(wordi in 1:length(optimal_word_placement_inf$words_on_lines[[linei]])){
      
      words_written <- words_written + 1
      word_to_write <- txt[optimal_word_placement_inf$words_on_lines[[linei]][wordi]]
      
      #adjust for sticky-outy bits
      word_expression <- correct_l2xp_vec(word_to_write)
      
      if(multicolor_words){
        text_cols(x = curr_x, y = curr_y, cex = cex,
                  string = word_expression, pos = 4, cols = cols_list[[words_written]])
      } else {
        # text2(x = curr_x, y = curr_y, cex = cex,
        #       labels = word_expression, pos = 4, col = col, drect = T)
        text3(x = curr_x, y = curr_y, cex = cex,
              labels = word_to_write, pos = 4, col = col, drect = T)
        abline(h=curr_y - strheight(correct_l2xp_vec("GIs"), cex = cex) / 2, lwd = 0.5)
      }
      if(justified & (linei != nlines)){
        curr_x <- curr_x + strwidth(word_expression) * cex + justified_space_between_words[linei]
      } else {
        curr_x <- curr_x + strwidth(word_expression) * cex + optimal_word_placement_inf$space_width_min  
      }
      
    }
    curr_x <- rect_coords$x0 - abs(rect_coords$x0 - rect_coords$x1) * str_width_lefter_start_ratio
    curr_y <- curr_y - optimal_word_placement_inf$vertical_space * (1 + str_height_lower_start_ratio) - extra_leading
  }
  
}

swap_one <- function(x, y_i){
  list(append(x[[1]], values = y_i[[1]], after = y_i[[2]] + x[[2]])[-(y_i[[2]] + x[[2]])], x[[2]] + length(y_i[[1]]) - 1)
}

swap <- function(x = vector(), y = list(), inds = vector()){
  #swaps y items into vector x for elements of x at locations given by inds
  x <- list(x, 0)
  y <- y[order(inds)]
  inds <- sort(inds)
  out <- Reduce(f = swap_one, x = c(list(x), lapply(seq_along(y), function(i) list(y[[i]], inds[i]))))
  out[[1]]
}

string_to_tokens <- function(txt_string){
  
  #get basic string
  txt <- strsplit(txt_string, split = " ")[[1]]
  
  #first adjust for new characters
  has_newlines <- grep(pattern = "\n", x = txt)
  newline_swap <- lapply(txt[has_newlines], function(x){
    y <- strsplit(x, "\n")[[1]]
    num_nl <- lengths(regmatches(x, gregexpr("\n", x)))
    y[1:num_nl] <- paste0(y[1:num_nl], "\n")
    y
  })
  txt <- swap(txt, newline_swap, has_newlines)
  
  #find relevant LaTeX notation & modify tokens accordingly
  math_pairs <- do.call(rbind, lapply(seq_along(txt), function(i){
    mathi = which(strsplit(txt[i], "")[[1]] == "$")
    if(length(mathi) != 0){
      return(data.frame(token = i, math_i = mathi))
    } else {
      return(integer(0))
    }
  }))
  
  if(length(math_pairs) != 0){
    math_pairs <- data.frame(cbind(matrix(math_pairs$token, ncol = 2, byrow = T), 
                                   matrix(math_pairs$math_i, ncol = 2, byrow = T)))
    colnames(math_pairs) <- c("open_i", "close_i", "char_i_open", "char_i_close")
    for(i in 1:nrow(math_pairs)){
      if(math_pairs$open_i[i] == math_pairs$close_i[i]){
        next()
      } else {
        txt[math_pairs$open_i[i]] <- paste0(txt[math_pairs$open_i[i]:math_pairs$close_i[i]], collapse = " ")
        txt <- txt[-((math_pairs$open_i[i]+1):math_pairs$close_i[i])]
      }
    }
  }
  
  
  bracket_pairs <- list(
    do.call(rbind, lapply(seq_along(txt), function(i){
      openi <- which(strsplit(txt[i], "")[[1]] == "{")
      if(length(openi) != 0){
        slashi <- which(sapply(strsplit(txt[i], "")[[1]], grepl, pattern = '\t'))
        if(length(slashi) == 0){
          return(data.frame(token_open = i, 
                            open_i = openi,
                            slash_i = "",
                            style = ""))
        } else {
          styles <- sapply(seq_along(slashi), function(ci) substr(txt[i], slashi[ci], openi[ci]))
          return(data.frame(token_open = i, 
                            open_i = openi,
                            slash_i = slashi,
                            style = styles))  
        }
      } else{
        return(integer(0))
      }
    })
    ),
    
    do.call(rbind, lapply(seq_along(txt), function(i){
      closei <- which(strsplit(txt[i], "")[[1]] == "}")
      if(length(closei) != 0){
        return(data.frame(token_close = i, 
                          close_i = closei))
      } else{
        return(integer(0))
      }
    })
    )
  )
  
  bracket_pos <- rbind(data.frame(type = "o", 
                                  pos = bracket_pairs[[1]]$token_open + 
                                    bracket_pairs[[1]]$open_i / 
                                    nchar(txt[bracket_pairs[[1]]$token_open]),
                                  ind = 1:nrow(bracket_pairs[[1]])),
                       data.frame(type = "c", 
                                  pos = bracket_pairs[[2]]$token_close + 
                                    bracket_pairs[[2]]$close_i / 
                                    nchar(txt[bracket_pairs[[2]]$token_close]),
                                  ind = 1:nrow(bracket_pairs[[2]]))
  )
  bracket_pos <- bracket_pos[order(bracket_pos$pos), c("type", "ind")]
  bracket_pos$score <- as.numeric(bracket_pos$type == "o") - as.numeric(bracket_pos$type == "c")
  bracket_pos$cumscore <- cumsum(bracket_pos$score)
  
  bracket_match <- data.frame(do.call(rbind, lapply(1:max(bracket_pos$cumscore), function(mcs){
    cbind(o = bracket_pos$ind[which(bracket_pos$type == "o" & bracket_pos$cumscore == mcs)], 
          c = bracket_pos$ind[which(bracket_pos$type == "c" & bracket_pos$cumscore == (mcs - 1))]) 
  })))
  bracket_match <- bracket_match[order(bracket_match$o),]
  
  bracket_pairs <- do.call(rbind, apply(bracket_match, 1, function(i) cbind(bracket_pairs[[1]][i[1],], bracket_pairs[[2]][i[2],])))
  
  #adjust tokens to reflect text modifications
  for(i in 1:nrow(bracket_pairs)){
    if(bracket_pairs$token_open[i] == bracket_pairs$token_close[i]){
      next()
    } else {
      inds <- bracket_pairs$token_open[i]:bracket_pairs$token_close[i]
      txt[inds] <- sapply(inds, function(j){
        new_txt <- txt[j]
        if(j != bracket_pairs$token_open[i]){
          new_txt <- paste0(bracket_pairs$style[i], new_txt)
        }
        if(j != bracket_pairs$token_close[i]){
          new_txt <- paste0(new_txt, "}")
        }
        new_txt
      })
    }
  }
  
  txt <- gsub(pattern = "\t", replacement = "\\t", x = txt, fixed = T)
  
  list(tokens = gsub(x = txt, pattern = "\n", replacement = ""),
       newlines = grepl(txt, pattern = "\n"))
  
}

correct_l2xp_vec <- function(x){
  # latex2exp::TeX(x) <- this doesn't handle bolditalic correctly
  empties <- which(x == "")
  x[empties] <- "placeholder"
  out <- lapply(x, latex2exp::TeX)
  out_str <- sapply(seq_along(out), function(i) as.character(out[[i]]))
  bis <- intersect(grep(pattern = "bold", out_str), grep(pattern = "italic", out_str))
  new_out_str <- sapply(out_str[bis], function(i) paste0("bolditalic(\"", strsplit(i, "\"")[[1]][2], "\")"))
  new_out_str <- swap(out_str, new_out_str, bis)
  new_out <- sapply(new_out_str, function(i){
    parse(text = i)
  })
  new_out[empties] <- ""
  new_out
}

find_compatible_start <- function(mat, p, q, target_rows = NULL, target_cols = NULL, allow_left_extension = F,
                                  left_extension_prop = 0.25) {
  
  n <- nrow(mat)
  m <- ncol(mat)
  
  # Append a matrix of dimension p x m filled with TRUE to the bottom
  extended_mat <- rbind(mat, matrix(TRUE, p, m))
  
  # Calculate row sums and identify ineligible rows
  row_sums <- apply(extended_mat, 1, sum)
  ineligible_rows <- which(row_sums < q)
  
  # Handling the case where no rows are ineligible
  max_start_row <- n + 1
  if(length(ineligible_rows) == 0) {
    eligible_start_rows <- 1:max_start_row
  } else {
    eligible_start_rows <- setdiff(1:max_start_row, ineligible_rows)
  }
  
  # Apply target row and column constraints if provided
  if (!is.null(target_rows)) {
    eligible_start_rows <- intersect(eligible_start_rows, target_rows)
  }
  eligible_cols <- if (is.null(target_cols)) 1:m else intersect(1:m, target_cols)
  
  # Search for the top-left most submatrix of size p x q
  for (start_row in eligible_start_rows) {
    for (j in eligible_cols) {
      if (start_row + p - 1 <= n + p && j + q - 1 <= m) {
        submatrix <- extended_mat[start_row:(start_row + p - 1), j:(j + q - 1)]
        if (all(submatrix)) {
          
          shortfall <- max(0, start_row + p - n - 1)
          
          if(allow_left_extension){
            leftmost_extension <- min(j, max(1, j - ceiling(left_extension_prop * m)))
            left_extended_submat <- extended_mat[start_row:(start_row + p - 1), 1:m, drop = FALSE]
            all_true_cols <- apply(left_extended_submat, 2, all)
            rle_all_true_cols <- rle(all_true_cols)
            run_inds_all_true_cols <- cbind(c(1, cumsum(rle_all_true_cols$lengths)[-length(rle_all_true_cols$lengths)] + 1), 
                              cumsum(rle_all_true_cols$lengths))
            run_inds_all_true_cols <- run_inds_all_true_cols[rle_all_true_cols$values, , drop = FALSE]
            max_bounds_all_true_cols <- run_inds_all_true_cols[run_inds_all_true_cols[,1] <= j & 
                                                   run_inds_all_true_cols[,2] >= j,]
            j <- max(leftmost_extension, max_bounds_all_true_cols[1])
          }
          
          # Find the maximum width using apply and all along columns
          potential_max_width_submat <- extended_mat[start_row:(start_row + p - 1), j:m, drop = FALSE]
          all_true_cols <- apply(potential_max_width_submat, 2, all)
          
          #max block has to overlap with submatrix -- already starts with it, can maybe also extend left
          max_width <- which(!all_true_cols)[1] - 1
          
          if (is.na(max_width)) { # In case all columns are TRUE
            max_width <- ncol(potential_max_width_submat)
          }
          
          return(list(row = start_row, col = j, shortfall = shortfall, max_width = max_width))
        }
      }
    }
  }
  
  # Return NULL if no matching submatrix is found
  return(NULL)
}


neighbors <- function(grid, pts, val = 0, distance = 1){
  dims <- dim(grid)
  offset_grid <- expand.grid(y = -distance:distance, 
                             x = -distance:distance)
  offset_grid <- offset_grid[!apply(offset_grid, 1, function(x) all(x==0)),]
  neighboring_pts <- pts[rep(1:nrow(pts), each = nrow(offset_grid)), ] + offset_grid[rep(1:nrow(offset_grid), times = nrow(pts)), ]
  neighboring_pts <- neighboring_pts[!duplicated(neighboring_pts),]
  neighboring_pts <- neighboring_pts[neighboring_pts[,1] >= 1 & neighboring_pts[,1] <= dims[1] &
                                       neighboring_pts[,2] >= 1 & neighboring_pts[,2] <= dims[2],]
  return(as.matrix(neighboring_pts[grid[as.matrix(neighboring_pts)] == val,]))
}


has_neighbors <- function(grid, pts, val = 0, distance = 1){
  
  #make useful variables
  dims <- dim(grid)
  offset_grid <- expand.grid(y = -distance:distance, 
                             x = -distance:distance)
  offset_grid <- offset_grid[!apply(offset_grid, 1, function(x) all(x==0)),]
  
  #find neighboring points
  pt_inds <- rep(1:nrow(pts), each = nrow(offset_grid))
  neighboring_pts <- pts[pt_inds, ] + offset_grid[rep(1:nrow(offset_grid), times = nrow(pts)),]
  
  #subset to only valid (inside of grid) points
  in_grid <- neighboring_pts[,1] >= 1 & neighboring_pts[,1] <= dims[1] &
    neighboring_pts[,2] >= 1 & neighboring_pts[,2] <= dims[2]
  pt_inds <- pt_inds[in_grid]
  neighboring_pts <- neighboring_pts[in_grid,]
  
  #find desired points
  target_pt_inds <- unique(pt_inds[grid[as.matrix(neighboring_pts)] == val])
  
  #return pt values
  return(as.matrix(pts[target_pt_inds,]))
}

expand_blocks <- function(grid, pt_sets = NULL, prop_expansion_rate = F, max_expansion_rate = 3) {
  
  if(is.null(pt_sets)){
    grid_vals <- sort(unique(c(0, grid)))
    pt_sets <- lapply(grid_vals, function(i){
      has_neighbors(grid, which(grid == i, T), 0)
    })
    names(pt_sets) <- grid_vals
  }
  
  if(prop_expansion_rate){
    n_border_pts <- sapply(pt_sets[as.character(grid_vals[grid_vals != 0])], nrow)
    expansion_distances <- ceiling(n_border_pts / max(n_border_pts) * max_expansion_rate)
  } else {
    expansion_distances <- setNames(rep(1, length(grid_vals) - 1), grid_vals[grid_vals != 0])
  }
  
  exp_grid <- grid
  while(nrow(pt_sets[["0"]]) != 0){
    expansion_order <- sample(setdiff(grid_vals, 0))
    for(i in expansion_order){
      
      if(nrow(pt_sets[[as.character(i)]]) == 0){
        next()
      }
      
      #find zero-valued neighbors
      neighboring_pts <- neighbors(grid = exp_grid,
                                   pts = pt_sets[[as.character(i)]],
                                   val = 0,
                                   distance = expansion_distances[as.character(i)])
      
      if(nrow(neighboring_pts) == 0){
        pt_sets[[as.character(i)]] <- neighboring_pts
        next()
      } else {
        #swap in to pt_sets and mark in grid
        exp_grid[neighboring_pts] <- i
        pt_sets[[as.character(i)]] <- has_neighbors(exp_grid, neighboring_pts, 0)
        
        #remove from 0-valued pt_set
        zero_vals_comb <- rbind(neighboring_pts, pt_sets[["0"]])
        pt_sets[["0"]] <- pt_sets[["0"]][!duplicated(zero_vals_comb)[-(1:nrow(neighboring_pts))],]
      }
      
    }
  }
  
  return(exp_grid)
}

place_word <- function(i, word_used, word_w, ws_w, v_adj, v_adj_s, word_h, word_h_buffered,
                       grid, word_locations, fast_invalid, verbose, n_words_per_sentence){
  
  if(all(word_used)){
    return(list(word_used = word_used, 
                grid = grid, 
                word_locations = word_locations))
  }
  
  # if(i == 29 & !is.null(word_locations[[i]]) && nrow(word_locations[[i]]) == 2){
  #   keep_going <- F
  #   break
  # }
  
  #process current state of grid
  start_position_info <- find_compatible_start(mat = grid == 0 | grid == i, 
                                               p = ceiling(word_h_buffered),
                                               q = ceiling(word_w[!word_used][1] + ws_w), 
                                               target_rows = ifelse2(is.null(word_locations[[i]]),
                                                                     NULL,
                                                                     tail(word_locations[[i]], 1)$ri + ceiling(word_h)),
                                               # target_rows = ifelse2(is.null(word_locations[[i]]), 
                                               #                       NULL, 
                                               #                       max(which(grid == i, arr.ind=T)[,1])),
                                               target_cols = ifelse2(is.null(word_locations[[i]]),
                                                                     NULL,
                                                                     (tail(word_locations[[i]], 1)$ci - ceiling(word_w[!word_used][1] / 2)):
                                                                       (tail(word_locations[[i]], 1)$ci + 
                                                                          ceiling(tail(word_locations[[i]], 1)$width) - ceiling(word_w[!word_used][1] / 2)
                                                                       )
                                               )
                                               # target_cols = NULL
  )
  
  # grid[max(1, (start_position_info$row)) :
  #        min(nrow(grid), (start_position_info$row + ceiling(word_h) - 1)),
  #      max(1, (start_position_info$col)) :
  #        min(ncol(grid), (start_position_info$col + start_position_info$max_width - 1))]
  
  # grid[max(1, (start_position_info$row-1)) :
  #        min(nrow(grid), (start_position_info$row + ceiling(word_h) - 0)),
  #      max(1, (start_position_info$col-1)) :
  #        min(ncol(grid), (start_position_info$col + start_position_info$max_width - 0))]
  
  
  #try another starting location if this failed to find anything
  if(is.null(start_position_info)){
    
    if(verbose){cat(paste0(" (bs) "))}
    
    if(fast_invalid){
      #mark all currently filled sentence locations as invalid
      grid[grid == i] <- -1
      
    } else {
      #empty filled spots for sentence
      grid[grid == i] <- 0
      
      #mark only first word location as invalid
      string_info <- head(word_locations[[i]], 1)
      end_col <- min(ceiling(string_info$ci + string_info$width - 1), ncol(grid))
      end_row <- ceiling(string_info$ri + string_info$height - 1)
      cols_filled <- string_info$ci : end_col
      rows_filled <- string_info$ri : end_row
      
      grid[rows_filled, cols_filled] <- -1
    }
    
    
    #clear other tracking variables
    word_used <- rep(F, n_words_per_sentence[i])
    word_locations[i] <- list(NULL)
    
    return(list(word_used = word_used, 
                grid = grid, 
                word_locations = word_locations))
    
  }
  
  if(verbose){cat(paste0(" (row: ", start_position_info$row, ")"))}
  
  #process sentence to be added
  words_needing_to_be_placed <- which(!word_used)
  cumulative_word_width <- cumsum(word_w[words_needing_to_be_placed]) + 
    cumsum(c(rep(ws_w, length(words_needing_to_be_placed)))) #also include ws at the end of the line
  words_that_fit <- words_needing_to_be_placed[cumulative_word_width < start_position_info$max_width]
  string_to_place <- paste0(names(dims[[i]]$w[words_that_fit]), collapse = " ")
  if(verbose){cat(paste0(" (words: ", string_to_place, ")"))}
  
  
  #compile new sentence info
  string_info <- data.frame(string = string_to_place,
                            ci = start_position_info$col,
                            ri = start_position_info$row,
                            width = as.numeric(cumulative_word_width[length(words_that_fit)]),
                            v_adj = v_adj,
                            height = word_h,
                            buffered_height = word_h_buffered,
                            v_adj_s = v_adj_s
  )
  
  
  #left or right justify text if possible, and if not, get as close as possible
  #otherwise prefer left justification in left half and right justification in right half of image
  justification <- NA
  if(!is.null(word_locations[[i]])){
    
    # if(i == 4 & nrow(word_locations[[i]]) == 2){
    #   keep_going <- F
    #   break
    # }
    
    previous_string_info <- tail(word_locations[[i]], 1)
    prev_on_the_right <- (previous_string_info$ci + previous_string_info$width / 2) > (ncol(grid) / 2)
    
    row_bounds <- c(tr = string_info$ri,
                    br = min(ceiling(string_info$ri + word_h_buffered - 1), nrow(grid)))
    
    lj_bounds <- c(lc = previous_string_info$ci,
                   rc = previous_string_info$ci + ceiling(string_info$width) - 1)
    
    rj_bounds <- c(rc = previous_string_info$ci + ceiling(previous_string_info$width) - 1)
    rj_bounds <- c(lc = as.numeric(rj_bounds["rc"]) - ceiling(string_info$width) + 1, rj_bounds)
    
    c_bounds <- c(cc = previous_string_info$ci + ceiling(previous_string_info$width / 2) - 1)
    c_bounds <- c(lc = as.numeric(c_bounds["cc"]) - ceiling(string_info$width / 2) + 1)
    c_bounds <- c(c_bounds["lc"],
                  rc = as.numeric(c_bounds["lc"]) + ceiling(string_info$width) - 1)
    
    
    if(string_info$ri <= nrow(grid)){
      
      #find maximum extent of horizontal stretch that contains start pos
      subgrid <- grid[row_bounds["tr"] : row_bounds["br"], , drop = F]
      rle_cols <- rle(apply(subgrid == 0 | subgrid == i, 2, all))
      run_inds <- cbind(c(1, cumsum(rle_cols$lengths)[-length(rle_cols$lengths)] + 1), 
                        cumsum(rle_cols$lengths))
      run_inds <- run_inds[rle_cols$values, , drop = FALSE]
      max_bounds <- run_inds[run_inds[,1] <= string_info$ci & run_inds[,2] >= string_info$ci,]
      names(max_bounds) <- c("lc", "rc")
      
      #check left justified
      if(all(lj_bounds >= 1 & lj_bounds <= ncol(grid))){
        subgrid <- grid[row_bounds["tr"] : row_bounds["br"],
                        lj_bounds["lc"] : lj_bounds["rc"]]
        left_possible <- all(subgrid == 0 | subgrid == i)
      } else {
        left_possible <- F
      }
      overlap_bounds <- c(lc = max(max_bounds["lc"], lj_bounds["lc"]),
                          rc = min(max_bounds["rc"], lj_bounds["rc"]))
      lj_closest_bounds <- c(lc = as.numeric(overlap_bounds["lc"]),
                             rc = as.numeric(overlap_bounds["lc"]) + ceiling(string_info$width) - 1)
      if(lj_closest_bounds["rc"] > max_bounds["rc"]){
        lj_closest_bounds <- lj_closest_bounds - (lj_closest_bounds["rc"] - max_bounds["rc"])
      }
      lj_closest_possible <- list(p = diff(overlap_bounds) / ceiling(string_info$width),
                                  closest_bounds = lj_closest_bounds)
      
      #check right justified
      if(all(rj_bounds >= 1 & rj_bounds <= ncol(grid))){
        subgrid <- grid[row_bounds["tr"] : row_bounds["br"],
                        rj_bounds["lc"] : rj_bounds["rc"]]
        right_possible <- all(subgrid == 0 | subgrid == i)
      } else {
        right_possible <- F
      }
      overlap_bounds <- c(lc = max(max_bounds["lc"], rj_bounds["lc"]),
                          rc = min(max_bounds["rc"], rj_bounds["rc"]))
      rj_closest_bounds <- c(lc = as.numeric(overlap_bounds["lc"]),
                             rc = as.numeric(overlap_bounds["lc"]) + ceiling(string_info$width) - 1)
      if(rj_closest_bounds["rc"] > max_bounds["rc"]){
        rj_closest_bounds <- rj_closest_bounds - (rj_closest_bounds["rc"] - max_bounds["rc"])
      }
      rj_closest_possible <- list(p = diff(overlap_bounds) / ceiling(string_info$width),
                                  closest_bounds = rj_closest_bounds)
      
      #check centered
      if(all(c_bounds >= 1 & c_bounds <= ncol(grid))){
        subgrid <- grid[row_bounds["tr"] : row_bounds["br"],
                        c_bounds["lc"] : c_bounds["rc"]] 
        center_possible <- all(subgrid == 0 | subgrid == i)
      } else {
        center_possible <- F
      }
      overlap_bounds <- c(lc = max(max_bounds["lc"], c_bounds["lc"]),
                          rc = min(max_bounds["rc"], c_bounds["rc"]))
      c_closest_bounds <- c(lc = as.numeric(overlap_bounds["lc"]),
                            rc = as.numeric(overlap_bounds["lc"]) + ceiling(string_info$width) - 1)
      if(c_closest_bounds["rc"] > max_bounds["rc"]){
        c_closest_bounds <- c_closest_bounds - (c_closest_bounds["rc"] - max_bounds["rc"])
      }
      c_closest_possible <- list(p = diff(overlap_bounds) / ceiling(string_info$width),
                                 closest_bounds = c_closest_bounds)
      
    } else { #if we are off the bottom of the grid
      
      max_bounds <- c(lc = 1, rc = ncol(grid))
      
      #check left justified
      if(all(lj_bounds >= 1 & lj_bounds <= ncol(grid))){
        left_possible <- T
      } else {
        left_possible <- F
        overlap_bounds <- c(lc = max(max_bounds["lc"], lj_bounds["lc"]),
                            rc = min(max_bounds["rc"], lj_bounds["rc"]))
        lj_closest_bounds <- c(lc = as.numeric(overlap_bounds["lc"]),
                               rc = as.numeric(overlap_bounds["lc"]) + ceiling(string_info$width) - 1)
        if(lj_closest_bounds["rc"] > max_bounds["rc"]){
          lj_closest_bounds <- lj_closest_bounds - (lj_closest_bounds["rc"] - max_bounds["rc"])
        }
        lj_closest_possible <- list(p = diff(overlap_bounds) / ceiling(string_info$width),
                                    closest_bounds = lj_closest_bounds)
      }
      
      #check right justified
      if(all(rj_bounds >= 1 & rj_bounds <= ncol(grid))){
        right_possible <- T
      } else {
        right_possible <- F
        overlap_bounds <- c(lc = max(max_bounds["lc"], rj_bounds["lc"]),
                            rc = min(max_bounds["rc"], rj_bounds["rc"]))
        rj_closest_bounds <- c(lc = as.numeric(overlap_bounds["lc"]),
                               rc = as.numeric(overlap_bounds["lc"]) + ceiling(string_info$width) - 1)
        if(rj_closest_bounds["rc"] > max_bounds["rc"]){
          rj_closest_bounds <- rj_closest_bounds - (rj_closest_bounds["rc"] - max_bounds["rc"])
        }
        rj_closest_possible <- list(p = diff(overlap_bounds) / ceiling(string_info$width),
                                    closest_bounds = rj_closest_bounds)
      }
    }
    
    #now shift the indices over appropriately
    
    #both possible
    if(left_possible & right_possible){
      
      if(prev_on_the_right){
        #right justify
        string_info$ci <- rj_bounds["lc"]
        justification <- "right"
      } else{
        #left justify
        string_info$ci <- lj_bounds["lc"]
        justification <- "left"
      }
      
      #only one possible
    } else if(left_possible){
      #left justify
      string_info$ci <- lj_bounds["lc"]
      justification <- "left"
      
    } else if(right_possible){
      #right justify
      string_info$ci <- rj_bounds["lc"]
      justification <- "right"
      
      #neither possible, but one much gets closer
      #can also consider chopping off words
    } else if(abs(rj_closest_possible[["p"]] - lj_closest_possible[["p"]]) > 0.25){
      if(rj_closest_possible[["p"]] > lj_closest_possible[["p"]]){
        #right justify
        string_info$ci <-rj_closest_possible[["closest_bounds"]]["lc"]
        justification <- "more_right"
      } else {
        #left justify
        string_info$ci <- lj_closest_possible[["closest_bounds"]]["lc"]
        justification <- "more_left"
      }
      
      #neither possible, but both are approx. equally close
    } else {
      if(prev_on_the_right){
        #right justify
        string_info$ci <- rj_closest_possible[["closest_bounds"]]["lc"]
        justification <- "more_right"
      } else{
        #left justify
        string_info$ci <- lj_closest_possible[["closest_bounds"]]["lc"]
        justification <- "more_left"
      }
    }
    
  }
  
  string_info$just <- justification
  
  #mark words as having been placed
  word_used[words_that_fit] <- T
  
  #process grid locations
  end_col <- ceiling(string_info$ci + string_info$width - 1)
  end_col <- min(end_col, ncol(grid))
  # end_row <- ceiling(string_info$ri + ifelse(all(word_used), word_h_buffered, word_h) - 1)
  end_row <- ceiling(string_info$ri + word_h_buffered - 1)
  string_info$end_ci <- end_col
  string_info$end_ri <- end_row
  
  cols_filled <- floor(string_info$ci) : end_col
  rows_filled <- string_info$ri : end_row
  
  #expand grid if necessary
  if(end_row > nrow(grid)){
    grid <- rbind(grid, matrix(0, ncol = ncol(grid), nrow = end_row - nrow(grid)))
    if(verbose){cat("\n\nexpand")}
  }
  
  #check all proposed spots are unoccupied
  if(!all(grid[rows_filled, cols_filled] == 0 |
          grid[rows_filled, cols_filled] == i)){
    keep_going <- F
    cat("\n\nsomething is wrong")
    break
  }
  
  #add string to image and record locations
  word_locations[[i]] <- rbind(word_locations[[i]], string_info)
  grid[rows_filled, cols_filled] <- i
  
  #return updated data objects
  return(list(word_used = word_used, 
              grid = grid, 
              word_locations = word_locations))
  
}

place_sentence <- function(dims, grid, v_adj_factor, v_adj_factor_sentences,
                           mean_word_h, fast_invalid = T, used, word_locations,
                           verbose = F, n_words_per_sentence, cex_scale){
  
  if(all(used)){
    return(list(grid = grid, 
                used = used, 
                word_locations = word_locations, 
                keep_going = keep_going))
  }
  
  #retrieve next word index
  i <- min(which(!used))
  word_used <- rep(F, n_words_per_sentence[i])
  
  #retrieve string metadata
  word_w <- dims[[i]]$w * cex_scale
  ws_w <- dims[[i]]$ws * cex_scale
  v_adj <- dims[[i]]$h * cex_scale * v_adj_factor
  v_adj_s <- dims[[i]]$h * cex_scale * v_adj_factor_sentences
  word_h <- dims[[i]]$h * cex_scale + v_adj
  word_h_buffered <- word_h + v_adj_s
  
  if(verbose){cat(paste0("\n\n", i))}else{cat(paste0("(", i, ") "))}
  
  while(!all(word_used)){
    word_out <- place_word(i = i, 
                           word_used = word_used, 
                           word_w = word_w, 
                           ws_w = ws_w, 
                           v_adj = v_adj, 
                           v_adj_s = v_adj_s, 
                           word_h = word_h, 
                           word_h_buffered = word_h_buffered,
                           grid = grid, 
                           word_locations = word_locations, 
                           fast_invalid = fast_invalid, 
                           verbose = verbose, 
                           n_words_per_sentence = n_words_per_sentence)
    
    word_used <- word_out$word_used
    grid <- word_out$grid
    word_locations <- word_out$word_locations
    
  }
  
  #reset temporary forbidden spots in grid
  grid[grid == -1] <- 0
  
  #mark sentence as used
  used[i] <- T
  if(all(used)){
    cat("\n\nWOO-HOO finished placing sentences")
  }
  
  return(list(grid = grid, 
              used = used, 
              word_locations = word_locations, 
              keep_going = keep_going))
  
  
}

#### simulation parameters ####
n_sentences <- 40
n_words_per_sentence <- sample(5:12, n_sentences, T)
sentences <- sapply(n_words_per_sentence, function(n_words){
  out <- paste0(sample(LETTERS, 1), paste0(replicate(n_words, paste0(sample(letters, sample(2:10, 1)), collapse = "")),
          collapse = " "),  sample(c(".", "?", "!"), 1), collapse = "")
  
})

#### real parameters ####
#read text in from disk
sentences <- readLines("~/Documents/Documents - nikolai/cam_commendations.txt", warn = F)
sentences <- trimws(sentences[1:(which(sentences == "")-1)])
sentences <- sample(sentences)
n_sentences <- length(sentences)
n_words_per_sentence <- as.numeric(sapply(sentences, function(x) {
  return(length(strsplit(x, " ")[[1]]))
}))


#### start loop ####
library(foreach)
library(doParallel)
library(parallel)

#initialize parallelization
foreach_parallel <- T
if(foreach_parallel){
  if(!exists("cl")){
    cl <- makeCluster(12, outfile="")
    registerDoParallel(cl)
  }
  getDoParWorkers()
}

# for(run_index in 1:100){
# out <- foreach(run_index=1:500, .packages = c("extrafont"), .combine = 'c') %dopar% {
run_index <- 1

#### metadata ####
rotation_bounds <- c(-10,10)
rotation_angles <- rbeta(n_sentences, 5, 5) * diff(rotation_bounds) + rotation_bounds[1]
img_width <- 5
img_height <- 5

#find smarter estimate for grid resolution:
#fixed value
grid_res <- 50
#or approx one cell per character? assuming cex_scale_scale = 1
cex_scale_scale <- 0.7
grid_res <- ceiling(sqrt(sum(sapply(sentences, nchar)) / sqrt(img_width * img_height)))

width <- img_width * grid_res
height <- img_height * grid_res
# color_options <- colors()
# color_options <- c(
#   "white", "lightgray", "lavender", "lightpink", "lightgreen",
#   "limegreen", "orange", "hotpink", "gold",
#   "yellow", "cyan", "magenta", "coral", "turquoise"
# )
font_options <- sample(fonts(), n_sentences, T)
font_options <- c(
  "Helvetica", 
  "Times New Roman", 
  "Arial Unicode MS", 
  "Courier", 
  "Georgia", 
  "Verdana", 
  "Palatino", 
  "Garamond", 
  "Lucida Grande",
  "Trebuchet MS", 
  "Futura", 
  "Arial", 
  "Baskerville",
  "Gill Sans", "Optima", "Hoefler Text", "Caslon", 
  "Charter"
)
font_families <- sample(font_options, n_sentences, T)

font_styles <- sample(1:4, n_sentences, T)
font_sizes <- runif(n_sentences, 1, 3)
font_sizes <- font_sizes / mean(font_sizes)

empty_plot <- function(width, height){
  par(xpd = NA, mar = c(0,0,0,0))
  plot(NULL, xlim = c(0, width), ylim = c(0, height))  
}

#### dims ####

dims <- lapply(1:length(sentences), function(i){
  
  # Open a temporary plotting device to measure text
  cairo_pdf(file = tempfile(), width = img_width, height = img_height)
  empty_plot(width, height)
  par(family = font_families[i])
  
  # Measure dimensions
  indiv_words <- strsplit(sentences[i], " ")[[1]]
  w <- setNames(
    strwidth(indiv_words, units = "user", font = font_styles[i], cex = font_sizes[i]), 
    indiv_words)
  tw <- strwidth(sentences[i], units = "user", font = font_styles[i], cex = font_sizes[i])
  h <- strheight(sentences[i], units = "user", 
                 font = font_styles[i], cex = font_sizes[i])
  ws <- strwidth(" ", units = "user", 
                 font = font_styles[i], cex = font_sizes[i])
  
  #close graphical device
  dev.off()
  
  return(list(h = h, w = w, ws = ws, total = tw * h))
})

# dims <- lapply(1:length(sentences), function(i){
# 
#   # Open a temporary plotting device to measure text
#   cairo_pdf(file = tempfile(), width = img_width, height = img_height)
#   empty_plot(width, height)
#   
#   # Measure dimensions
#   bb <- text3(x = img_width / 2, y = img_height / 2, labels = sentences[i],
#               cex = font_sizes[i], col = 1, family = font_families[i],
#               font = font_styles[i], pos = c(1,4), drect = T, return_bb = T)
#   tw <- bb$xright - bb$xleft
#   h <- (bb$ytop - bb$ybottom) * 1.05
#   
#   w <- sapply(strsplit(sentences[i], " ")[[1]], function(wi){
#     bb <- text3(x = img_width / 2, y = img_height / 2, labels = wi,
#                 cex = font_sizes[i], col = 1, family = font_families[i],
#                 font = font_styles[i], pos = c(1,4), drect = T, return_bb = T)
#     return(bb$xright - bb$xleft)
#   })
#   
#   bb_ws <- text3(x = img_width / 2, y = img_height / 2, labels = " ",
#                  cex = font_sizes[i], col = 1, family = font_families[i],
#                  font = font_styles[i], pos = c(1,4), drect = T, return_bb = T)
#   ws <- bb_ws$xright - bb_ws$xleft
#   
#   #close graphical device
#   dev.off()
# 
#   return(list(h = h, w = w, ws = ws, total = tw * h))
# })

#process a few graphical parameters
total_img_area <- height * width
total_sentence_area <- sum(sapply(dims, function(x) x$total))
cex_scale <- sqrt(total_img_area / total_sentence_area) * cex_scale_scale

ifelse2 <- function(test, yes, no) if(test) yes else no

#### process text locations ####

#initialize text parameters
v_adj_factor <- 0.2 #this is the spacing between lines as a proportion of char height
v_adj_factor_sentences <- 0.4 #this is the spacing between sentences as a proportion of char height
mean_word_h <- mean(sapply(dims, function(x) x$h)) #have a variable and a fixed component here?
fast_invalid <- T
used <- rep(F, n_sentences)
grid <- matrix(data = 0, nrow = height, ncol = width)
word_locations <- lapply(1:n_sentences, function(i) NULL)
verbose <- T
keep_going <- T

#run loop to place words
while(!all(used)){
  
  out <- place_sentence(dims = dims, grid = grid, v_adj_factor = v_adj_factor, v_adj_factor_sentences = v_adj_factor_sentences, 
                        mean_word_h = mean_word_h, used = used, fast_invalid = T, word_locations = word_locations,
                 verbose = F, n_words_per_sentence = n_words_per_sentence, cex_scale = cex_scale)
  grid <- out$grid
  used <- out$used
  word_locations <- out$word_locations
  keep_going <- out$keep_going
  
  if(!keep_going) break
  
  
}

# grid <- expand_blocks(grid)
concave_hulls <- lapply((1:n_sentences)[!sapply(word_locations, is.null)], function(i){
  sentence_inds <- which(grid == i, arr.ind = T)
  
  # lapply(1:nrow(word_locations[[i]]), function(li){
  #   string_info <- word_locations[[i]][li,]
  #   cols_filled <- floor(string_info$ci) : string_info$end_ci
  #   rows_filled <- string_info$ri : string_info$end_ri
  # })
  
  sentence_inds[,1] <- nrow(grid) - sentence_inds[,1] + 1
  concave_hull <- concaveman::concaveman(sentence_inds, concavity = -1)
  kde <- MASS::kde2d(sentence_inds[,1], y = sentence_inds[,2], n = 50)
  label_ind <- which(kde$z == max(kde$z), arr.ind = TRUE)
  list(concave_hull = concave_hull, label_loc = c(x = kde$y[label_ind[2]], y = kde$x[label_ind[1]]))
})


#### draw picture ####

#sample colors
background_color <- "#191b30"
color_options <- as.character(wesanderson::wes_palette("Zissou1"))
color_options <- as.character(wesanderson::wes_palette("FantasticFox1"))
color_options <- c(as.character(wesanderson::wes_palette("Royal1")),
                   as.character(wesanderson::wes_palette("Royal2")))
# color_options <- as.character(wesanderson::wes_palette("Darjeeling1"))
colors <- sample(color_options, n_sentences, T)

#sample inverse propto adjacent colors
label_locs <- do.call(rbind, lapply(concave_hulls, function(x) x$label_loc))
chull_dists <- as.matrix(dist(label_locs))
colors <- rep(NA, n_sentences)
softmax <- function(x) exp(x) / sum(exp(x))
for(i in 1:length(concave_hulls)){
  
  if(i == 1){
    colors[i] <- sample(color_options, 1, T)
  } else {
    dists <- setNames(chull_dists[i, 1:(i-1)], colors[1:(i-1)])
    unlisted_cols <- setdiff(color_options, names(dists))
    dists <- c(dists, setNames(rep(max(dists) * 1E3, length(unlisted_cols)), unlisted_cols))
    dists <- sapply(split(dists, names(dists)), min)
    dists <- dists / mean(dists)
    colors[i] <- sample(names(dists), 1, T, prob = softmax(dists))
  }
}


cairo_pdf(file = paste0("~/ccc/ccc_inds_", paste0(rep(0, 3-nchar(run_index)), collapse = ""), run_index, ".pdf"), 
          width = img_width, height = img_height * nrow(grid) / height)
par(xpd = NA)

empty_plot(ncol(grid), nrow(grid))
for(i in 1:length(concave_hulls)){
  polygon(concave_hulls[[i]]$concave_hull[,2:1], col = adjustcolor(colors[i], 0.5))
  text(x = concave_hulls[[i]]$label_loc[["x"]], y = concave_hulls[[i]]$label_loc[["y"]], 
       labels = i, font = 2)
    
}

dev.off()

cairo_pdf(file = paste0("~/ccc/ccc_", paste0(rep(0, 3-nchar(run_index)), collapse = ""), run_index, ".pdf"), 
          width = img_width, height = img_height * nrow(grid) / height)
empty_plot(ncol(grid), nrow(grid))
rect(par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4],
     col = background_color)
for(i in 1:length(word_locations)){
  par(family = font_families[i])
  wl <- word_locations[[i]]
  if(is.null(wl)) next
  for(j in 1:nrow(wl)){
    xloc <- wl$ci[j]
    yloc <- nrow(grid) - wl$ri[j] + 1 - wl$height[j] / 2 - wl$v_adj_s[j] / 2
    
    #POSTPROCESSING
    #align word right if we are super close to border
    if((ncol(grid) - (xloc + ceiling(wl$width[j]))) / ncol(grid) < 0.05){
      subgrid <- grid[wl$ri[j] : (wl$ri[j] + ceiling(wl$height[j])),
                      xloc : ncol(grid)]
      if(all(subgrid == 0 | subgrid == i)){
        xloc <- ncol(grid) - ceiling(wl$width[j]) + 1
      }
    }
    
    #draw text
    text3(x = xloc, 
          y = yloc, 
          labels = wl$string[j],
         cex = font_sizes[i] * cex_scale, col = colors[i],
         family = font_families[i],
         font = font_styles[i], pos = c(4), drect = F, dline = F, vexp = T)
  }

}

dev.off()

# save(list = c("grid", "word_locations"), 
#      file = paste0("~/ccc/ccc_", paste0(rep(0, 3-nchar(run_index)), collapse = ""), run_index, ".RData"), 
#      envir = .GlobalEnv)

sink("~/ccc/ccc_vals.csv", append = T)
cat(paste0(run_index, ", ", round(mean(grid == 0), 3), ", ", round(5 * nrow(grid) / ncol(grid), 3),"\n"))
sink()

list(grid = grid,
     word_locations = word_locations,
     sentences = sentences,
     width = width,
     height = height,
     grid_res = grid_res,
     font_families = font_families,
     font_styles = font_styles,
     font_sizes = font_sizes,
     cex_scale = cex_scale,
     concave_hulls = concave_hulls)

# }

summs <- read.csv("~/ccc/ccc_vals.csv")
summs <- summs[order(summs[,2]),]

head(summs, 20)

# mean(grid == 0)
# paste0("5 x ", 5 * nrow(grid) / ncol(grid))

#TODO
# For angle functionality, don’t think about rotating the text — 
#rotate <<grid>> and pad it into a rectangle with ‘-2’s, 
#and then unrotate / unpad it once the word has been written
# this can also accommodate random shapes — 
#just initialize grid with forbidden locations set to -2

#during initial rounds, block off sections of page at random to simulate random scrawling 
#also try rotation thing maybe +/- 10 degrees
# Change size of font throughout message
# Seed initial more orderly blocks of text in random locations (median font?)
# Give shorter messages a larger font, on average