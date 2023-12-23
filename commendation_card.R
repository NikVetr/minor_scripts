# Install and load necessary packages
library(extrafont)

#### functions ####
remove_top <- function(x, replacement){
  notop <- gsub("\\^s*\\{[^\\)]+\\}", replacement, x)
  notop <- gsub("\\^[a-z|0-9|A-Z]{1}", replacement, notop)
  notop
}

text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex, ...) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex, ...) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(labels, cex = cex, ...) / 2, 
         xright = adj_x + strwidth(labels, cex = cex, ...) / 2, 
         ybottom = adj_y - strheight(labels, cex = cex, ...) / 2, 
         ytop = adj_y + strheight(labels, cex = cex, ...) / 2)
    # abline(h = adj_y - strheight(labels, cex = cex) / 2, lwd = 0.5)
  }
}

remove_tb <- function(x, replacement){
  remove_top(remove_bottom(x, replacement), replacement)
}

remove_bottom <- function(x, replacement){
  nobot <- gsub("g|j|p|q|y|,|\\(|\\)|Q", replacement, x)
  nobot <- gsub("\\_s*\\{[^\\)]+\\}", replacement, nobot) #underscore in brackets
  nobot <- gsub("_[a-z|0-9|A-Z]{1}", replacement, nobot) #underscore w/ one letter following
  nobot
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

correct_l2xp_vec <- function(x){
  # latex2exp::TeX(x) <- this doesn't handle bolditalic correctly
  empties <- which(x == "")
  x[empties] <- "placeholder"
  out <- lapply(x, latex2exp::TeX)
  out_str <- sapply(seq_along(out), function(i) as.character(out[[i]]))
  bis <- intersect(grep(pattern = "bold", out_str), grep(pattern = "italic", out_str))
  new_out_str <- sapply(out_str[bis], function(i) paste0("bolditalic(\"", strsplit(i, "\"")[[1]][2], "\")"))
  new_out_str <- swap(out_str, new_out_str, bis)
  new_out <- sapply(new_out_str, function(i) parse(text = i))
  new_out[empties] <- ""
  new_out
}

text3 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, col = 1, replacement = "a", ...){
  
  #convert text label to expression
  if(all(class(labels) == "character")){
    word_expression <- correct_l2xp_vec(labels)  
  } else {
    word_expression <- labels
  }
  
  #find general params
  strw <- strwidth(word_expression, cex = cex)
  strh <- strheight(word_expression, cex = cex)
  base_strh <- strheight(correct_l2xp_vec("GIs"), cex = cex)
  
  #adjust base location
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strw / 2, 0)
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
  
  #print the text itself
  text(x = adj_x, y = adj_y + adj_ledding, labels = word_expression, pos = NULL, cex = cex, col = col, ...)
  
  #draw a box around it if desired
  if(drect){
    rect(xleft = adj_x - strw / 2, 
         xright = adj_x + strw / 2, 
         ybottom = adj_y - strh / 2 + adj_ledding, 
         ytop = adj_y + strh / 2 + adj_ledding, border = col)
    abline(h=y - strheight(latex2exp::TeX("GIs"), cex = cex) / 2, lwd = 0.5)
  }
}


#### simulation parameters ####
n_sentences <- 3
n_words_per_sentence <- sample(5:12, n_sentences, T)
sentences <- sapply(n_words_per_sentence, function(n_words){
  out <- paste0(sample(LETTERS, 1), paste0(replicate(n_words, paste0(sample(letters, sample(2:10, 1)), collapse = "")),
          collapse = " "),  sample(c(".", "?", "!"), 1), collapse = "")
  
})

#metadata
img_width <- 14
img_height <- 7
grid_res <- 10
width <- img_width * grid_res
height <- img_height * grid_res
colors <- sample(colors(), n_sentences, T)
font_families <- sample(fonts(), n_sentences, T)
font_families <- sample(c("Helvetica", "Times New Roman", "Arial Unicode MS")[3], n_sentences, T)
font_styles <- sample(1:4, n_sentences, T)
font_sizes <- runif(n_sentences, 1, 3)

empty_plot <- function(width, height){
  plot(NULL, xlim = c(0, width), ylim = c(0, height))  
}

dims <- lapply(1:length(sentences), function(i){
  
  # Open a temporary plotting device to measure text
  cairo_pdf(file = tempfile(), width = img_width, height = img_height)
  empty_plot(width, height)
  par(family = font_families[i])
  plot(1, 1, xlim = c(0, width), ylim = c(0, height))
  
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

#process a few graphical parameters
total_img_area <- height * width
total_sentence_area <- sum(sapply(dims, function(x) x$total))
cex_scale <- sqrt(total_img_area / total_sentence_area)

#initiate writing in of text:
used <- rep(F, n_sentences)
grid <- matrix(data = 0, nrow = img_height * grid_res, ncol = img_width * grid_res)
word_locations <- lapply(1:n_sentences, function(i) NULL)
print_progress <- T
while(!all(used)){
    
  # i <- sample(which(!used), 1)
  i <- min(which(!used))
  word_used <- rep(F, n_words_per_sentence[i])
  
  if(print_progress){cat(paste0("\n\n", i))}
  
  while(!all(word_used)){
    
    #process current state of grid
    free_cells <- which(grid == 0, arr.ind = T)
    n_free_in_rows <- table(free_cells[,1])
    n_free_in_rows <- c(n_free_in_rows, setNames(rep(0, nrow(grid) - length(n_free_in_rows)), 
             setdiff(1:nrow(grid), as.numeric(names(n_free_in_rows)))))
    n_free_in_rows <- n_free_in_rows[order(as.numeric(names(n_free_in_rows)))]
    
    #constrain new words in sentence to always fall below old ones
    #TODO enforce vertical adjacency constraint, 
    #and keep trying until you successfully place whole sentence
    if(!is.null(word_locations[[i]])){
      last_row_filled <- max(which(grid == i, arr.ind = T)[,1])
      n_free_in_rows[1:last_row_filled] <- 0
    }
    
    word_widths <- dims[[i]]$w * cex_scale
    whitespace <- dims[[i]]$ws * cex_scale
    
    #process line to start writing
    first_free_row <- min(as.numeric(names(n_free_in_rows)[which(n_free_in_rows > word_widths[!word_used][1])]))
    if(first_free_row == Inf){
      first_free_row <- nrow(grid) + 1
      n_free_in_rows[first_free_row] <- ncol(grid)
    }
    width_available <- n_free_in_rows[first_free_row]
    if(print_progress){cat(paste0(" (putting in row: ", first_free_row, ")"))}
    
    #TODO check that there is space underneath the row, not just width
    
    #process sentence to be added
    words_needing_to_be_placed <- which(!word_used)
    cumulative_word_width <- cumsum(word_widths[words_needing_to_be_placed]) + 
      cumsum(c(0, rep(whitespace, length(words_needing_to_be_placed) - 1)))
    words_that_fit <- words_needing_to_be_placed[cumulative_word_width < width_available]
    string_to_place <- paste0(names(dims[[i]]$w[words_that_fit]), collapse = " ")
    if(print_progress){cat(paste0(" (putting in words: ", string_to_place, ")"))}
    
    #add sentence to image and record locations
    string_info <- data.frame(string = string_to_place,
                              ci = (width - n_free_in_rows[first_free_row]) + 1,
                              ri = first_free_row,
                              width = as.numeric(cumulative_word_width[length(words_that_fit)])
    )
    word_locations[[i]] <- rbind(word_locations[[i]], 
                                 string_info)
    word_used[words_that_fit] <- T
    
    #fill in grid
    cols_filled <- string_info$ci : min(ceiling(string_info$ci + string_info$width), ncol(grid))
    rows_filled <- string_info$ri : ceiling(string_info$ri + dims[[i]]$h * cex_scale)
    
    #expand grid if necessary
    if(max(rows_filled) > nrow(grid)){
      grid <- rbind(grid, matrix(0, ncol = ncol(grid), nrow = max(rows_filled) - nrow(grid) ))
    }
    
    #fill in grid for next round
    grid[rows_filled, cols_filled] <- i
  }
  
  used[i] <- T

}

#### draw picture ####
cairo_pdf(file = "~/ccc_inds.pdf", width = img_width, height = img_height)
par(xpd = NA)

empty_plot(ncol(grid), nrow(grid))

for(i in 1:nrow(grid)){
  for(j in 1:ncol(grid))
    if(grid[i, j] != 0){
      # points(x = j, y = i, pch = 15, col = colors[grid[i, j]])
      text(x = j, y = nrow(grid) - i + 1, pch = 15, labels = grid[i, j])
    }
}

dev.off()

cairo_pdf(file = "~/ccc.pdf", width = img_width, height = img_height)
par(xpd = NA)

empty_plot(ncol(grid), nrow(grid))
for(i in 1:length(word_locations)){
  wl <- word_locations[[i]]
  if(is.null(wl)) break
  for(j in 1:nrow(wl)){
    text2(x = wl$ci[j], y = nrow(grid) - wl$ri[j] + 1, labels = wl$string[j],
         cex = font_sizes[i] * cex_scale, col = colors[i],
         #family = font_families[i],
         font = font_styles[i], pos = c(1,4), drect = T)
  }

}

dev.off()
