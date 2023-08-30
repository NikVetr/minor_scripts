diff2tex <- function(target, current, 
                     latex_wrappers = list(target = "hlg", current = "hlr"),
                     neutral_terms = c("autocite")){
  # Get the diff in ansi8 format
  diff_out <- as.character(diffobj::diffChr(target = target, current = current, 
                                            format="ansi8", mode = "unified",
                                            context = -1L, disp.width	= 1E4))
  d <- paste0(diff_out[grepl("90m", diff_out)], collapse = "")
  #periods are handled differently from other punctuation :/
  
  # Remove extra characters and split into before and after
  #whoops this fails if tehre exist > in the text already
  # d <- gsub("\033\\[|39m|90m", "", d)
  # d <- gsub("33m< ", "", d)
  # d <- strsplit(d, ">")[[1]]
  
  #let's try this way... can prob be cleaned up but w/e
  d <- strsplit(d, "\033\\[34m>")[[1]]
  d <- sapply(d, gsub, pattern = "\033\\[|39m|90m", replacement = "")
  d <- sapply(d, gsub, pattern = "33m< ", replacement = "")
  
  d <- trimws(as.character(sapply(d, gsub, pattern = "34m$", replacement = "")))
  d <- unlist(lapply(d, strsplit, split = " "), recursive = F)
  
  #remove neutral term labels (cos it breaks the LaTeX)
  d <- lapply(d, function(x){
    nx <- x
    for(ni in neutral_terms){
      nx[grepl(ni, x)] <- gsub("33m|34m", "", nx[grepl(ni, x)])
    }
    nx
  }) 
  
  m <- lapply(d, function(x) rle(grepl("33m", x) - grepl("34m", x)))
  d <- lapply(d, gsub, pattern = "33m|34m", replacement = "")
  
  #process before and after
  d <- lapply(setNames(1:2, c("target", "current")), function(i){
    dn <- rep("", length(m[[i]]$lengths))
    l_inds <- rep(m[[i]]$lengths, m[[i]]$lengths)
    cumum <- cumsum(m[[i]]$lengths)
    concats <- cbind(cumum[m[[i]]$lengths != 1] - m[[i]]$lengths[m[[i]]$lengths != 1] + 1,
                     cumum[m[[i]]$lengths != 1])
    dn[m[[i]]$lengths == 1] <- d[[i]][l_inds == 1]
    dn[m[[i]]$lengths != 1] <- apply(concats, 1, function(x) paste0(d[[i]][x[1]:x[2]], collapse = " "))
    list(tokens = dn, code = m[[i]]$values)
  })
  
  #make LaTeX strings
  d$target$LaTeX[d$target$code == 1] <- paste0("\\", latex_wrappers$target, "{", d$target$tokens[d$target$code == 1], "}")
  d$target$LaTeX[d$target$code == 0] <- d$target$tokens[d$target$code == 0]
  d$target$LaTeX <- paste0(d$target$LaTeX, collapse = " ")
  
  d$current$LaTeX[d$current$code == -1] <- paste0("\\", latex_wrappers$current, "{", d$current$tokens[d$current$code == -1], "}")
  d$current$LaTeX[d$current$code == 0] <- d$current$tokens[d$current$code == 0]
  d$current$LaTeX <- paste0(d$current$LaTeX, collapse = " ")
  
  #return result
  list(target = d$target$LaTeX, current = d$current$LaTeX)
}

diff2tex(target = "this is the target sentence? this one too...", current = "and this one is current")

b4a <- readLines("~/Documents/Documents - nikolai/b4_n_after.txt", warn = F)
starts <- c(1,grep("smoldivider", b4a), length(b4a)+1)
out <- lapply(1:(length(starts) - 1), function(i){
  
  x <- b4a[starts[i]:(starts[i+1]-1)]
  
  orig <- x[grepl("Original:", x)]
  amen <- x[grepl("Amended:", x)]
  
  orig_txt <- trimws(gsub(".*\\\\begin\\{quote\\}", "", orig))
  amen_txt <- trimws(gsub(".*\\\\begin\\{quote\\}", "", amen))
  
  LaTeX_inserted <- diff2tex(target = amen_txt, current = orig_txt)
  
  list(gsub(orig_txt, LaTeX_inserted$current, x = orig, fixed = T),
       "\n",
       gsub(amen_txt, LaTeX_inserted$target, x = amen, fixed = T),
       "\n\n")
  
})
out <- unlist(out)

sink("~/Documents/Documents - nikolai/b4_n_after_revised.txt")
cat(out)
sink()

cat(out)
