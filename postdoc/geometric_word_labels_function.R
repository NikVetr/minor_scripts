xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

rgb2col <- function(x){ rgb(x[1,], x[2,], x[3,])}

priority_shift_vector <- function(vpw, eps = 1E-6){
  #problem is that pseq does staggered drop in priority even when drop is instant
  #need to look ahead of each initial drop to see how far to go
  rp_shift_index <-  c(TRUE, abs(diff(vpw)) >= eps)
  rlerp <- rle(rp_shift_index)
  rlerp <- data.frame(val = rlerp$values, len = rlerp$lengths)
  rlerp$incr <- rlerp$extra <- 0
  rlerp$incr[rlerp$val] <- 1 + rlerp$len[which(rlerp$val) + 1]
  rlerp$incr[is.na(rlerp$incr)] <- 1
  rlerp$extra[rlerp$val] <- rlerp$len[rlerp$val] - 1
  rlerp <- rlerp[rlerp$val,]
  repvec <- cumsum(rlerp$incr + rlerp$extra) #need to insert extras before cumsuming
  repvec <- Reduce(function(acc, i) 
    append(acc, cumsum(rep(1, rlerp$extra[i])) + ifelse(i==1, 0, repvec[i-1]), after = i-1), 
    seq_along(rlerp$extra), 
    init = repvec)
  reptimes <- Reduce(function(acc, i) 
    append(acc, rep(1, rlerp$extra[i]), after = i-1), 
    seq_along(rlerp$extra), 
    init = rlerp$incr)
  #adjust vector of shifts by half the width of the passing segments (because they've been shifted)
  repvec <- repvec - (reptimes - 1) / 2
  rp_incr_vec <- rep(repvec, times = reptimes)
  return(rp_incr_vec)
}

stacked_word_labels <- function(words, wcols, x0, y0, nro = 2, variable_priority = F, offx, 
                                ffamily = "Arial", eps = 1E-6, offset_y = 0, offset_x = 0, wsh_scale = 1){
  
  #get word count and correct for insufficiency
  nwords <- length(words)
  nro <- min(nro, nwords)
 
  #initialize word metadata
  dw <- data.frame(lab = words)
  dw$x0 <- x0
  dw$y0 <- y0
  dw <- dw[order(x0),]
  wcols <- wcols[order(x0)]
  wcols <- setNames(wcols, dw$lab)
  dw$row <- rep(1:nro, ceiling(nwords / nro))[1:nwords]
  midp <- mean(range(dw$x0)) + offset_x
  
  #continue word metadata
  dw$w <- strwidth(dw$lab, units = "user", family = ffamily)
  dw$h <- strheight(dw$lab, units = "user", family = ffamily)
  wsw <- strwidth("  ", units = "user", family = ffamily)
  wsh <- max(dw$h) * 1 * wsh_scale
  dw$y1 <- max(dw$y0) + diff(par("usr")[3:4])/5 + (dw$row - 1) * (max(dw$h) + wsh) + offset_y
  dw$x1 <- midp
  dw$rb <- dw$lb <- 0
  
  #handle single word case
  if(nwords == 1){
    text(x = dw$x0, y = dw$y1 + wsh, labels = dw$lab, family = ffamily, col = wcols)
    segments(x0 = dw$x0, y0 = dw$y0, x1 = dw$x0, y1 = dw$y1 - dw$h / 2 + wsh, col = wcols)
    return()
  }
  
  #set parameter for spacing between paths in vertical direction
  offy <- offx / xyrat()
  
  #find locations of words
  for(i in 1:nro){
    
    dws <- dw[dw$row == i,]
    
    #find bounds for words
    dws$lb <- dws$x1 - dws$w/2
    dws$rb <- dws$x1 + dws$w/2
    
    #stagger the words
    dws$x1 <- dws$x1 + cumsum(c(0, dws$w[-nrow(dws)]/2 + dws$w[-1]/2) + wsw)
    dws$lb <- dws$x1 - dws$w/2
    dws$rb <- dws$x1 + dws$w/2
    
    #center on midpoint of row
    shift_x <- - (max(dws$rb) + min(dws$lb)) / 2 + midp
    dws$x1 <- dws$x1 + shift_x
    dws$lb <- dws$lb + shift_x
    dws$rb <- dws$rb + shift_x
    
    #reassign
    dw[dw$row == i,] <- dws
  }
  
  #get line metadata
  row_inds <- setNames(1:nro, paste0("rt_", 1:nro))
  row_targets <- lapply(row_inds, function(ri){
    rt <- setNames(numeric(nwords), dw$lab)
    rt[dw$row < ri] <- NA
    rt[dw$row == ri] <- dw$x1[dw$row == ri]
    #retrieve midpoint of words' rb and lb one row below
    rt[dw$row > ri] <- dw$rb[which(dw$row > ri) - 
                               (dw$row[which(dw$row > ri)] - ri)] + wsw / 2
    return(rt)
  })
  
  #adjust spacing between words + row targets according to number of vertical and horiz strands
  row_targets <- lapply(setNames(seq_along(row_targets), names(row_targets)), function(ri){
    nrti <- rti <- row_targets[[ri]][!is.na(row_targets[[ri]])]
    trti <- table(rti)
    expand_locs <- as.numeric(names(trti[trti > 1]))
    expand_amounts <- (trti[trti > 1] - 1) * offx
    for(ei in seq_along(expand_locs)){
      nrti[rti < (expand_locs[ei] - eps)] <- nrti[rti < (expand_locs[ei] - eps)] - expand_amounts[ei] / 2
      nrti[rti > (expand_locs[ei] + eps)] <- nrti[rti > (expand_locs[ei] + eps)] + expand_amounts[ei] / 2
    }
    nrti <- nrti - mean(range(nrti)) + mean(range(rti))
    out <- row_targets[[ri]]
    out[names(nrti)] <- nrti
    return(out)
  })
  
  #update horizontal destinations in data frame
  dw$x1 <- do.call(rbind, row_targets)[cbind(dw$row, 1:nwords)]
  
  #establish priority queue to resolve overlaps
  #higher priority means closer to horiz_line_max
  row_priority_out <- lapply(row_inds, function(ri){
    if(ri > 1){
      rt0 <- setNames(row_targets[[ri-1]], dw$lab)
    } else {
      rt0 <- setNames(dw$x0, dw$lab)
    }
    rt1 <- setNames(row_targets[[ri]], dw$lab)
    
    #find partial overlaps
    omat <- (outer(rt0, rt1, "<=") & outer(rt1, rt0, ">=")) |
      (outer(rt1, rt0, "<=") & outer(rt0, rt1, ">="))
    colnames(omat) <- rownames(omat) <- dw$lab
    diag(omat) <- F
    omat[is.na(omat)] <- F
    ocounts <- rowSums(omat, na.rm = T)
    olaps_w <- apply(omat, 1, which, simplify = F)
    
    #find containment / full overlaps
    comat <- outer(rt0, rt0, "<=") & outer(rt1, rt1, ">=")
    diag(comat) <- F
    
    #find directions
    dirs <- sign(rt1 - rt0)
    opp_dirs <- outer(dirs, dirs, "*") == -1
    
    #check if completion is even possible
    all(!(opp_dirs & omat))
    all(!(comat))
    
    #the segment whose destination is closer to the midpoint 
    #needs to be on the inside and receive higher priority
    #when there are ties (segs heading to same target, eg a channel up), 
    #whoever started farther away gets priority
    #also, need to respect chains of priority / transitivity
    #if a < b and b < c then a < c
    
    #compare each pair of current entries
    #whichever one terminates first in the direction they're traveling gets higher priority
    rp <- sapply(1:nwords, function(i){
      sapply(1:nwords, function(j){
        
        if(i==j){return(F)}
        if(!omat[i,j]){return(F)}
        
        drt <- rt1[j] - rt1[i]
        direc <- dirs[i]
        
        #find the final row where at least one lives
        #and use that to determine local priority
        #when they part ways, whose destination is on the left vs right?
        #this determines priority at current level
        #(if in final destination A is negative / left wrt B, 
        #and traveling left currently, A has lower priority
        #if traveling right currently, A has higher priority
        
        incr <- 0
        while(abs(drt) < 1E-6){
          incr <- incr + 1
          rt1_incr <- setNames(row_targets[[ri+incr]], dw$lab)
          drt <- rt1_incr[j] - rt1_incr[i]
        }
        
        return(sign(drt) == direc)
        
      })
    })
    colnames(rp) <- rownames(rp) <- dw$lab
    
    #accommodate recursive priority (from staggered overlaps)
    #by treating as directed graph and computing transitive closure
    trans_rp <- ggm::transClos(rp)
    # trans_rp <- trans_rp * omat
    pr <- rowSums(trans_rp)
    
    #this should only apply at the point when a new strand joins, though
    #currently, it applies along the start of the strand
    
    #find where priority shifts occur due to higher priority strands terminating
    #basically, count the times along each strand (rt1 - rt0) when another rt1 is passed
    pr_shiftnames <- setNames(names(rt1[!is.na(rt1)]), names(rt1[!is.na(rt1)]))
    pr_shifts <- lapply(pr_shiftnames, function(tn){
      setd_rt1 <- rt1[setdiff(names(rt1), tn)]
      setd_rt1 <- setd_rt1[!is.na(setd_rt1)]
      rt_passed <- setd_rt1[setd_rt1 > rt0[tn] & setd_rt1 < rt1[tn] | 
                              setd_rt1 < rt0[tn] & setd_rt1 > rt1[tn]]
      return(rt_passed)
    })
    pr_shifts <- pr_shifts[sapply(pr_shifts, length) > 0]
    
    return(list(pr = pr, pr_shifts = pr_shifts))
    
  })
  row_priority <- lapply(row_priority_out, function(rp) rp[["pr"]])
  row_priority_shifts <- lapply(row_priority_out, function(rp) rp[["pr_shifts"]])
  
  
  #adjust vertical spacing to reflect number of horizontal rows
  max_lines_per_row <- sapply(row_priority, max)
  for(ri in 1:nro){
    dw$y1[dw$row >= ri] <- dw$y1[dw$row >= ri] + max(row_priority[[ri]]) * offy
  }
  
  #find top of horizontal lines for each row
  horiz_line_max <- sapply(1:nro, function(ri) 
    mean(dw$y1[dw$row == ri] - dw$h[dw$row == ri] * wsh_scale)) 
  
  #continue plotting to connect words with points using lines
  
  #plot the words
  text(x = dw$x1, y = dw$y1, labels = dw$lab, family = ffamily, col = wcols)
  
  #specify offsets
  dirs <- sign(row_targets[[1]] - dw$x0)
  x_offset <- row_priority[[1]] * dirs * offx
  y_offset <- row_priority[[1]] * offy
  
  #re-center x offsets on target
  offset_groups <- split(names(row_targets[[1]]), as.character(row_targets[[1]]))
  names(offset_groups) <- NULL
  x_offset_midp <- unlist(lapply(offset_groups, function(x) {
    setNames(rep(mean(range(x_offset[x])), length(x)), x)
  }))
  x_offset[names(x_offset_midp)] <- x_offset[names(x_offset_midp)] - x_offset_midp
  
  #ignore x-offsets if terminating here
  twords <- dw$row == 1
  x_offset[twords] <- 0
  
  #vertical line from point
  segments(x0 = dw$x0, 
           x1 = dw$x0, 
           y0 = dw$y0, 
           y1 = horiz_line_max[1] - y_offset,
           col = wcols)
  
  if(variable_priority){
    
    #parse words with variable priority
    vp <- row_priority_shifts[[1]] #variable priority
    sp <- !(dw$lab %in% names(vp)) #static priority
    stp <- fip <- row_priority[[1]]
    if(length(vp) > 0){
      fip[names(vp)] <- fip[names(vp)] - sapply(vp, length)  
    }  
    #horizontal line to word target or passageway
    segments(x0 = dw$x0[sp], 
             x1 = (row_targets[[1]] + x_offset)[sp], 
             y0 = (horiz_line_max[1] - y_offset)[sp], 
             y1 = (horiz_line_max[1] - y_offset)[sp],
             col = wcols[sp])
    
    #vertical line up to the terminating words
    segments(x0 = row_targets[[1]][twords & sp] + x_offset[twords & sp], 
             x1 = row_targets[[1]][twords & sp] + x_offset[twords & sp], 
             y0 = horiz_line_max[1] - y_offset[twords & sp], 
             y1 = dw$y1[twords & sp] - dw$h[twords & sp] / 2,
             col = wcols[twords & sp])
    
    for(wi in names(vp)){
      # wi = "qnM8"
      
      #segment info
      vpw <- sort(vp[[wi]])
      rp_incr_vec <- priority_shift_vector(vpw)
      
      pseq <- row_priority[[1]][wi] - rp_incr_vec + 1
      if(dirs[wi] == -1){
        vpw <- rev(vpw)
      }
      ns <- length(vpw)
      corner_offset <- pseq * offx * dirs[wi]
      
      #horizontal connecting lines
      segments(x0 = c(setNames(dw$x0, dw$lab)[wi], vpw + corner_offset), 
               x1 = c(vpw + corner_offset, row_targets[[1]][wi] + x_offset[wi]), 
               y0 = horiz_line_max[1] - y_offset[wi] + 
                 c(0, seq_along(vpw)) * offy, 
               y1 = horiz_line_max[1] - y_offset[wi] + 
                 c(0, seq_along(vpw)) * offy,
               col = wcols[wi])
      
      #vertical connecting lines
      segments(x0 = vpw + corner_offset,
               x1 = vpw + corner_offset,
               y0 = horiz_line_max[1] - y_offset[wi] +
                 seq_along(vpw) * offy,
               y1 = horiz_line_max[1] - y_offset[wi] +
                 c(0, seq_along(vpw))[-(ns+1)] * offy,
               col = wcols[wi])
      
      #vertical line up to the terminating words
      segments(x0 = row_targets[[1]][twords] + x_offset[twords], 
               x1 = row_targets[[1]][twords] + x_offset[twords], 
               y0 = horiz_line_max[1] - fip[twords] * offy, 
               y1 = dw$y1[twords] - dw$h[twords] / 2,
               col = wcols[twords])
      
    }
    
    #record previous values for offsets
    x_offset_prev <- x_offset
    y_offset_prev <- fip * offy
    
  } else {
    
    #horizontal line to word target or passageway
    segments(x0 = dw$x0, 
             x1 = row_targets[[1]] + x_offset, 
             y0 = horiz_line_max[1] - y_offset, 
             y1 = horiz_line_max[1] - y_offset,
             col = wcols)
    
    #vertical line up to the terminating words
    segments(x0 = row_targets[[1]][twords] + x_offset[twords], 
             x1 = row_targets[[1]][twords] + x_offset[twords], 
             y0 = horiz_line_max[1] - y_offset[twords], 
             y1 = dw$y1[twords] - dw$h[twords] / 2,
             col = wcols[twords])
    
    #record previous values for offsets
    x_offset_prev <- x_offset
    y_offset_prev <- y_offset
  }
  
  ####
  
  
  #iterate through rows connecting starts and ends
  for(ri in 2:nro){
    
    #which words need lines drawn?
    rt_val <- !is.na(row_targets[[ri]])# & row_priority[[ri]] == 0
    
    #compute offsets
    dirs <- sign(row_targets[[ri]] - row_targets[[ri-1]])
    x_offset <-  row_priority[[ri]] * dirs * offx
    y_offset <-  row_priority[[ri]] * offy
    
    #re-center x offsets on target
    offset_groups <- split(names(row_targets[[ri]]), as.character(row_targets[[ri]]))
    names(offset_groups) <- NULL
    x_offset_midp <- unlist(lapply(offset_groups, function(x) {
      setNames(rep(mean(range(x_offset[x])), length(x)), x)
    }))
    x_offset[names(x_offset_midp)] <- x_offset[names(x_offset_midp)] - x_offset_midp
    
    #ignore x-offsets if terminating here
    twords <- dw$row == ri
    x_offset[twords] <- 0
    
    if(variable_priority){
      
      #parse words with variable priority
      vp <- row_priority_shifts[[ri]] #variable priority
      sp <- !(dw$lab %in% names(vp)) #static priority
      stp <- fip <- row_priority[[ri]] #starting and final priority
      if(length(vp) > 0){
        fip[names(vp)] <- fip[names(vp)] - sapply(vp, length)  
      }
      
      #vertical lines to next level
      segments(x0 = row_targets[[ri-1]][rt_val] + x_offset_prev[rt_val],
               x1 = row_targets[[ri-1]][rt_val] + x_offset_prev[rt_val],
               y0 = horiz_line_max[ri-1] - y_offset_prev[rt_val],
               y1 = horiz_line_max[ri] - y_offset[rt_val],
               col = wcols[rt_val])
      
      #horizontal lines to next row target
      segments(x0 = row_targets[[ri-1]][rt_val & sp] + x_offset_prev[rt_val & sp], 
               x1 = row_targets[[ri]][rt_val & sp] + x_offset[rt_val & sp], 
               y0 = horiz_line_max[ri] - y_offset[rt_val & sp], 
               y1 = horiz_line_max[ri] - y_offset[rt_val & sp],
               col = wcols[rt_val & sp])
      
      #vertical line up to the terminating words
      segments(x0 = row_targets[[ri]][twords & sp] + x_offset[twords & sp], 
               x1 = row_targets[[ri]][twords & sp] + x_offset[twords & sp], 
               y0 = horiz_line_max[ri] - y_offset[twords & sp], 
               y1 = dw$y1[twords & sp] - dw$h[twords & sp] / 2,
               col = wcols[twords & sp])
      
      for(wi in names(vp)){
        
        #segment info
        vpw <- sort(vp[[wi]])
        rp_incr_vec <- priority_shift_vector(vpw)
        pseq <- row_priority[[ri]][wi] - rp_incr_vec + 1
        if(dirs[wi] == -1){
          vpw <- rev(vpw)
        }
        ns <- length(vpw)
        corner_offset <- pseq *  offx * dirs[wi]
        
        #horizontal connecting lines
        segments(x0 = c(row_targets[[ri-1]][wi] + x_offset_prev[wi], vpw + corner_offset), 
                 x1 = c(vpw + corner_offset, row_targets[[ri]][wi] + x_offset[wi]), 
                 y0 = horiz_line_max[ri] - y_offset[wi] + 
                   c(0, seq_along(vpw)) * offy, 
                 y1 = horiz_line_max[ri] - y_offset[wi] + 
                   c(0, seq_along(vpw)) * offy,
                 col = wcols[wi])
        
        #vertical connecting lines
        segments(x0 = vpw + corner_offset,
                 x1 = vpw + corner_offset,
                 y0 = horiz_line_max[ri] - y_offset[wi] +
                   seq_along(vpw) * offy,
                 y1 = horiz_line_max[ri] - y_offset[wi] +
                   c(0, seq_along(vpw))[-(ns+1)] * offy,
                 col = wcols[wi])
        
        #vertical line up to the terminating words
        segments(x0 = row_targets[[ri]][twords] + x_offset[twords], 
                 x1 = row_targets[[ri]][twords] + x_offset[twords], 
                 y0 = horiz_line_max[ri] - fip[twords] * offy, 
                 y1 = dw$y1[twords] - dw$h[twords] / 2,
                 col = wcols[twords])
        
      }
      
      #retain the last round's info
      x_offset_prev <- x_offset
      y_offset_prev <- fip * offy
      
    } else {
      
      #vertical lines to next level
      segments(x0 = row_targets[[ri-1]][rt_val] + x_offset_prev[rt_val],
               x1 = row_targets[[ri-1]][rt_val] + x_offset_prev[rt_val],
               y0 = horiz_line_max[ri-1] - y_offset_prev[rt_val],
               y1 = horiz_line_max[ri] - y_offset[rt_val],
               col = wcols[rt_val])
      
      #horizontal lines to next row target
      segments(x0 = row_targets[[ri-1]][rt_val] + x_offset_prev[rt_val], 
               x1 = row_targets[[ri]][rt_val] + x_offset[rt_val], 
               y0 = horiz_line_max[ri] - y_offset[rt_val], 
               y1 = horiz_line_max[ri] - y_offset[rt_val],
               col = wcols[rt_val])
      
      #vertical line up to the terminating words
      segments(x0 = row_targets[[ri]][twords] + x_offset[twords], 
               x1 = row_targets[[ri]][twords] + x_offset[twords], 
               y0 = horiz_line_max[ri] - y_offset[twords], 
               y1 = dw$y1[twords] - dw$h[twords] / 2,
               col = wcols[twords])
      
      #retain the last round's info
      x_offset_prev <- x_offset
      y_offset_prev <- y_offset
      
    }
    
  }
    
}

