labs <- "\\textbf{$h^2_{SNP}$ Estimation}\nQuantify both\nSNP heritability\nenrichment and\nexpression $h_{SNP}$"
labs <- "\\textbf{GWAS}\n114 Genome-\nWide Association\nStudies across 12\ntrait categories"
# labs <- paste0(unlist(strsplit("SNP heritability", "")), "\n", collapse = "")


xyrat <- function(){
  prop <- c(diff(par("usr")[1:2]), diff(par("usr")[3:4])) / par("pin")
  prop[1] / prop[2]
}

polar2cart <- function(t, r){
  return(c(r*cos(t), r * sin(t)))
}

text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(labels, cex = cex) / 2, 
         xright = adj_x + strwidth(labels, cex = cex) / 2, 
         ybottom = adj_y - strheight(labels, cex = cex) / 2, 
         ytop = adj_y + strheight(labels, cex = cex) / 2)
  }
}

latext <- function(labs, x, y, cex = 1, LETTing = T, boxh = NA, first_line_col = 1, first_line_hadj = NA, col = 1, pos = NULL, ...){
  new_labs <- strsplit(labs, split = "\n")[[1]]
  new_labs_no_uflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u", "a", l)))
  new_labs_no_oflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("\\^", "a", l)))
  new_labs <- lapply(new_labs, function(l) latex2exp::TeX(l))
  
  wsh <- strheight("M\nM", cex = cex) - strheight("M", cex = cex) * 2
  lineh <- sapply(new_labs, strheight, cex = cex)  
  lineh_no_uflow <- sapply(new_labs_no_uflow, strheight, cex = cex)
  lineh_no_oflow <- sapply(new_labs_no_oflow, strheight, cex = cex)
  ebot <- lineh_no_uflow - lineh
  etop <- lineh_no_oflow - lineh
  uflow_adj <- (lineh_no_uflow - lineh) / 2 - (lineh - lineh_no_oflow) / 2
  uflow_adj <- (lineh_no_uflow - lineh)
  flow_adj <- ebot
  flow_adj[etop < -1E-6] <- flow_adj[etop < -1E-6] / 2 + etop[etop < -1E-6] / 2
  
  charh <- rep(strheight("A", cex = cex), length(new_labs)-1)
  yadj <- -cumsum(c(0, charh + wsh)) + lineh/2
  
  if(!is.na(boxh)){
    yadj[-1] <- yadj[-1] - (boxh + tail(yadj,1) + yadj[2]) / 5
  }
  
  for(i in 1:length(new_labs)){
    # print(ifelse(i==1 & !is.na(first_line_hadj), paste0(new_labs[[i]], ": ", (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2), 0))
    text(x = x + ifelse(i==1 & !is.na(first_line_hadj), (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2, 0), 
         y = y + yadj[i] + flow_adj[i], labels = new_labs[[i]], cex = cex, col = ifelse(i==1, first_line_col, col),
         pos = ifelse(i==1 & !is.na(first_line_hadj), NULL, pos), ...)
    # points(x = x, y = y+yadj[i]-lineh[i]/2)
    # segments(x0 = x, x1 = x + strwidth(new_labs[[i]], cex = cex) * 1.5, y0 = y+yadj[i]-lineh[i]/2+0.0, y1 = y+yadj[i]-lineh[i]/2+0.0)
    # rect(xleft = x, xright = x + strwidth(new_labs[[i]], cex = cex), ybottom = y+yadj[i]-lineh[i]/2, ytop = y+yadj[i]+lineh[i]/2)
  }
}

rrect <- function(loc, w, h, pe = 0.25, npts = 50, rot = 0, hat_prop = 0.15, bold_border = 0,
                  col = 0, lwd = 1, border = 1, background_col = 0, ...){

  #get arc where yval has height 1, xval has height xyrat
  xyr <- xyrat()
  tr <- polar2cart(t = seq(0, pi/2, length.out = round(npts/4)), r = 1) %*% diag(c(xyr, 1))
  
  #scale to smaller edge
  tallboi <- w < (h * xyr)
  tr <- tr * ifelse(tallboi, pe * w / 2 / xyr, pe * h / 2)
  
  #shift to circumscribed centered rectangle quadrant
  tr[,1] <- (tr[,1] + w / 2 - tr[1,1])
  tr[,2] <- (tr[,2] + h / 2 - tr[nrow(tr),2])
  # rect(xleft = loc[1] - w/2, ybottom = loc[2] - h/2, xright = loc[1] + w/2, ytop = loc[2] + h/2)
  
  #find polygon points and recenter
  outer_coords_nc <- rbind(tr, cbind(-rev(tr[,1]), rev(tr[,2])), cbind(-tr[,1], -tr[,2]), cbind(rev(tr[,1]), -rev(tr[,2])))
  outer_coords <- outer_coords_nc + rep(loc, each = nrow(outer_coords_nc))

  #first draw outer border if desired & inner polygon
  if(bold_border > 1E-6){
    
    #draw solid outer polygon
    polygon(outer_coords, col = border, border = border, lwd = 1E-6)
    
    #adjust border to bounds
    bold_border <- max(c(min(c(bold_border, 1)), 0))
    
    if(tallboi){
      nw <- (1-bold_border) * w
      nh <- h - bold_border * w / xyr
    } else { #shortboi
      nw <- w - bold_border * h * xyr
      nh <- (1-bold_border) * h
    }
    
    tr <- polar2cart(t = seq(0, pi/2, length.out = round(npts/4)), r = 1) %*% diag(c(xyr, 1))
    tr <- tr * ifelse(tallboi, pe * nw / 2 / xyr, pe * nh / 2)
    tr[,1] <- (tr[,1] + nw / 2 - tr[1,1])
    tr[,2] <- (tr[,2] + nh / 2 - tr[nrow(tr),2])
    
    inner_coords_nc <- rbind(tr, cbind(-rev(tr[,1]), rev(tr[,2])), cbind(-tr[,1], -tr[,2]), cbind(rev(tr[,1]), -rev(tr[,2])))
    inner_coords <- inner_coords_nc + rep(loc, each = nrow(inner_coords_nc))
    # rect(xleft = loc[1] - nw/2, ybottom = loc[2] - nh/2, xright = loc[1] + nw/2, ytop = loc[2] + nh/2)

    # inner_coords <- outer_coords_nc %*% diag(ifelse2(tallboi,
    #                                                  c(1-bold_border, 1 - bold_border * w * xyr / h),
    #                                                  c(1 - bold_border * h / xyr / w, 1-bold_border))
    #                                          ) +
    #   rep(loc, each = nrow(inner_coords))

    polygon(inner_coords, col = background_col, lwd = 1E-6, border = background_col)
    polygon(inner_coords, col = col, lwd = 1E-6, border = border)
  } else {
    polygon(outer_coords, col = col, lwd = lwd, border = border)
  }

  #add a hat if desired
  if(hat_prop > 1E-6){
    sub_coords <- outer_coords[outer_coords[,2] > loc[2] - h/2 + h*(1-hat_prop),]
    if(bold_border > 1E-6){
      polygon(sub_coords, col = background_col, lwd = lwd, border = background_col)
      polygon(sub_coords, col = col, lwd = 1E-6, border = border)
    } else {
      polygon(sub_coords, col = col, lwd = lwd, border = border)  
    }
  }
  
}


# plot(NULL, xlim = c(-0.125,9), ylim = c(0,3.5))
# rrect(loc = c(4,2), w = 1, h = 2, border = adjustcolor(2, 1), col = adjustcolor(2, 0.1),
#       lwd = 1, pe = 1, hat_prop = 0, bold_border = 0)

loc = c(2,4)
w = 3
h = 6
plot(NULL, xlim = c(0,5), ylim = c(0,8))
rrect(loc = loc, w = w, h = h, border = adjustcolor(2, 1), col = adjustcolor(2, 0.1),
      lwd = 2, pe = 0.25, hat_prop = 0, bold_border = 0.1)

latext2 <- function(labs, x, y, cex = 1, LETTing = T, boxh = NA, first_line_col = 1, first_line_hadj = NA, col = 1, pos = NULL, ...){
  new_labs <- strsplit(labs, split = "\n")[[1]]
  new_labs_no_uflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("g|j|p|q|y|,|_|\\(|\\)|Q|u", "a", l)))
  new_labs_no_oflow <- lapply(new_labs, function(l) latex2exp::TeX(gsub("\\^", "a", l)))
  new_labs <- lapply(new_labs, function(l) latex2exp::TeX(l))
  
  wsh <- strheight("M\nM", cex = cex) - strheight("M", cex = cex) * 2
  lineh <- sapply(new_labs, strheight, cex = cex)  
  lineh_no_uflow <- sapply(new_labs_no_uflow, strheight, cex = cex)
  lineh_no_oflow <- sapply(new_labs_no_oflow, strheight, cex = cex)
  ebot <- lineh_no_uflow - lineh
  etop <- lineh_no_oflow - lineh
  flow_adj <- ebot
  flow_adj[etop < -1E-6] <- flow_adj[etop < -1E-6] / 2 + etop[etop < -1E-6] / 2
  
  charh <- rep(strheight("A", cex = cex), length(new_labs)-1)
  yadj <- -cumsum(c(0, charh + wsh)) + lineh/2
  
  if(!is.na(boxh)){
    yadj[-1] <- yadj[-1] - (boxh + tail(yadj,1) + yadj[2]) / 5
  }
  
  for(i in 1:length(new_labs)){
    # print(ifelse(i==1 & !is.na(first_line_hadj), paste0(new_labs[[i]], ": ", (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2), 0))
    text2(x = x + ifelse(i==1 & !is.na(first_line_hadj), (first_line_hadj - strwidth(new_labs[[i]], cex = cex)) / 2, 0), 
         y = y + yadj[i] + flow_adj[i], labels = new_labs[[i]], cex = cex, col = ifelse(i==1, first_line_col, col),
         pos = pos, ...)
    # points(x = x, y = y+yadj[i]-lineh[i]/2)
    # segments(x0 = x, x1 = x + strwidth(new_labs[[i]], cex = cex) * 1.5, y0 = y+yadj[i]-lineh[i]/2+0.0, y1 = y+yadj[i]-lineh[i]/2+0.0)
    # rect(xleft = x, xright = x + strwidth(new_labs[[i]], cex = cex), ybottom = y+yadj[i]-lineh[i]/2, ytop = y+yadj[i]+lineh[i]/2)
  }
}

points(loc[1], loc[2], pch = 19)
latext2(labs = labs,
       x = loc[1],
       y = loc[2],
       pos = 3, cex = text_cex[i], LETTing = T,
       boxh = ifelse(cols[i] == "#FFFFFF", NA, hs[i]),
       first_line_col = ifelse(cols[i] == "#FFFFFF", 1, cols[i]),
       first_line_hadj = NA, drect = T)



