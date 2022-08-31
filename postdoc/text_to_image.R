#### load libraries ####
library(png)

#functions
pythag <- function(x) sqrt(sum(x^2))

#read image and text
img_name <- c("motrpac_M", "rat_anatomy", "treadmill_rat", "treadmill_rat_horizontal", "rethink_logo", "treadmill_rat_horizontal_straightuptop",
              "love-song", "treadmill_rat_horizontal_straightuptop_thickerlines")[5]
path_to_image <- paste0("~/Pictures/", img_name, ".png")
img <- readPNG(source = path_to_image)

#add transparent layer if none exists
if(dim(img)[3] == 3){
  transparent_color <- "black"
  tc_dist <- sqrt(3) / 3
  tprgb <- c(col2rgb(transparent_color))
  pixrgb <- cbind(c(img[,,1]), c(img[,,2]), c(img[,,3]))
  distrgb <- apply(t(t(pixrgb) - tprgb), 1, pythag)
  img <- abind::abind(img, img[,,3], along = 3)
  img[,,4] <- as.numeric(distrgb > tc_dist)
}

img[,,4] <- img[,,4]^(0.5)

#get text
# txt <- sample(c(letters, LETTERS, rep(" ", 10)), 1E4, replace = T)
# txt <-  c(rbind(strsplit(paste0(1:8E2, " ", collapse = ""), " ")[[1]], " "))
# txt <- readLines("~/Documents/Documents - nikolai/rp_intro.txt", warn = F)
txt_name <- c("pass1b_landscape", "nicole_thesis", "nicole_thesis_short", "nicole_thesis_extrashort", "love-song", "rp_intro")[6]
path_to_txt <- paste0("~/Documents/Documents - nikolai/", txt_name, ".txt")
txt <- readLines(path_to_txt)
txt <- paste0(txt, collapse = " ")
txt <- gsub(txt, pattern = "\t", replacement = " ")
txt <- gsub(txt, pattern = "\\s+", replacement = " ")
grep("\t", txt)
grep("  ", txt)

txt <- strsplit(txt, "")[[1]]
# txt <- paste0(strsplit(txt, " ")[[1]])
# txt <-  c(rbind(txt, " "))
# txt <- c("Rethink Priorities", txt[-(1:nchar("Rethink Priorities"))])
# txt <- c(txt, txt)
# txt <- txt[1:1E3]
# txt <- strsplit(paste0(1:1E3, " ", collapse = ""), "")[[1]]

#check text for illegal characters & evaluate fill proportion
font = 2
family = "Arial"
uniqchars <- unique(txt)
w2h <- c(strwidth(paste0(uniqchars, collapse = ""), font = font, family = family, units = "inches"),
         strheight(paste0(uniqchars, collapse = ""), font = font, family = family, units = "inches"))
w2h <- w2h / min(w2h) * 100
png(paste0("~/uniqchars-", txt_name,".png"), width = w2h[1]*0.5, height = w2h[2], units = "px")
par(mar = c(0,0,0,0))
plot(NA, NA, xlim = c(0,w2h[1] / 100), ylim = c(-0.5,1.5), frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
uniqcharcex <- 1 / strheight(uniqchars[1], units = "user", font = font, family = family)
uniqcharwidths <- strwidth(uniqchars, units = "user", cex = uniqcharcex, font = font, family = family)
cs_uniqcharwidths <- cumsum(c(0,uniqcharwidths))
# text(labels = paste0(uniqchars, collapse = ""), x = cs_uniqcharwidths[1], y = 0.5,
#      col = 1, cex = cex, font = font, pos = 4, family = family)
for(i in 1:length(uniqchars)){
  text(labels = uniqchars[i], x = cs_uniqcharwidths[i], y = 0.5,
       col = 1, cex = uniqcharcex, font = font, pos = 4, family = family)
  rect(xleft = cs_uniqcharwidths[i], ybottom = 0, ytop = 1, xright = cs_uniqcharwidths[i] + uniqcharwidths[i], border = 1)
}
dev.off()

# 
# ucimg <- readPNG("~/uniqchars.png")
# file.remove("~/uniqchars.png")
# ucimghits <- ucimg[,,1] + ucimg[,,2] + ucimg[,,3]
# ucimghits_i <- which((abs(ucimghits - 3) > 1E-9) & ucimg[,,4] > 0.01) 
# ucimghits[ucimghits_i] <- 1
# ucimghits[-ucimghits_i] <- 0
# ucimghits_indloc <- as.data.frame(which(ucimghits == 1, arr.ind = T))
# bottom_row <- max(ucimghits_indloc$row[ucimghits_indloc$col == ucimghits_indloc$col[1]])
# ucimghits_indloc$col[ucimghits_indloc$row == (bottom_row-2)]
# 
# cs_uniqcharwidths_pix <- c(round(cs_uniqcharwidths * 100), ncol(ucimghits))
# propspace_uniqchar <- sapply(1:length(uniqchars), function(i) mean(ucimghits[,cs_uniqcharwidths_pix[i]:cs_uniqcharwidths_pix[i+1]]))

get_char_propspace <- function(uc, font, family){
  png(paste0("~/tmp.png"), width = 1000, height = 250, units = "px")
  par(mar = c(0,0,0,0))
  plot(NA, NA, xlim = c(0,4), ylim = c(0,1), frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  rect(xleft = -5, ybottom = -5, ytop = 5, xright = 5, border = 1, col = 1)
  uccex <- 1 / strheight(uc, units = "user", font = font, family = family) * 0.5
  text(x = 0, y = 0.5, pos = 4, cex = uccex, col = rgb(0,0,1,1), labels = uc)
  boxw <- strwidth(uc, units = "user", cex = uccex, font = font, family = family)
  boxh <- strheight(uc, units = "user", cex = uccex, font = font, family = family)
  rect(xleft = 2, ybottom = 0, ytop = 0 + boxh, xright = 2 + boxw, border = rgb(1,0,0,1), col = rgb(1,0,0,1))
  dev.off()
  ucmg <- png::readPNG("~/tmp.png")
  file.remove("~/tmp.png")
  sum(ucmg[,,3] > 0.9) / sum(ucmg[,,1] > 0.9)
}

uniqchars <- unique(txt)
propspace_uniqchar <- sapply(uniqchars, function(uc) get_char_propspace(uc, font = font, family = family))

#### rotate image if desired ###
rotation_angle <- 0
img <- EBImage::rotate(img, -rotation_angle) #let's just do it with an existing package :S
rotmat_00 <- function(t){r <- t / 360 * 2 * pi; matrix(c(cos(r), -sin(r), sin(r), cos(r)), 2, 2, byrow = T)}
# interpolate <- T
# if(abs(rotation_angle - 0) > 1E-6){
#   diagpix <- ceiling(sqrt(ncol(img)^2 + nrow(img)^2))
#   rotimg <- array(0, dim = c(diagpix,diagpix,4))
#   orig_pixcoords <- as.data.frame(expand.grid(1:nrow(img), 1:ncol(img)))
#   colnames(orig_pixcoords) <- c("row", "col")
#   orig_pixcoords$row <- nrow(img) - orig_pixcoords$row + 1
#   pixcoords <- t(t(orig_pixcoords) - c(nrow(img), ncol(img)) / 2)
#   
#   rotmat <- rotmat_00(rotation_angle)
#   new_pixcoords <- unrounded_pix_coords <- t(rotmat %*% t(pixcoords))
#   new_pixcoords <- round(t(t(new_pixcoords) + c(diagpix, diagpix) / 2))
#   
#   red <- img[,,1]; newred <- matrix(0, diagpix, diagpix)
#   green <- img[,,2]; newgreen <- matrix(0, diagpix, diagpix)
#   blue <- img[,,3]; newblue <- matrix(0, diagpix, diagpix)
#   alpha <- img[,,4]; newalpha <- matrix(0, diagpix, diagpix)
#   
#   #2x2 interpolation
#   if(interpolate){
#     alt_pixcoords <- unique(rbind(
#       cbind(ceiling(unrounded_pix_coords[,1]), ceiling(unrounded_pix_coords[,2])),
#       cbind(ceiling(unrounded_pix_coords[,1]), floor(unrounded_pix_coords[,2])),
#       cbind(floor(unrounded_pix_coords[,1]), ceiling(unrounded_pix_coords[,2])),
#       cbind(floor(unrounded_pix_coords[,1]), floor(unrounded_pix_coords[,2]))
#     ))
#     matching_alt_pixcoords <- t(rotmat_00(-rotation_angle) %*% t(alt_pixcoords))
#     rmap <- round(matching_alt_pixcoords)
#     indmatch <- match(apply(rmap, 1, paste0, collapse = ""), apply(pixcoords, 1, paste0, collapse = ""))
#     alt_pixcoords <- alt_pixcoords[-is.na(indmatch),]
#     indmatch <- indmatch[-is.na(indmatch)]
#     alt_pixcoords <- round(t(t(alt_pixcoords) + c(diagpix, diagpix) / 2))
#     
#     newred[as.matrix(alt_pixcoords)] <- red[as.matrix(orig_pixcoords[indmatch,])]
#     newgreen[as.matrix(alt_pixcoords)] <- green[as.matrix(orig_pixcoords[indmatch,])]
#     newblue[as.matrix(alt_pixcoords)] <- blue[as.matrix(orig_pixcoords[indmatch,])]
#     newalpha[as.matrix(alt_pixcoords)] <- alpha[as.matrix(orig_pixcoords[indmatch,])]
#     
#   } else {
#     # new_pixcoords[,1] <- diagpix - new_pixcoords[,1] + 1
#     newred[as.matrix(new_pixcoords)] <- red[as.matrix(orig_pixcoords)]
#     newgreen[as.matrix(new_pixcoords)] <- green[as.matrix(orig_pixcoords)]
#     newblue[as.matrix(new_pixcoords)] <- blue[as.matrix(orig_pixcoords)]
#     newalpha[as.matrix(new_pixcoords)] <- alpha[as.matrix(orig_pixcoords)]
# 
#   }
# 
#   rotimg[,,1] <- newred
#   rotimg[,,2] <- newgreen
#   rotimg[,,3] <- newblue
#   rotimg[,,4] <- newalpha
#   img <- rotimg
#   
# }

#get hits matrix
hits <- img[,,1] + img[,,2] + img[,,3]
hits_i <- which((abs(hits - 3) > 1E-9) & img[,,4] > 0.2) 
hits[hits_i] <- 1
hits[-hits_i] <- 0

#do easy fill
# sf <- sqrt(sum(hits) / length(txt))
# subimg <- img[round(seq(1, nrow(img), by = sf)), round(seq(1, ncol(img), by = sf)),]
# subhits <- hits[round(seq(1, nrow(hits), by = sf)), round(seq(1, ncol(hits), by = sf))]
# indsubhits <- which(subhits == 1, arr.ind = T)
# indsubhits <- indsubhits[order(indsubhits[,1]),]
# cols <- t(sapply(1:nrow(indsubhits), function(i) subimg[indsubhits[i,1], indsubhits[i,2],]))
# cols <- sapply(1:nrow(cols), function(i) rgb(cols[i,1], cols[i,2], cols[i,3], cols[i,4]))

# # paint image
# png("~/Pictures/out.png", width = ncol(hits) * 5, height = nrow(hits) * 5, units = "px")
# plot(NA, NA, xlim = c(1,ncol(subhits)), ylim = c(1,nrow(subhits)),
#      frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# for(i in 1:length(cols)){
#   if(i %% round(length(cols) / 100) == 0) cat(i / round(length(cols) / 100), " ")
#   text(labels = txt[i], y = nrow(subhits) - indsubhits[i,1], x = indsubhits[i,2],
#        col = cols[i], cex = 1, font = 2)
# }
# dev.off()

#### do harder fill
pxlw <- 1 / ncol(hits)
pxlh <- 1 / nrow(hits)
indhits <- which(hits == 1, arr.ind = T)
indhits <- indhits[order(indhits[,1]),]
indlocs <- data.frame(x = indhits[,2] / ncol(hits), y = 1 - indhits[,1] / nrow(hits))

indlocs$tsw <- NA
indlocs$bi <- NA
rlepixy <- rle(indlocs$y)
rlepixy_locs <- cumsum(c(1,rlepixy$lengths))

# sum_consec_ones <- function(x){
#   if(!any(x == 1)){return(x)}
#   if(!any(x != 1)){return(NA)}
#   inds_1 <- which(x == 1)
#   inds_n1 <- which(x != 1)
#   if(length(inds_1) == 1){
#     if(inds_1 != length(x)){
#       x[inds_1+1] <- x[inds_1+1] + 1  
#     } else {
#       x[inds_1-1] <- x[inds_1-1] + 1
#     }
#     x <- x[-inds_1]
#     return(x)
#   }
#   matches1n1 <- which(sapply(inds_n1, function(inds_n1i) inds_n1i - inds_1) == 1, arr.ind = T)
#   x[inds_n1[matches1n1[,2]]] <- x[inds_n1[matches1n1[,2]]] + 1
#   x <- x[-inds_1[matches1n1[,1]]]
#   if(any(x[1:max(which(x != 1))] == 1)){
#     x <- sum_consec_ones(x)
#   }
#   maxn1i <- max(which(x != 1))
#   x[maxn1i] <- x[maxn1i] + (length(x) - maxn1i)
#   x <- x[1:maxn1i]
#   x
# }

for(i in (1:length(rlepixy$values))){
  curr_inds <- ifelse(i==1,1,rlepixy_locs[i]):(rlepixy_locs[i+1]-1)
  # rlepixsuby <- rle(c(diff(indlocs[curr_inds,]$x, pxlw)) - pxlw < (1*pxlw))
  # rlengths <- rlepixsuby$lengths
  # #adjust rlengths to not count the skip as a run
  # #when there are lots of '1's in a row, add them up and put them in the next non-1 space
  # rlengths <- sum_consec_ones(rlengths)
  # twosuby <- rlengths * pxlw #total width on each subrow 'y'
  # # print(i)
  # 
  #bah can't figure out this rle thing :/
  
  v = indlocs[curr_inds,]$x
  vb <- rep(NA, length(v))
  vbi <- rep(NA, length(v))
  for(vi in 1:length(v)){
    if(vi == 1){
      cvl <- 1
      bi <- 1
    } else {
      dv <- (v[vi] - v[vi-1] - pxlw) < 1E-6
      if(dv){
        cvl <- cvl + 1
        if(vi == length(v)){
          vb[(vi):(vi-cvl+1)] <- cvl
          vbi[(vi):(vi-cvl+1)] <- bi
        }
      } else {
        vb[(vi-1):(vi-cvl)] <- cvl
        vbi[(vi-1):(vi-cvl)] <- bi
        cvl <- 1
        bi <- bi + 1
        if(vi == length(v)){
          vb[vi] <- cvl
          vbi[vi] <- bi
        }
      }
    }
  }
  
  indlocs$tsw[curr_inds] <- vb * pxlw
  indlocs$bi[curr_inds] <- vbi
  
  # unlist(lapply(1:length(rlepixsuby$values), function(rlei){
  #   if(rlepixsuby$values[rlei]){
  #     rep(rlepixsuby$lengths[rlei], rlepixsuby$lengths[rlei])
  #   } else {
  #     rep(1, rlepixsuby$lengths[rlei])
  #   }
  # })) * pxlw
  # 
  # if(i != length(rlepixy$values)){
  #   indlocs$tsw[curr_inds] <- rep(twosuby, rlengths)  
  #   if(length(curr_inds) != length(rep(twosuby, rlengths))){print(i); print(rlepixsuby$lengths)}
  # } else {
  #   curr_inds <- curr_inds[-length(curr_inds)]
  #   indlocs$tsw[curr_inds] <- rep(twosuby, rlengths)[1:length(curr_inds)]
  # }
}


#fill in remaining parameters of indlocs
twoy <- rlepixy$lengths * pxlw #total width on each row 'y'
indlocs$tw <- rep(twoy, rlepixy$lengths)
indlocs$li <- rep(seq_along(twoy), rlepixy$lengths)

#get colors for indlocs
pixcols <- t(sapply(1:nrow(indhits), function(i) img[indhits[i,1], indhits[i,2],]))
pixcols <- sapply(1:nrow(pixcols), function(i) rgb(pixcols[i,1], pixcols[i,2], pixcols[i,3], pixcols[i,4]))
indlocs$col <- pixcols

#see if we have the correct info
png(paste0("~/Pictures/", img_name, "_test-indlocs.png"), width = ncol(hits) * 5, height = nrow(hits) * 5, units = "px")
par(mar = c(0,0,0,0))
plot(NA, NA, xlim = c(0,1), ylim = c(0,1), frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(indlocs$x, indlocs$y, col = indlocs$col, pch = 15, cex = 1)
dev.off()

#### start preprocessing drawing ####
npixels <-c(ncol(hits) * 5, nrow(hits) * 5)
npixels <- c(12, 12) * 300
png(paste0("~/Pictures/", img_name, "_", txt_name, "_textfill_bcol-1.png"), width = npixels[2], height = npixels[1], units = "px")
par(mar = c(0,0,0,0))
plot(NA, NA, xlim = c(0,1), ylim = c(0,1), frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")

#specify or retrieve some graphical parameters
n_rounds <- 10
cex_history <- matrix(NA, nrow = n_rounds, ncol= 2)
flanking <- F
wiggle <- F
max_ind <- length(txt)
# max_ind <- length(strsplit(paste0(1:84, " ", collapse = ""), "")[[1]])-2
background_cols <- c("black","white")
render_image <- T
adjust_ws_kerning <- T
adjust_colors <- F

#initialize parameters
paintable_area <- sum(hits) / length(hits)
charh_multiplier <- 1.2
txta <- strwidth(paste0(txt, collapse = ""), units = "user", font = font, family = family) * 
  strheight(paste0(txt, collapse = ""), units = "user", font = font, family = family) * charh_multiplier
cex <- sqrt(paintable_area / txta)
prev_subline_i <- 1

#start iterative optimization of cex
for(round in 1:n_rounds){
  
  toolong <- F
  cat(paste0("\n", round, ": "))
  charwidths <- setNames(strwidth(uniqchars, units = "user", cex = cex, font = font, family = family), uniqchars)
  charheights <- setNames(strheight(uniqchars, units = "user", cex = cex, font = font, family = family) * charh_multiplier, uniqchars)
  charheight <- as.numeric(charheights[1])
  txtw <- cumsum(charwidths[txt])
  
  #figure out charater locations and colors
  ind_indlocs <- 1
  # ind_indlocs <- min(which(indlocs$tw > 0.5)) #if you want to start on a more encompassing line
  charlocx <- indlocs$x[ind_indlocs]
  charlocy <- indlocs$y[ind_indlocs]
  tw <- indlocs$tw[ind_indlocs]
  tsw <- indlocs$tsw[ind_indlocs]
  twtd <- 0
  tswtd <- 0
  charmat <- matrix(NA, ncol = 6, nrow = length(txt))
  charmat <- as.data.frame(charmat)
  colnames(charmat) <- c("x", "y", "col", "line", "block", "char")
  charmat$char <- txt
  line_i <- 1
  curr_block <- 1
  go_to_next_line <- F
  
  for(i in 1:max_ind){
    if(i %% round(length(txt) / 10) == 0) cat(i / round(length(txt) / 10) * 10, " ")
    
    charmat[i,1:5] <- c(charlocx, charlocy, NA, line_i, curr_block)
    charlocx <- charlocx + charwidths[txt[i]]
    
    twtd <- twtd + charwidths[txt[i]]
    tswtd <- tswtd + charwidths[txt[i]]
    # print(paste0(i, ": ", paste0(charmat[i,], collapse =" ")))
    if((tswtd > tsw) & !(twtd > tw)){
      twtd <- twtd + tsw - tswtd #to not falsely accumulate extra subwidths when they are small
      
      closest_indlocy <- indlocs$y[which.min(abs(charlocy - indlocs$y))]
      valid_indlocs <- indlocs[which(abs(closest_indlocy - indlocs$y) < 1E-9),]
      # valid_indlocs <- indlocs[which(abs(charlocy - indlocs$y) < 1E-9),]
      if(min(which(charlocx < valid_indlocs$x)) == Inf){
        go_to_next_line <- T
        # print("inf val"); break
      } else {
        block_inds <- which(charmat$block == curr_block & charmat$line == line_i)
        # new_ind_indlocs <- which(charmat[,5] == (curr_block+1) & charmat[,4] == line_i)
        new_ind_indlocs <- min(which((charlocx + 0.1*pxlw) < valid_indlocs$x)) #machine precision hack
        charlocx <- valid_indlocs$x[new_ind_indlocs]
        
        #spread whitespace so that you're on the edges of the current area
        if(round == n_rounds & adjust_ws_kerning){
          block_center <- charmat[block_inds[1],1] + tsw / 2
          
          if(length(block_inds) > 1){
            if(txt[max(block_inds)] == " "){
              # charmat[max(block_inds),1] <- charmat[max(block_inds)-1,1]
              block_inds <- block_inds[-length(block_inds)]
            }
            if(txt[min(block_inds)] == " "){
              # charmat[min(block_inds),1] <- charmat[min(block_inds)+1,1]
              block_inds <- block_inds[-1]
            }
            if(sum(txt[block_inds] != " ") == 1){
              # charmat[block_inds,1] <- charmat[block_inds,1] + (tsw - tswtd) / 2 #centers single character
              non_ws_ind <- block_inds[which(txt[block_inds] != " ")]
              charmat[non_ws_ind,1] <- block_center - charwidths[txt[non_ws_ind]] / 1.5 #to accommodate pos = 4
            } else {
              missing_space <- tsw - (tswtd)
              which_ws <- which(txt[block_inds] == " ")
              per_ws_incr <- missing_space / length(which_ws)
              per_ws_incr <- max(per_ws_incr, -charwidths[" "] / 2) #never remove more than half a space
              incr_vec <- rep(0, length(block_inds))
              incr_vec[which_ws] <- per_ws_incr
              incr_vec <- cumsum(incr_vec)
              
              #not centering characters properly -- find block center and char center and recenter chars in block
              char_center <- charmat[block_inds[1],1] + (diff(range(charmat[block_inds,1] + incr_vec)) + charwidths[charmat$char[block_inds[length(block_inds)]]]) / 2
              center_disp <- block_center - char_center
              
              #finally, apply the adjustment
              charmat[block_inds,1] <- charmat[block_inds,1] + incr_vec + center_disp
            }
          }
        }
        
        prev_subline_i <- i + 1
        curr_block <- curr_block + 1
        tswtd <- 0
        tsw <- valid_indlocs$tsw[new_ind_indlocs]
        
      }
    }
    if(twtd > tw | go_to_next_line){
      go_to_next_line <- F
      block_inds <- which(charmat$block == curr_block & charmat$line == line_i)
      
      #spread whitespace so that you're on the edges of the current area
      if(round == n_rounds & adjust_ws_kerning){
        block_center <- charmat[block_inds[1],1] + tsw / 2
        
        if(length(block_inds) > 1){
          if(txt[max(block_inds)] == " "){
            # charmat[max(block_inds),1] <- charmat[max(block_inds)-1,1]
            block_inds <- block_inds[-length(block_inds)]
          }
          if(txt[min(block_inds)] == " "){
            # charmat[min(block_inds),1] <- charmat[min(block_inds)+1,1]
            block_inds <- block_inds[-1]
          }
          if(sum(txt[block_inds] != " ") == 1){
            # charmat[block_inds,1] <- charmat[block_inds,1] + (tsw - tswtd) / 2 #centers single character
            non_ws_ind <- block_inds[which(txt[block_inds] != " ")]
            charmat[non_ws_ind,1] <- block_center - charwidths[txt[non_ws_ind]] / 1.5 #to accommodate pos = 4
          } else {
            missing_space <- tsw - (tswtd)
            which_ws <- which(txt[block_inds] == " ")
            per_ws_incr <- missing_space / length(which_ws)
            per_ws_incr <- max(per_ws_incr, -charwidths[" "] / 2) #never remove more than half a space
            incr_vec <- rep(0, length(block_inds))
            incr_vec[which_ws] <- per_ws_incr
            incr_vec <- cumsum(incr_vec)
            
            #not centering characters properly -- find block center and char center and recenter chars in block
            char_center <- charmat[block_inds[1],1] + (diff(range(charmat[block_inds,1] + incr_vec)) + charwidths[charmat$char[block_inds[length(block_inds)]]]) / 2
            center_disp <- block_center - char_center
            
            #finally, apply the adjustment
            charmat[block_inds,1] <- charmat[block_inds,1] + incr_vec + center_disp
          }
        }
      }
      
      
      twtd <- 0
      tswtd <- 0
      curr_block <- 1
      line_i <- line_i + 1
      charlocy <- charlocy - charheight
      ind_indlocs <- which.min(abs(charlocy - indlocs$y))
      tw <- indlocs$tw[ind_indlocs]
      tsw <- indlocs$tsw[ind_indlocs]
      charlocx <- indlocs$x[ind_indlocs]
      if(charlocy > min(indlocs$y)){
        # charlocy <- indlocs$y[ind_indlocs] #this is why I am getting final line palimpset  
      } else {
        toolong <- T
        got_through_prop <- i / length(txt)
        cat(paste0("ending at ", round(got_through_prop, 3) * 100, "% "))
        break
      }
    }
    
    #end of character loop
    # charmat[i,1:5] <- c(charlocx, charlocy, NA, line_i, curr_block)
    # charlocx <- charlocx + charwidths[txt[i]]
    
  }
  
  #find relation between last character and last pixel
  closest_ind <- which.min(apply(t(t(indlocs[,1:2]) - as.numeric(tail(charmat, 1)[1:2])), 1, pythag))
  prop_through <- closest_ind / nrow(indlocs)
  
  #adjust cex 
  newcexprop <- 0.95
  if(T){
    if(!toolong){
      newcex <- cex / prop_through^0.5 + cex / 400
      print(newcex / cex)
      cex_history[round,] <- c(cex, newcex / cex)
    } else {
      newcex <- cex * got_through_prop^0.5 - cex / 400
      print(newcex / cex)
      cex_history[round,] <- c(cex, newcex / cex)
    }
    flanking <- (!all(cex_history[1:round,2] > 1) & !all(cex_history[1:round,2] < 1)) & round != 1
    if(!flanking){
      cex <- newcex * newcexprop + cex * (1-newcexprop)
    } else {
      curr_cex_history <- cex_history[1:round,]
      cvs <- curr_cex_history[,1]
      rats <- curr_cex_history[,2]
      closest_below <- which(rats < 1) #too much text
      closest_below <- cvs[closest_below[which.max(rats[closest_below])]]
      closest_above <- which(rats > 1) #too little text
      closest_above <- cvs[closest_above[which.min(rats[closest_above])]]
      # cex_weights <- (1 / abs(1-cex_history[1:round,2]))^0.5
      # cex_weights <- cex_weights / sum(cex_weights)
      # cex <- sum(cex_history[1:round,1] * cex_weights) * 0.9 + propcex * 0.1
      if(round != n_rounds){
        cex <- (closest_above * 0.25 + closest_below * 0.75)
        # cex <- cex * ifelse(abs(rats[round] - rats[round-1]) < 1E-4 & wiggle, exp(rnorm(1, 0, 0.01)), 1)
      }
      if(round == (n_rounds-1)){   
          cex <- curr_cex_history[which.min(abs(curr_cex_history[,2]-1)),1]
          
          #get just the ones where all the text is there
          cvs <- curr_cex_history[,1]
          rats <- curr_cex_history[,2]
          closest_above <- which(rats > 1)
          cex <- cvs[closest_above[which.min(rats[closest_above])]]
      }
    }
  }

}

dev.off()


#figure out color info
nlines <- length(unique(charmat[,4]))
rlelines <- rle(charmat[,4])
line_inds <- cumsum(c(0,rlelines$lengths))
background_col_rgb <- t(col2rgb(background_cols[1]))
for(i in 1:nlines){
  cmatinds <- (line_inds[i]+1):(line_inds[i+1])
  subcm <- charmat[cmatinds,]
  yloc <- as.numeric(subcm[1,2])
  closest_indlocy <- indlocs$y[which.min(abs(yloc - indlocs$y))]
  subindlocs <- indlocs[abs(indlocs$y - closest_indlocy) < 1E-9,]
  xlocs <- subcm$x + charwidths[subcm$char] / 2
  
  nearest_matches <- findInterval(xlocs, c(-Inf, subindlocs$x[-1]-diff(subindlocs$x)/2, Inf))
  native_cols <- subindlocs$col[nearest_matches]
  
  if(adjust_colors){
    charprops <- propspace_uniqchar[charmat$char[cmatinds]]^0.5
    native_col_rgb <- t(col2rgb(native_cols))
    needed_cols <- data.frame(red = (native_col_rgb[,1] -  background_col_rgb[1] * (1-charprops)) / charprops,
                         blue = (native_col_rgb[,2] -  background_col_rgb[2] * (1-charprops)) / charprops,
                         green = (native_col_rgb[,3] -  background_col_rgb[3] * (1-charprops)) / charprops)
    needed_cols <- needed_cols / 255
    needed_cols <- needed_cols / apply(cbind(needed_cols, 1), 1, max)
    charmat[cmatinds,3] <- sapply(1:nrow(needed_cols), function(coli) ifelse(!any(is.na(needed_cols[coli,])), 
                                                      rgb(needed_cols$red[coli], needed_cols$blue[coli], needed_cols$green[coli], alpha = 1),
                                                      NA))
    
  } else {
    charmat[cmatinds,3] <- native_cols
  }
  
  
}

#### draw the characters ####

if(abs(rotation_angle) > 1E-3){
  rotmat <- rotmat_00(rotation_angle)
  new_pixcoords <- t(rotmat %*% t(cbind(charmat$x - 0.5, charmat$y - 0.5)))
  charmat$x <- new_pixcoords[,1] + 0.5
  charmat$y <- new_pixcoords[,2] + 0.5
}


for(bcoli in 1:length(background_cols)){
  if(render_image){
    png(paste0("~/Pictures/", img_name, "_", txt_name, "_textfill_bcol-", bcoli, ".png"), width = npixels[2], height = npixels[1], units = "px")
    par(mar = c(0,0,0,0))
    plot(NA, NA, xlim = c(0,1), ylim = c(0,1), frame.plot = F, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    rect(-10, -10, 10, 10, col = background_cols[bcoli], border = background_cols[bcoli], xpd = NA)
    
    for(i in (1:max_ind)){
      if(i %% round(length(txt) / 10) == 0) cat(i / round(length(txt) / 10) * 10, " ")
      text(labels = charmat$char[i], x = charmat$x[i] - ifelse(txt[i] == "Rethink Priorities", 0.02, 0), y = charmat$y[i],
           col = charmat$col[i], cex = cex, font = font, pos = 4, family = family, srt = rotation_angle)
      
    }
  }
  dev.off()
}

