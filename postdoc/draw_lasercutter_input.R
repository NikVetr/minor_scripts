library(sf)
library(igraph)
library(extrafont)
library(extrafontdb)

#functions
source("~/repos/polylines/R/functions.R")

#import fonts
# font_import()
# loadfonts()
names(pdfFonts())[grepl("oboto", names(pdfFonts()))]

#read in and process names
labnames <- trimws(readLines("~/montgomery_lab_names.txt", warn = F))
space_locs <- gregexpr(" ", labnames)
n_max_words <- 2
labnames_list <- lapply(1:length(labnames), function(ni){
  
  #process spaces
  sls <- c(space_locs[[ni]])
  if(n_max_words == 2){
    sls <- sls[1]
  } else if(n_max_words == 3){
    sls <- unique(c(sls[1], sls[length(sls)]))
  } else {
    sls <- sls[1:max(min(length(sls), n_max_words-1), 1)]
  }
  
  #divide up words
  starts <- c(1, sls + 1)
  stops <- c(sls - 1, nchar(labnames[ni]))
  labname_split <- sapply(1:length(starts), function(si) 
    substr(labnames[ni], start = starts[si], stop = stops[si]))
  labname_split
})

#specify dimensions of plotting (in inches)
dim_grid <- c(8,8)
dim_unit <- c(w = 2.19, h = 1.22)
dim_spacing <- c(w = 0.72, h = 0.19)
start_loc <- c(x = 1.7, y = 1.56)
dim_plot <- c(w = 23.77, h = 13.0)
dim_name <- dim_unit - dim_unit["w"] / 5
degree_name_overlap <- 0.3
name_cex_ratio_base <- 0.55
name_cex_ratio_factor <- 1.1
cex_scale <- 0.95
vws_scale <- 0.5
vws_adj_yqpgj <- 1.75

#construct matrix of centers
centers <- as.matrix(expand.grid(1:dim_grid[1], 1:dim_grid[2])) - 1
centers <- centers %*% diag(dim_unit + dim_spacing)
centers <- t(t(centers) + start_loc)
colnames(centers) <- c("x", "y")

#specify font information
font_name <- "Source Serif Pro"

#construct list of indiv name positions
name_metadata <- lapply(labnames_list, function(name){
  
  #extract important string information
  nw <- strwidth(name, "inches", family = font_name)
  nhs <- strheight(name, "inches", family = font_name)
  nh <- strheight(paste0(name, "AygT", collapse = "\n"), "inches")
  has_yqpgj <- any(strsplit(name[1], "")[[1]] %in% strsplit("yqpgj", "")[[1]])
  vws <- (nh - sum(nhs)) * vws_scale * ifelse(has_yqpgj, vws_adj_yqpgj, 1)
  nn <- length(name)
  
  if(nchar(name[2]) > 8){
    name_cex_ratio <- name_cex_ratio_base * name_cex_ratio_factor / (nchar(name[2]) / 8)
  } else {
    name_cex_ratio <- name_cex_ratio_base
  }
  
  #find the appropriate cex for each word
  new_cex_base <- min(dim_name["w"] / max(nw), dim_name["h"] / nh)
  new_cex <- rep(new_cex_base, nn) * c(2 / (name_cex_ratio + 1), 
                                       rep(2 * name_cex_ratio / (name_cex_ratio + 1) , nn-1))
  new_cex_rescale <- min(dim_name["w"] / max(nw * new_cex), 
                         dim_name["h"] / sum((nhs + vws * (nn - 1)) * new_cex))
  new_cex <- new_cex_rescale * new_cex * cex_scale
  
  #find new dimensions of words
  new_nw <- new_cex * nw
  new_nhs <- new_cex * nhs
  new_vws <- mean(new_cex * vws)
  new_nh <- sum(new_nhs) + new_vws * (nn - 1)
  
  #find displacements to these dimensions
  out <- data.frame(matrix(NA, ncol = 4, nrow = nn))
  colnames(out) <- c("text", "dx", "dy", "cex")
  out$text <- name
  out$cex <- new_cex
  
  #first for first word
  out$dx[1] <- new_nw[1] / 2 - dim_name["w"] / 2
  out$dy[1] <- new_nh / 2 - new_nhs[1] / 2
  
  #then for last word
  out$dx[nn] <- dim_name["w"] / 2 - new_nw[nn] / 2
  out$dy[nn] <- new_nhs[nn] / 2 - new_nh / 2
  
  #leave middle word alone? (can work out proportionate disp later)
  if(nn == 3){
    out$dx[2] <- 0
    out$dy[2] <- 0
  }
  
  #scooch names in if there is room
  out$w <- new_nw
  out$lb <- out$dx - out$w / 2
  out$rb <- out$dx + out$w / 2
  overlap <- out$rb[1] - out$lb[2]
  overlap_prop <- overlap / out$w
  
  if(max(overlap_prop) < degree_name_overlap){
    shift_in_by <- out$w[which.max(abs(overlap_prop))] * degree_name_overlap - overlap
    room_available <- dim_name["w"] / 2 - c(out$rb[1], -out$lb[2])
    out_of_room <- room_available < (shift_in_by / 2)
    shift_each_by <- rep(shift_in_by / 2, nn)
    shift_each_by[out_of_room] <- room_available[out_of_room]
    shift_each_by <- shift_each_by * c(1, -1)
    out$dx <- out$dx + shift_each_by
  }
  
  #recalculate boundaries
  out$w <- new_nw
  out$lb <- out$dx - out$w / 2
  out$rb <- out$dx + out$w / 2
  
  #record other useful info
  out$h <- new_nhs
  out$ws <- new_vws
  
  return(out)
  
})

#### plot ####

#now do the plotting
svg(filename = "~/lab_names_grid.svg", width = dim_plot["w"], height = dim_plot["h"])

#initiate plotting window
par(mar = c(0,0,0,0))
plot(1,1, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = F,
     xlim = c(0, dim_plot["w"]),
     ylim = c(0, dim_plot["h"])
)

for(i in 1:prod(dim_grid)){
  
  cat(paste0(i, " "))
  
  #plot rects for troubleshooting
  rect(xleft = centers[i,"x"] - dim_unit["w"]/2, 
       ybottom = centers[i,"y"] - dim_unit["h"]/2, 
       xright = centers[i,"x"] + dim_unit["w"]/2, 
       ytop = centers[i,"y"] + dim_unit["h"]/2, border = "grey20")
  # 
  # rect(xleft = centers[i,"x"] - dim_name["w"]/2, 
  #      ybottom = centers[i,"y"] - dim_name["h"]/2, 
  #      xright = centers[i,"x"] + dim_name["w"]/2, 
  #      ytop = centers[i,"y"] + dim_name["h"]/2, lty = 2, border = "grey20")
  # 

  #plot names
  if(i <= length(labnames_list)){
    
    #write name
    for(j in 1:length(labnames_list[[i]])){
      text(x = centers[i,1] + name_metadata[[i]]$dx[j], 
           y = centers[i,2] + name_metadata[[i]]$dy[j], 
           labels = name_metadata[[i]]$text[j], family = font_name, 
           cex = name_metadata[[i]]$cex[j], xpd = NA)
    }
    
    #draw DNA embellishment
    
    #need to calculate necessary height of double helix
    
    #then calculate number of twists to fill the space of that height
    #amplitude is 1 by default, so total thickness of DH is 2
    
    #then calculate length of line to underline the surname
    
    if(i %in% c(1:100)){
      nf_width <- name_metadata[[i]]$rb[2] - name_metadata[[i]]$lb[1] #nf = namefield
      prop_surname <- name_metadata[[i]]$w[2] * 
        (nchar(name_metadata[[i]]$text[2]) + 0.75) / nchar(name_metadata[[i]]$text[2]) / 
        nf_width
      helix_width <- (1-prop_surname) * nf_width * 1.1
      helix_height <- (name_metadata[[i]]$h[2] / 2 + name_metadata[[i]]$ws[2]) * 1.1
      if(helix_width > nf_width / 20){
        ntwists <- max(floor(helix_width / helix_height / (pi) * 2), 1)
        dh_bounds <- c(-pi/2, -pi/2 + pi * (ntwists))
        extend_straight_by <- diff(dh_bounds) * prop_surname / (1-prop_surname) + nf_width / 5
        box_dim <- c(nf_width,
                     name_metadata[[i]]$h[2] + name_metadata[[i]]$ws[2])
        realized_helix_height <- nf_width / (diff(dh_bounds) + extend_straight_by) * 2
        tweaked_amplitude <- helix_height / realized_helix_height
        dna_center <- centers[i,] + 
          c(0, 
            name_metadata[[i]]$dy[2] - 
            name_metadata[[i]]$h[2]/10)
        draw_DNA(target_center = dna_center, dh_bounds = dh_bounds, amplitude = tweaked_amplitude,
                 box_dim = box_dim, rot = 90, strand_thickness = 0.01, 
                 extend_straight = ifelse(ntwists %% 2 == 0, 2, 1), 
                 extend_straight_by = extend_straight_by)  
      }
      # if(i==60) break
    }
      
  }
  
  
  
    
}

dev.off()

