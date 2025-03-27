library(sf)
library(igraph)
library(extrafont)
library(extrafontdb)

#### functions ####

#from previous scripts
source("~/repos/polylines/R/functions.R")
source("~/scripts/minor_scripts/postdoc/true_text_dim_functions.R")
source("~/scripts/minor_scripts/postdoc/rect_constraint.R")

#for this specific script
transform_points <- function(x, orig_span, orig_range, relative_range) {
  new_range <- orig_span[1] + relative_range * diff(orig_span)
  mapped_x <- new_range[1] + (x - orig_range[1]) * (diff(new_range) / diff(orig_range))
  return(mapped_x)
}

find_gaps <- function(df, y, widen_by = 0.0) {
  x <- df[,2]
  below <- x < y
  changes <- diff(c(FALSE, below, FALSE))
  starts <- which(changes == 1)
  ends <- which(changes == -1) - 1
  bound_inds <- data.frame(start = starts, end = ends)
  # bounds <- data.frame(start = df[starts,1], end = df[ends,1])
  # bounds <- do.call(rbind, apply(bounds, 1, sort, simplify = F))
  bounds <- do.call(rbind, apply(bound_inds, 1, function(bis) 
    range(df[bis[1]:bis[2],1]), simplify = F))
  bounds[,1] <- bounds[,1] - widen_by
  bounds[,2] <- bounds[,2] + widen_by
  return(bounds)
}

poly2glyph <- function(concat_poly, orig_poly){
  glyph_inds <- rep(1:length(orig_poly), sapply(orig_poly, nrow))
  return(split(data.frame(concat_poly), glyph_inds))
}

fix_polys <- function(glyphs, glyphs_idx, close_glyphs = F){
  lapply(seq_along(glyphs), function(glyph_i){
    glyph <- glyphs[[glyph_i]]
    glyph_idx <- glyphs_idx[[glyph_i]]
    subglyphs <- lapply(split(glyph, glyph_idx), as.matrix) #split glyph on IDs
    subglyphs <- subglyphs[sapply(subglyphs, nrow) > 1] #clean up move idx
    
    #close these if not closed
    if(close_glyphs){
      subglyphs <- lapply(seq_along(subglyphs), function(gi){
        sg <- subglyphs[[gi]]
        if(!all(head(sg, 1) == tail(sg, 1))){
          sg <- rbind(sg, head(sg, 1))
        }
        return(sg)
      })  
    }
    
    # subglyphs <- subglyphs[c(2,3,1)] #order from outermost to innermost?
    glyph_poly <- sf::st_polygon(subglyphs)
    # glyph_poly <- ?sf::st_make_valid(glyph_poly)
    glyph_sf <- sf::st_sfc(glyph_poly)
    return(glyph_sf)
  })
}

dr <- function(x){diff(range(x))}

drrat <- function(x){dr(x[,1]) / dr(x[,2])}

find_closest_distance <- function(glyph1, glyph2) {
  min_distance <- sf::st_distance(glyph1, glyph2, by_element = FALSE)
  return(min(min_distance))
}

estimate_stroke_width <- function(glyph) {
  A <- as.numeric(sf::st_area(glyph))  # Get polygon area
  P <- as.numeric(sf::st_length(sf::st_boundary(glyph)))  # Get polygon perimeter
  
  #use a rectangular approximation
  #so if I imagine a tall, skinny rectangle, the width is the area / length, 
  #because the area is the length * the width, A = W * L. 
  #The perimeter is 2*width + 2*length, or P = 2*W + 2*L. 
  #So to get the length I need to do L = P/2 - W). 
  #I can rearrange to get eg L = P/2 - A/L, and substitute in to get A = WP/2-W^2. 
  #Rearranging, we have W^2 - WP/2 + A = 0, 
  #which is a quadratic equation with solution
  
  a <- 1
  b <- - (P / 2)
  c <- A
  discriminant <- b^2 - 4 * a * c
  LW_values <- (-b + c(-1,1) * sqrt(discriminant)) / (2*a)
  
  return(min(LW_values))
}


convert_to_plot <- function(svgdf){
  elem_idx <- svgdf$elem_idx
  coords <- data.frame(x = svgdf$x, y = svgdf$y)
  coords$y <- max(coords$y) - coords$y
  split_coords <- split(coords, elem_idx)
  idx <- split(svgdf$idx, elem_idx)
  
  #check for weird artefacting
  #sometimes there is a point and it is a separate shape but all the same pt
  keep <- sapply(split_coords, nrow) - sapply(sapply(split_coords, duplicated), sum) != 1
  split_coords <- split_coords[keep]
  idx <- idx[keep]
  
  #convert to better behaved polys
  fixed_polys <- split_coords
  for(pli in 1:length(split_coords)){
    subpolys <- split(split_coords[[pli]], idx[[pli]])
    subpolys <- lapply(subpolys, as.matrix)
    subpolys <- lapply(seq_along(subpolys), function(gi){
      sg <- subpolys[[gi]]
      if(!all(head(sg, 1) == tail(sg, 1))){
        sg <- rbind(sg, head(sg, 1))
      }
      return(sg)
    })  
    fixed_polys[[pli]] <- sf::st_make_valid(sf::st_polygon(subpolys))
  }
  
  return(list(coords = coords, 
              fixed_polys = fixed_polys))
  
}

plot_polys <- function(fixed_polys, col = 1, border = NA, add = T, disp = c(0,0), ...){
  for(pli in 1:length(fixed_polys)){
    plot(fixed_polys[[pli]] + disp, col = col, border = border, add = add)
  }  
}

#### import data ####

#import fonts
# font_import()
# loadfonts()
names(pdfFonts())[grepl("oboto", names(pdfFonts()))]

#render text with polygons, not svg tags
showtext::showtext_auto(T)

#read in and process names
labnames <- trimws(readLines("~/montgomery_lab_names.txt", warn = F))
nlabn <- length(labnames)
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


#### parameters ####
#specify dimensions of plotting (in inches)
dim_grid <- c(8,8)
dim_unit <- c(w = 2.19, h = 1.22)
dim_spacing <- c(w = 0.705, h = 0.1926)
start_loc <- c(x = 1.7, y = 1.56)
name_center_adj_prop <- c(x = 0.01, y = 0.01)
name_center_adj <- name_center_adj_prop * dim_unit
logo_center_adj_prop <- c(x = 0, y = -0.005)
logo_center_adj <- logo_center_adj_prop * dim_unit
dim_plot <- c(w = 23.77, h = 13.0)
dim_name <- dim_unit - dim_unit["w"] / 5
dim_logo <- dim_unit - dim_unit["w"] / 10
degree_name_overlap <- 0.4
name_cex_ratio_base <- 0.55
name_cex_ratio_factor <- 1.1
cex_scale <- 0.95
vws_scale <- 0.6
vws_adj_yqpgj <- 1.5
twist_ratio_xy <- 3
prop_dist <- 0.2
use_smallcaps <- T
use_bolded_glyphs_for_smallcaps <- F
raise_J <- T
lower_Q <- T
thin_smallcaps_firstletter <- T
disp_y_CoM <- F

#construct matrix of centers
centers <- as.matrix(expand.grid(1:dim_grid[1], 1:dim_grid[2])) - 1
centers <- centers %*% diag(dim_unit + dim_spacing)
centers <- t(t(centers) + start_loc)
colnames(centers) <- c("x", "y")

#specify font information
font_name <- "Source Serif Pro"

#change first names for nicer layout if desired
if(use_smallcaps){
  first_names <- do.call(rbind, labnames_list)[,1]
  last_names <- do.call(rbind, labnames_list)[,2]
  new_first_names <- sapply(first_names, function(name){
    has_yqpgj <- any(strsplit(name, "")[[1]] %in% strsplit("yqpgj", "")[[1]])  
    if(has_yqpgj){
      return(toupper(name))
    } else {
      return(name)
    }
  })
  new_last_names <- sapply(seq_along(first_names), function(ni){
    name <- first_names[ni]
    surname <- last_names[ni]
    has_yqpgj <- any(strsplit(name, "")[[1]] %in% strsplit("yqpgj", "")[[1]])  
    if(has_yqpgj){
      return(toupper(surname))
    } else {
      return(surname)
    }
  })
  labnames_list <- apply(cbind(new_first_names, new_last_names), 1, function(x) as.character(x), simplify = F)
}

#construct indices for when we have more names than cells
nplots <- ceiling(nlabn / 64)
plot_inds <- 1:nlabn %% 64
plot_inds[plot_inds == 0] <- 64
plot_inds <- split(plot_inds, rep(1:nplots, each = 64)[1:nlabn])
plot_locs <- lapply(plot_inds, function(x){
  y <- x
  if(length(x) > 1){
    y[2] <- 64
    y[length(x)] <- 2
  }
  return(y)
})
plot_name_inds <- split((1:nlabn), rep(1:nplots, each = 64)[1:nlabn])
unlisted_pli <- unlist(plot_inds)

#### initial bpolys ####
#retrieve true location / bounding info for names
labnames_list_hash <- digest::digest(labnames_list)
name_info_path <- paste0("~/montgomery-lab-names_", 
                         gsub(" ", "-", font_name), "_", 
                         labnames_list_hash,".RData")
if(file.exists(name_info_path)){
  load(name_info_path)
} else {
  name_info <- lapply(labnames_list, function(name){
    cat(paste0(paste0(name, collapse = " "), ", "))
    lapply(name, function(subname){
      text_distortion(string = subname, font_name = font_name, device = "svg")
    })
  })
  bolded_name_info <- lapply(labnames_list, function(name){
    cat(paste0(paste0(name, collapse = " "), ", "))
    lapply(name, function(subname){
      if(subname == toupper(subname)){
        return(true_text_params(string = subname, 
                                target_center = c(mean(par("usr")[1:2]), mean(par("usr")[3:4])), 
                                cex = 1, device = "svg", font_name = font_name, 
                                return_td_res = T, font = 2))
      } else {
        return(NA)
      }
    })
  })
  names(name_info) <- names(bolded_name_info) <- sapply(labnames_list, paste0, collapse = " ")
  save(name_info, bolded_name_info, file = name_info_path)  
}


#get the bpolys for the bolded names too, for kerning purposes


#### initial positions ####
#construct list of indiv name positions
name_metadata <- lapply(labnames_list, function(name){
  
  #extract important string information
  nw <- strwidth(name, "inches", family = font_name)
  nhs <- strheight(name, "inches", family = font_name)
  nh <- strheight(paste0(name, "AygT", collapse = "\n"), "inches")
  has_yqpgj <- any(strsplit(name[1], "")[[1]] %in% strsplit("yqpgj", "")[[1]])
  vws <- (nh - sum(nhs)) * vws_scale * ifelse(has_yqpgj, vws_adj_yqpgj, 1)
  hws <- strwidth(" ", "inches", family = font_name)
  
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
  new_hws <- mean(new_cex * hws)
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
  out$vws <- new_vws
  out$hws <- new_hws
  
  #now let's try to get the DNA metadata in there too?
  
  return(out)
  
})

#### correct metadata ####
#let's now reprocess these metadata in the context of their true metadata
if(!("bpoly_info" %in% names(name_metadata[[1]]))){
  
  temp_filename <- paste0(tempfile(), ".svg")
  svg(filename = temp_filename, width = dim_plot["w"], height = dim_plot["h"])
  
  par(mar = c(0,0,0,0))
  plot.new()
  plot.window(
    xlim = c(0, dim_plot["w"]),
    ylim = c(0, dim_plot["h"])
  )
  
  if(use_smallcaps){
    
    #retrieve bolded and regular glyphs
    all_glyphs <- sort(unique(strsplit(toupper(paste0(unlist(labnames_list), collapse = "")), "")[[1]]))
    all_glyphs <- setdiff(all_glyphs, " ")
    all_glyphs_cat <- paste0(all_glyphs, collapse = "")
    bold_bpoly <- true_text_params(string = all_glyphs_cat, 
                                   target_center = c(mean(par("usr")[1:2]), mean(par("usr")[3:4])), 
                                   cex = 1, device = "svg", font_name = font_name, 
                                   return_td_res = T, font = 2)
    reg_bpoly <- true_text_params(string = all_glyphs_cat,
                                  target_center = c(mean(par("usr")[1:2]), mean(par("usr")[3:4])),
                                  cex = 1, device = "svg", font_name = font_name,
                                  return_td_res = T, font = 1)
    smallcaps_bpolys <- setNames(bold_bpoly$bpoly, all_glyphs)
    smallcaps_svgdf <- setNames(split(bold_bpoly$td_res$svgdf, bold_bpoly$td_res$svgdf$elem_idx)[-1], all_glyphs)
    
    #compute scale factor (to keep stroke width constant-ish)
    reg_areas <- sapply(fix_polys(reg_bpoly$bpoly, reg_bpoly$td_res$svgdf$idx[-c(1:4)], T), st_area)
    bold_areas <- sapply(fix_polys(bold_bpoly$bpoly, bold_bpoly$td_res$svgdf$idx[-c(1:4)], T), st_area)
    area_ratios <- reg_areas / bold_areas
    
    #bounding box ratios look very similar!
    reg_h <- sapply(reg_bpoly$bpoly, function(glyph) dr(glyph[,2]))
    reg_w <- sapply(reg_bpoly$bpoly, function(glyph) dr(glyph[,1]))
    
    bold_h <- sapply(bold_bpoly$bpoly, function(glyph) dr(glyph[,2]))
    bold_w <- sapply(bold_bpoly$bpoly, function(glyph) dr(glyph[,1]))
    
    h_ratios <- reg_h / bold_h
    w_ratios <- reg_w / bold_w
    
    bbox_area_ratios <- h_ratios * w_ratios
    
    #so most of the difference is probably due to stroke thickness
    #we can scale by that value to keep stroke width constant
    scale_by <- median(area_ratios)
    # plot(do.call(rbind, bold_bpoly$bpoly), type = "l", asp = 1)
    # plot(do.call(rbind, reg_bpoly$bpoly), type = "l", asp = 1)
  }
  
  dna_info <- draw_DNA(dh_bounds = pi/2*c(-1,1), amplitude = 1,
                       rot = 90, strand_thickness = 0.01, 
                       extend_straight = 1, 
                       extend_straight_by = 0, col = 1, 
                       topleft_xy = c(0,0), box_DNA = F, 
                       return_info = T,
                       proportional_distance = prop_dist)
  
  for(i in 1:nlabn){
    
    cat(paste0(i, " "))
    
    nm <- name_metadata[[i]]
    for(j in 1:length(labnames_list[[i]])){
      bpoly_info <- true_text_params(string = nm$text[j], 
                                     target_center = c(centers[unlisted_pli[i],1] + nm$dx[j], 
                                                       centers[unlisted_pli[i],2] + nm$dy[j]), 
                                     cex = nm$cex[j], device = "svg", font_name = font_name, 
                                     return_td_res = T)
      
      #raise capital J up to median height of other letters
      if(raise_J){
        if(grepl("J", name_metadata[[i]]$text[j])){
          name_chars <- strsplit(name_metadata[[i]]$text[j], "")[[1]]
          J_ind <- which(name_chars == "J")
          baseline_ys <- sapply(setdiff(1:length(name_chars), J_ind), function(nji){
            min(bpoly_info$bpoly[[nji]][,2])})
          baseline_y <- median(baseline_ys)
          for(ji in J_ind){
            bpoly_info$bpoly[[ji]][,2] <- bpoly_info$bpoly[[ji]][,2] - 
              min(bpoly_info$bpoly[[ji]][,2]) + baseline_y
          }
          bpoly_info$bbox <- bpoly2bbox(do.call(rbind, bpoly_info$bpoly))
        }
      }
      
      #if making the first name into smallcaps, adjust accordingly
      if(use_smallcaps){
        
        if(name_metadata[[i]]$text[1] == toupper(name_metadata[[i]]$text[1])){
          
          #retrieve glyph information
          name_chars <- strsplit(name_metadata[[i]]$text[j], "")[[1]]
          sc_bpoly <- bpoly_info$bpoly
          idx <- bpoly_info$td_res$svgdf$idx[-c(1:4)]
          split_idx <- split(idx, rep(1:length(sc_bpoly), sapply(sc_bpoly, nrow)))
          
          if(thin_smallcaps_firstletter){
            first_glyph <- fix_polys(list(sc_bpoly[[1]]), split_idx[[1]], close_glyphs = T)[[1]]
            fg_asp <- apply(bpoly2bbox(sc_bpoly[[1]]), 2, dr)
            
            stroke_width <- estimate_stroke_width(first_glyph)
            prop_sw_to_remove <- 1 - (1 * 0.5 + scale_by * 0.5)
            #adjust for top and bottom also losing a bit of space?
            lost_height <- prop_sw_to_remove * stroke_width * 2 / fg_asp[2]
            lost_width <- prop_sw_to_remove * stroke_width * 2 / fg_asp[1]
            #or not, for now!
            thinned_first_glyph <- sf::st_buffer(first_glyph, dist = -prop_sw_to_remove*stroke_width, nQuadSegs = 1)
            #this does round corners, but maybe that is fine for now
            # thinned_first_glyph <- st_simplify(thinned_first_glyph, dTolerance = 0.00001)  # Adjust tolerance as needed
            fg_coords <- sf::st_coordinates(thinned_first_glyph)  # Get x, y, and indices
            fg_poly <- cbind(x = fg_coords[,"X"], y = fg_coords[,"Y"])
            
            #transform the polygon to the original bounding box
            fg_poly[,1] <- (fg_poly[,1] - min(fg_poly[,1])) / dr(fg_poly[,1]) * fg_asp[1] + min(sc_bpoly[[1]][,1])
            fg_poly[,2] <- (fg_poly[,2] - min(fg_poly[,2])) / dr(fg_poly[,2]) * fg_asp[2] + min(sc_bpoly[[1]][,2])
            
            # plot(first_glyph, asp = 1, col = 2); axis(1, line  = -3); axis(2, line  = -3)
            # plot(thinned_first_glyph, asp = 1, add = T)
            # lines(fg_poly[,1], fg_poly[,2], col = 3, lwd = 3)
            
            #update the polygon
            sc_bpoly[[1]] <- as.data.frame(rbind(fg_poly, tail(fg_poly, 1)))
            
            #update the svgdf
            svgdf_split <- split(bpoly_info$td_res$svgdf, bpoly_info$td_res$svgdf$elem_idx)
            fg_svgdf <- svgdf_split[[2]]
            if(nrow(fg_svgdf) < nrow(fg_poly)){
              buffer_svgdf <- as.data.frame(matrix(NA, nrow(fg_poly) - nrow(fg_svgdf), ncol(fg_svgdf)))
              colnames(buffer_svgdf) <- colnames(fg_svgdf)
              fg_svgdf <- rbind(fg_svgdf, buffer_svgdf)
              fg_svgdf$elem_idx <- fg_svgdf$elem_idx[1]
            } else {
              fg_svgdf <- fg_svgdf[1:nrow(fg_poly),]
            }
            fg_svgdf$x <- fg_poly[,1]
            fg_svgdf$y <- fg_poly[,2]
            fg_svgdf$idx <- fg_coords[,"L2"]
            svgdf_split[[2]] <- rbind(fg_svgdf, tail(fg_svgdf, 1))
            svgdf_split[[2]]$idx[length(svgdf_split[[2]]$idx)] <- svgdf_split[[2]]$idx[length(svgdf_split[[2]]$idx-1)] + 1
            bpoly_info$td_res$svgdf <- do.call(rbind, svgdf_split)
            idx <- bpoly_info$td_res$svgdf$idx[-c(1:4)]
            split_idx <- split(idx, rep(1:length(sc_bpoly), sapply(sc_bpoly, nrow)))
            
          }
          
          #retrieve useful parameters about names
          orig_left_b <- sapply(sc_bpoly, function(x) min(x[,1]))
          orig_right_b <- sapply(sc_bpoly, function(x) max(x[,1]))
          orig_mid <- (orig_left_b + orig_right_b) / 2
          orig_glyph_w <- orig_right_b - orig_left_b
          orig_glyph_kern <- c(0, orig_left_b[-1] - orig_right_b[-length(orig_right_b)])
          
          orig_bot_b <- sapply(sc_bpoly, function(x) min(x[,2]))
          orig_top_b <- sapply(sc_bpoly, function(x) max(x[,2]))
          orig_glyph_h <- orig_top_b - orig_bot_b
          
          #sub in bolded glyphs for unbolded glyphs
          bolded_poly <- bolded_name_info[[i]][[j]]$bpoly
          bolded_idx <- bolded_name_info[[i]][[j]]$td_res$svgdf$idx[-c(1:4)]
          bolded_split_idx <- split(bolded_idx, rep(1:length(bolded_poly), sapply(bolded_poly, nrow)))
          bolded_left_b <- sapply(bolded_poly, function(x) min(x[,1]))
          bolded_right_b <- sapply(bolded_poly, function(x) max(x[,1])) 
          if(use_bolded_glyphs_for_smallcaps){
            sc_bpoly[2:length(name_chars)] <- bolded_poly[2:length(name_chars)]
            
            #to rescale glyphs to correct height, get average correction
            bolded_resc_correction <- mean(sapply(2:length(sc_bpoly), function(gi){
              glyph <- sc_bpoly[[gi]]
              unbolded_h <- orig_glyph_h[gi]
              bolded_h <- dr(glyph[,2])
              resc_correction <- unbolded_h / bolded_h
            }))
            glyph <- do.call(rbind, bolded_poly[-1])
            glyph_inds <- rep(1:length(bolded_poly[-1]), sapply(bolded_poly[-1], nrow))
          } else {
            resc_correction <- 1
            glyph <- do.call(rbind, sc_bpoly[-1])
            glyph_inds <- rep(1:length(sc_bpoly[-1]), sapply(sc_bpoly[-1], nrow))
          }
          
          #scale and shift accordingly
          scale_glyph_by <- scale_by * resc_correction
          
          #scale and translate glyph
          glyph[,1] <- (glyph[,1] - min(glyph[,1])) * scale_glyph_by #scale by factor
          glyph[,1] <- glyph[,1] + orig_right_b[1] + (orig_left_b[2] - orig_right_b[1]) * scale_by #return to previous location
          glyph[,2] <- (glyph[,2] - min(glyph[,2])) * scale_glyph_by #scale by factor
          glyph[,2] <- glyph[,2] + orig_bot_b[1] + (orig_bot_b[2] - orig_bot_b[1]) * scale_by  #return to previous location
          
          #split back into individual letters
          sc_bpoly[2:length(sc_bpoly)] <- split(glyph, glyph_inds)
          
          #to figure out kerning, use the mean distance between all the properly kerned glyphs
          if(use_bolded_glyphs_for_smallcaps){
            curr_fixed_glyphs <- fix_polys(glyphs = sc_bpoly[-1], unlist(bolded_split_idx[-1]), close_glyphs = T)
          } else {
            curr_fixed_glyphs <- fix_polys(glyphs = sc_bpoly[-1], unlist(split_idx[-1]), close_glyphs = T)
          }
          curr_closest_dists <- sapply(2:length(curr_fixed_glyphs), function(ci){
            find_closest_distance(glyph1 = curr_fixed_glyphs[[ci-1]], glyph2 = curr_fixed_glyphs[[ci]])
          })
          mean_closest_dist <- mean(curr_closest_dists) * mean(c(1, 1/scale_by)) #set this mean distance intermediate to the scale factor
          
          #this is a 1D, nonlinear optimization problem.
          shift_glyph <- function(glyph, shift_x) {
            glyph_shifted <- glyph
            glyph_shifted <- glyph_shifted + c(shift_x, 0)  # Shift x-coordinates
            return(glyph_shifted[[1]])
          }
          
          distance_sq_err <- function(shift_x, glyph1, glyph2, target_dist) {
            glyph2_shifted <- shift_glyph(glyph2, shift_x)
            current_dist <- find_closest_distance(glyph1, glyph2_shifted)
            # print(c(curr_dist = current_dist, shift_x = shift_x))
            if (current_dist <= 1e-5) {
              # If overlapping, assign a large penalty so the optimizer avoids it
              return(1e6)
            }
            out <- (current_dist - target_dist)^2
            return(out)
          }
          
          max_left_shift <- (min(sc_bpoly[[1]][,1]) - min(sc_bpoly[[2]][,1]))
          find_optimal_shift_optimize <- function(glyph1, glyph2, mean_closest_dist, search_range=c(max_left_shift, 10)) {
            result <- optimize(
              f = distance_sq_err, 
              interval = search_range, 
              glyph1 = glyph1, 
              glyph2 = glyph2, 
              target_dist = mean_closest_dist
            )
            return(result$minimum)
          }
          
          if(use_bolded_glyphs_for_smallcaps){
            first_fixed_glyph <- fix_polys(glyphs = list(sc_bpoly[[1]]), bolded_split_idx[[1]], close_glyphs = T)[[1]]  
          } else {
            first_fixed_glyph <- fix_polys(glyphs = list(sc_bpoly[[1]]), glyphs_idx = list(split_idx[[1]]), close_glyphs = T)[[1]]  
          }
          
          second_fixed_glyph <- curr_fixed_glyphs[[1]]
          optimal_shift <- find_optimal_shift_optimize(first_fixed_glyph, second_fixed_glyph, mean_closest_dist)
          
          #shift all the glyphs accordingly
          shift_x_by <- c(0, rep(optimal_shift, length(sc_bpoly)-1))
          sc_bpoly <- lapply(seq_along(sc_bpoly), function(gi){
            shifted_glyph <- sc_bpoly[[gi]]
            shifted_glyph[,1] <- shifted_glyph[,1] + shift_x_by[gi]
            return(shifted_glyph)
          })
          
          # plot(do.call(rbind, sc_bpoly), type = "l", asp = 1)
          # lines(do.call(rbind, bpoly_info$bpoly)[,1],
          #       do.call(rbind, bpoly_info$bpoly)[,2], type = "l", asp = 1, col = 2)
          # plot(do.call(rbind, bpoly_info$bpoly), type = "l", asp = 1)
          bpoly_info$bpoly <- sc_bpoly
          bpoly_info$bbox <- bpoly2bbox(do.call(rbind, sc_bpoly))
          
          #update the td_res svgdf
          new_svgdf_f2l <- split(bpoly_info$td_res$svgdf, bpoly_info$td_res$svgdf$elem_idx)[1:2]
          if(use_bolded_glyphs_for_smallcaps){
            new_svgdf_ll <- split(bolded_name_info[[i]][[j]]$td_res$svgdf, bolded_name_info[[i]][[j]]$td_res$svgdf$elem_idx)[-c(1:2)]  
          } else {
            new_svgdf_ll <- split(bpoly_info$td_res$svgdf, bpoly_info$td_res$svgdf$elem_idx)[-c(1:2)]  
          }
          
          new_svgdf <- c(new_svgdf_f2l, new_svgdf_ll)
          new_svgdf <- do.call(rbind, new_svgdf)
          bpoly_info$td_res$svgdf <- new_svgdf
        }
      }
      
      #record to variable
      name_metadata[[i]]$bpoly_info[j] <- list(bpoly_info)
      
    }
    
  }
  
  dev.off()
  
}

#### optimize placement ####
#using these as initial values, let's optimize the placement of graphical elements
plot_stuff <- T
b2_h <- max(sapply(dna_info$bp, function(bp) diff(range(bp$y)))) * (1-prop_dist)
b2_w <- diff(range(dna_info$s1[[1]]$x))
mean_iar2f <- b2_w / b2_h # inverse aspect ratio for b2 (w / h), corresponding to the DNA strand top nucs
cell_metadata <- lapply(1:nlabn, function(i) NULL)
for(i in 1:nlabn){
  
  cat(paste0(i, " "))
  
  #retrieve name info
  nm <- name_metadata[[i]]
  
  #get abstracted relative dimensions in visual units (inches)
  total_width <- 1
  total_height <- as.numeric(dim_name[2] / dim_name[1])
  b1_bbox <- nm$bpoly_info[[1]]$bbox
  ar1 <- diff(range(b1_bbox[,2])) / diff(range(b1_bbox[,1])) #h / w
  b3_bbox <- nm$bpoly_info[[2]]$bbox
  max_letter_b3_bbox_ind <- which.max(sapply(seq_along(nm$bpoly_info[[2]]$bpoly), function(li){
    diff(range(bpoly2bbox(nm$bpoly_info[[2]]$bpoly[[li]])[,2]))}))
  max_letter_b3_bbox <- bpoly2bbox(nm$bpoly_info[[2]]$bpoly[[max_letter_b3_bbox_ind]])
  max_letter_b3_bbox$Var1 <- nm$bpoly_info[[2]]$bbox$Var1
  ar3 <- diff(range(max_letter_b3_bbox[,2])) / diff(range(max_letter_b3_bbox[,1])) #h / w
  
  #since we scale from the tallest letter of the surname, we need to correct for this when warping later
  b3_h_rat_corr <- dr(nm$bpoly_info[[2]]$bbox[,2]) / dr(max_letter_b3_bbox[,2])
  
  #optimize placement of elements in fully generic sense
  rectangles <- place_rectangles_bfgs(
    prop_dist   = prop_dist,
    total_width = total_width,
    total_height = total_height,
    ar1         = ar1,
    ar3         = ar3,
    start_w1    = 0.5,
    start_w2    = 0.4,
    start_h2    = 0.4 / mean_iar2f,
    big_penalty = 1e3,
    ratio_mean_b1b2_w  = 2.5,
    ratio_sd_b1b2_w    = 10,
    ratio_mean_b1b2_h  = 2.5,
    ratio_sd_b1b2_h    = 0.5,
    ratio_mean_b2b3_w = 1,
    ratio_sd_b2b3_w = 10,
    ratio_mean_b1_b2b3_w = 0.925,
    ratio_sd_b1_b2b3_w = 0.05,
    mean_iar2f = mean_iar2f,
    sd_iar2f = 0.2,
    mean_ntwist = 1.5,
    sd_ntwist = 0.35,
    iar2f_weight = 1,
    two_stage_iar2f = T, 
    nrep = 20,
    print_penalties = F
  )
  
  rect_info <- get_rect_info(rectangles)
  ntwists <- round(rect_info$ar2_factor)
  if(plot_stuff){
    plot_rectangles(rectangles)
  }
  
  #hmm ar2 factor not being processed properly
  #rectangles report the same factor, but box_dim[1] / box_dim[2] differs a lot!
  
  #plot locations
  if(plot_stuff){
    par(mar = c(0,0,0,0))
    plot(1,1, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = F,
         xlim = c(0, dim_unit["w"]),
         ylim = c(0, dim_unit["h"]),
         asp = 1
    )
  }
  
  #outer rect
  outer_bbox <- bpoly2bbox(rbind(c(dim_unit["w"]/2 - dim_unit["w"]/2, dim_unit["h"]/2 - dim_unit["h"]/2),
                                 c(dim_unit["w"]/2 + dim_unit["w"]/2, dim_unit["h"]/2 + dim_unit["h"]/2)))
  cell_metadata[[i]]$outer_bbox <- outer_bbox
  if(plot_stuff){
    lines(outer_bbox, col = "grey20")
  }
  
  #inner rect
  inner_bbox <- bpoly2bbox(rbind(c(dim_unit["w"]/2 - dim_name["w"]/2, dim_unit["h"]/2 - dim_name["h"]/2),
                                 c(dim_unit["w"]/2 + dim_name["w"]/2, dim_unit["h"]/2 + dim_name["h"]/2)))
  cell_metadata[[i]]$inner_bbox <- inner_bbox
  if(plot_stuff){
    lines(inner_bbox, col = "grey20", lty = 2)
  }
  
  #plot names
  resc_rect <- T
  
  #write name
  nm_polys <- list()
  nm_polys_glyphs <- list(NULL, NULL)
  nm_polys_bbox <- list()
  for(j in 1:length(labnames_list[[i]])){
    
    nm_poly <- do.call(rbind, nm$bpoly_info[[j]]$bpoly)
    nm_poly_plot <- cbind(nm_poly[,1] - centers[unlisted_pli[i],"x"] + dim_unit["w"]/2, 
                          nm_poly[,2] - centers[unlisted_pli[i],"y"] + dim_unit["h"]/2)
    
    if(resc_rect){
      resc_poly <- nm_poly_plot
      x_buffer <- (1 - rect_info$total_width) / 2
      y_buffer <- (total_height - rect_info$total_height) / 2 / total_height
      
      #to rescale first name
      if(j==1){
        x_r <- c(rectangles$b1$x, rectangles$b1$x + rectangles$b1$width)
        x_rr <- x_r + x_buffer
        resc_poly[,1] <- transform_points(x = resc_poly[,1],
                                          orig_span = range(inner_bbox[,1]),
                                          orig_range = range(resc_poly[,1]),
                                          relative_range = x_rr)
        
        y_r <- (c(rectangles$b1$y - rectangles$b1$height, rectangles$b1$y) + total_height) / total_height
        y_rr <- y_r - y_buffer
        resc_poly[,2] <- transform_points(x = resc_poly[,2],
                                          orig_span = range(inner_bbox[,2]),
                                          orig_range = range(resc_poly[,2]),
                                          relative_range = y_rr)
      } else {
        #to rescale second name
        x_r <- c(rectangles$b3$x, rectangles$b3$x + rectangles$b3$width)
        x_rr <- x_r + x_buffer
        resc_poly[,1] <- transform_points(x = resc_poly[,1],
                                          orig_span = range(inner_bbox[,1]),
                                          orig_range = range(resc_poly[,1]),
                                          relative_range = x_rr
        )
        
        y_rr <- c(
          ((total_height - rect_info$total_height) / 2 + rectangles$h2 * prop_dist / 2) / total_height,
          ((total_height - rect_info$total_height) / 2 + rectangles$h3 + rectangles$h2 * prop_dist / 2) / total_height
        )
        y_r <- (c(-rectangles$b3$y - rectangles$b3$height, -rectangles$b3$y) + total_height) / total_height
        y_rr <- y_r - y_buffer
        y_rr[1] <- y_rr[2] - diff(y_rr) * b3_h_rat_corr
        resc_poly[,2] <- transform_points(x = resc_poly[,2],
                                          orig_span = range(inner_bbox[,2]),
                                          orig_range = range(resc_poly[,2]),
                                          relative_range = y_rr)
      }
      
      #record the overall polygon
      nm_poly_plot <- nm_polys[[j]] <- resc_poly
      
      #split the polygon into individual glyphs (using the original poly indices)
      nm_polys_glyphs[[j]]$orig_poly <- poly2glyph(concat_poly = resc_poly, 
                                                   orig_poly = nm$bpoly_info[[j]]$bpoly)
      
      #convert these to valid polygons for plotting
      glyphs_idx <- poly2glyph(concat_poly = nm$bpoly_info[[j]]$td_res$svgdf$idx[-c(1:4)], 
                               orig_poly = nm$bpoly_info[[j]]$bpoly)
      nm_polys_glyphs[[j]]$fixed_poly <- fix_polys(glyphs = nm_polys_glyphs[[j]]$orig_poly, 
                                                   glyphs_idx = glyphs_idx)
    }
    
    #plot name
    if(plot_stuff){
      # lines(nm_poly_plot)
      for(gi in 1:length(nm_polys_glyphs[[j]]$fixed_poly)){
        plot(nm_polys_glyphs[[j]]$fixed_poly[[gi]], col = 1, border = NA, add = T)
      }
      
    }
    nm_poly_plot_bbox <- bpoly2bbox(nm_poly_plot)
    if(plot_stuff){
      lines(nm_poly_plot_bbox, lty = 2, col = 2)
    }
    nm_polys_bbox[[j]] <- nm_poly_plot_bbox
    
    #also around first letter's full name equiv box
    first_letter_bbox <- bpoly2bbox(nm_poly_plot[1:nrow(nm$bpoly_info[[j]]$bpoly[[1]]),])
    first_letter_bbox$Var1 <- bpoly2bbox(nm_poly_plot)$Var1
    first_letter_bbox_plot <- first_letter_bbox
    if(plot_stuff){
      lines(first_letter_bbox_plot, lty = 2, col = 2)
    }
    
  }
  
  #draw DNA embellishment
  dh_bounds <- c(-pi/2, -pi/2 + pi * ntwists)
  tl_xy_DNA <- c(min(nm_polys[[1]][,1]), max(nm_polys[[2]][,2]))
  box_dim <- c(w = min(first_letter_bbox[,1]) - tl_xy_DNA[1], 
               h = (tl_xy_DNA[2] - min(first_letter_bbox[,2])) / (1-prop_dist))
  box_dim[1] <- box_dim[1] - box_dim[2] * prop_dist
  
  dna_bbox <- bpoly2bbox(rbind(c(tl_xy_DNA[1], tl_xy_DNA[2] - box_dim[2]), 
                               c(tl_xy_DNA[1] + box_dim[1], tl_xy_DNA[2])))
  if(plot_stuff){
    lines(dna_bbox, lty = 2, col = 2)
  }
  # box_dim[1] / (box_dim[2] * (1-prop_dist))
  # mean_iar2f
  
  #find how far to extend straight line
  extend_straight_by <- max(first_letter_bbox[,1]) - (tl_xy_DNA[1] + box_dim[1])
  if(plot_stuff){
    draw_DNA(dh_bounds = dh_bounds, amplitude = 1,
             rot = 90, strand_thickness = 0.02, 
             extend_straight = ifelse(ntwists %% 2 == 0, 2, 1), 
             extend_straight_by = extend_straight_by, col = 1, 
             topleft_xy = tl_xy_DNA, forced_box = box_dim, box_DNA = T, 
             straight_extension_real_units = T, return_info = F,
             proportional_distance = prop_dist)
  }
  
  #record information
  cell_metadata[[i]]$nm_polys <- nm_polys
  cell_metadata[[i]]$nm_polys_glyphs <- nm_polys_glyphs
  cell_metadata[[i]]$nm_polys_bbox <- nm_polys_bbox
  cell_metadata[[i]]$dna_bbox <- dna_bbox
  cell_metadata[[i]]$dna_params$dh_bounds <- dh_bounds
  cell_metadata[[i]]$dna_params$extend_straight <- ifelse(ntwists %% 2 == 0, 2, 1)
  cell_metadata[[i]]$dna_params$extend_straight_by <- extend_straight_by
  cell_metadata[[i]]$dna_params$tl_xy_DNA <- tl_xy_DNA
  cell_metadata[[i]]$dna_params$box_dim <- box_dim
  cell_metadata[[i]]$dna_params$prop_dist <- prop_dist

  
}

#### plot names ####
svg_adj <- (96/72)
#the * (96/72) is to reflect R vs Inkscape DPI mismatch

#plot the names
for(pli in 1:nplots){
  
  svg(filename = paste0("~/lab_names_grid_", pli, ".svg"), 
      width = dim_plot["w"] * svg_adj, height = dim_plot["h"] * svg_adj) 
  
  
  #initialize plotting window
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, dim_plot["w"]), ylim = c(0, dim_plot["h"]), asp = 1)
  
  #set plotting settings
  plot_outer_rect <- T
  plot_inner_rect <- F
  plot_name_rect <- F
  plot_DNA_rect <- F
  
  for(i in plot_inds[[pli]]){
    # for(i in 63:63){
    
    cat(paste0(i, " "))
    
    #retrieve name info
    name_i <- plot_name_inds[[pli]][i]
    loc_i <- plot_locs[[pli]][i]
    nm <- name_metadata[[name_i]]
    cm <- cell_metadata[[name_i]]
    disp_x <- centers[loc_i,"x"] - dim_unit["w"]/2
    disp_y <- centers[loc_i,"y"] - dim_unit["h"]/2
    
    #plot rects for troubleshooting
    if(plot_outer_rect){
      lines(cm$outer_bbox[,1] + disp_x,
            cm$outer_bbox[,2] + disp_y, col = "grey20")
    }
    
    if(plot_inner_rect){
      lines(cm$inner_bbox[,1] + disp_x,
            cm$inner_bbox[,2] + disp_y, lty = 2, col = "grey20")
    }
    
    #adjust for plotting names
    disp_x <- disp_x + name_center_adj["x"]
    disp_y <- disp_y + name_center_adj["y"]
    
    #extra adjustment to accommodate center of mass?
    if(disp_y_CoM){
      centroids <- lapply(cm$nm_polys_glyphs[[1]]$fixed_poly, st_centroid)
      areas <- sapply(cm$nm_polys_glyphs[[1]]$fixed_poly, st_area)
      centroids_matrix <- do.call(rbind, lapply(centroids, st_coordinates))
      CoM <- c(x = sum(areas * centroids_matrix[, 1]) / sum(areas),
               y = sum(areas * centroids_matrix[, 2]) / sum(areas))
      bbox_center <- sapply(apply(do.call(rbind, cm$nm_polys_glyphs[[1]]$orig_poly), 
                                  2, range, simplify = F), mean)
      yspace_left <- max(cm$inner_bbox[,2]) - max(cm$nm_polys_bbox[[1]][,2])
      extra_disp_y_CoM <- (bbox_center[2] - CoM[2]) / 4
      extra_disp_y_CoM <- min(c(0, yspace_left, extra_disp_y_CoM))
      disp_y <- disp_y + extra_disp_y_CoM
    }
    
    
    #write name
    for(j in 1:nrow(nm)){
      
      extra_disp_y_Q <- 0
      if(lower_Q){
        if(substr(nm$text[j], 1, 1) == "Q")
          extra_disp_y_Q <- -dr(cm$nm_polys_glyphs[[j]]$orig_poly[[1]][,2]) / 8
      }
      
      for(gi in 1:length(cm$nm_polys_glyphs[[j]]$fixed_poly)){
        plot(cm$nm_polys_glyphs[[j]]$fixed_poly[[gi]] + c(disp_x, disp_y + extra_disp_y_Q), col = 1, border = NA,
             add = T)
      }
      
      if(plot_name_rect){
        lines(cm$nm_polys_bbox[[j]][,1] + disp_x,
              cm$nm_polys_bbox[[j]][,2] + disp_y, col = 2, lty = 2)
      }
    }
    
    #draw DNA embellishment
    if(plot_DNA_rect){
      lines(cm$dna_bbox[,1] + disp_x,
            cm$dna_bbox[,2] + disp_y, col = 2, lty = 2)  
    }
    
    #are there gaps in the DNA strand?
    surname_poly <- t(t(cm$nm_polys[[2]]) + c(disp_x, disp_y))
    if(any(surname_poly[,2] < min(cm$dna_bbox[,2]) + disp_y)){
      locs <- find_gaps(df = surname_poly, 
                        y = min(cm$dna_bbox[,2]) + disp_y, 
                        widen_by = 0.02)
      gaps <- list(dir = "x",
                   locs = locs)  
    } else {
      gaps <- NULL
    }
    
    strand_thickness <- cm$dna_params$box_dim["h"] * 0.08
    draw_DNA(dh_bounds = cm$dna_params$dh_bounds, amplitude = 1,
             rot = 90, strand_thickness = strand_thickness, 
             extend_straight = cm$dna_params$extend_straight, 
             extend_straight_by = cm$dna_params$extend_straight_by, col = 1, 
             topleft_xy = cm$dna_params$tl_xy_DNA + c(disp_x, disp_y), 
             forced_box = cm$dna_params$box_dim, box_DNA = T,
             straight_extension_real_units = T, return_info = F,
             proportional_distance = prop_dist, gaps = gaps, take_union = T, 
             pow_thickness = 1.25, overlap_bp_strand = T, move_bp = F)
    
  }
  
  dev.off()
  
}

#### plot logos ####

#read logos in
mlab_path <- "/Users/nikgvetr/Pictures/lab_logo-metal_background_lab-not-band_monocolor.svg"
mband_path <- "/Users/nikgvetr/Pictures/MD_logo-metal_background_band_monocolor.svg"
mlab <- svgparser::read_svg(mlab_path, obj_type = 'data.frame')
mband <- svgparser::read_svg(mband_path, obj_type = 'data.frame')

#resize to the desired dimension

#identify scale factor
mlab_bbox <- bpoly2bbox(mlab[,c("x", "y")])
mband_bbox <- bpoly2bbox(mband[,c("x", "y")])
mlab_dim <- apply(mlab_bbox, 2, dr)
mband_dim <- apply(mband_bbox, 2, dr)
mlab_scale <- min(dim_logo / mlab_dim)
mband_scale <- min(dim_logo / mband_dim)

#scale coordinates
mlab$x <- mlab$x * mlab_scale
mlab$y <- mlab$y * mlab_scale
mlab_dim <- mlab_dim * mlab_scale
  
mband$x <- mband$x * mband_scale
mband$y <- mband$y * mband_scale
mband_dim <- mband_dim * mband_scale

#shift centers
mlab_buffer <- dim_unit - mlab_dim
mband_buffer <- dim_unit - mband_dim
# mlab_buffer <- c(w = 0, h = 0)
# mband_buffer <- c(w = 0, h = 0)
mlab$x <- mlab$x - min(mlab$x)# + mlab_buffer["w"] / 2
mlab$y <- mlab$y - min(mlab$y)# + mlab_buffer["h"] / 2
mband$x <- mband$x - min(mband$x)# + mband_buffer["w"] / 2
mband$y <- mband$y - min(mband$y)# + mband_buffer["h"] / 2

#get coords and sf polys from file
mband_fp <- convert_to_plot(mband)
mlab_fp <- convert_to_plot(mlab)

#the actual plotting

#first test it out in the plot pane
plot.new()
plot.window(xlim = range(mband_fp$coords$x),
            ylim = range(mband_fp$coords$y), asp = 1)
plot_polys(mband_fp$fixed_polys)
plot.new()
# plot.window(xlim = range(mlab_fp$coords$x),
#             ylim = range(mlab_fp$coords$y), asp = 1)
plot.window(xlim = range(mband_fp$coords$x),
            ylim = range(mband_fp$coords$y), asp = 1)
plot_polys(mlab_fp$fixed_polys)

#now plot to disk
montgomery_1_seen <- F
for(pli in 1:nplots){
  
  svg(filename = paste0("~/lab_logos_grid_", pli, ".svg"), 
      width = dim_plot["w"] * svg_adj, height = dim_plot["h"] * svg_adj) 
  
  
  #initialize plotting window
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, dim_plot["w"]), ylim = c(0, dim_plot["h"]), asp = 1)
  
  #set plotting settings
  plot_outer_rect <- T
  plot_inner_rect <- F
  
  for(i in plot_inds[[pli]]){
    # for(i in 63:63){
    
    cat(paste0(i, " "))
    
    #retrieve name info for logo adjustment
    name_i <- plot_name_inds[[pli]][i]
    loc_i <- plot_locs[[pli]][i]
    nm <- name_metadata[[name_i]]
    disp_x <- centers[loc_i,"x"] - dim_unit["w"]/2
    disp_y <- centers[loc_i,"y"] - dim_unit["h"]/2
    
    #plot rects for troubleshooting
    if(plot_outer_rect){
      lines(cm$outer_bbox[,1] + disp_x,
            cm$outer_bbox[,2] + disp_y, col = "grey20")
    }
    
    if(plot_inner_rect){
      lines(cm$inner_bbox[,1] + disp_x,
            cm$inner_bbox[,2] + disp_y, lty = 2, col = "grey20")
    }
    
    
    #write in logo
    if(tolower(nm$text[2]) %in% c("dahl", "md") | (tolower(nm$text[2]) == "montgomery" & montgomery_1_seen) ){
      plot_polys(mband_fp$fixed_polys, disp = c(disp_x, disp_y) + mband_buffer / 2)
    } else {
      
      #adjust for plotting logos
      disp_x <- disp_x + logo_center_adj["x"]
      disp_y <- disp_y + logo_center_adj["y"]
      
      plot_polys(mlab_fp$fixed_polys, disp = c(disp_x, disp_y) + mlab_buffer / 2)
      if(tolower(nm$text[2]) == "montgomery"){
        montgomery_1_seen <- T
      }
    }
    
    
  }
  
  dev.off()
  
}
