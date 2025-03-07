library(sf)
library(igraph)
library(extrafont)
library(extrafontdb)

#functions
source("~/repos/polylines/R/functions.R")
source("~/scripts/minor_scripts/postdoc/true_text_dim_functions.R")
source("~/scripts/minor_scripts/postdoc/rect_constraint.R")

#import fonts
# font_import()
# loadfonts()
names(pdfFonts())[grepl("oboto", names(pdfFonts()))]

#render text with polygons, not svg tags
showtext::showtext_auto(T)

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
degree_name_overlap <- 0.4
name_cex_ratio_base <- 0.55
name_cex_ratio_factor <- 1.1
cex_scale <- 0.95
vws_scale <- 0.6
vws_adj_yqpgj <- 1.5
twist_ratio_xy <- 3
prop_dist <- 0.15

#construct matrix of centers
centers <- as.matrix(expand.grid(1:dim_grid[1], 1:dim_grid[2])) - 1
centers <- centers %*% diag(dim_unit + dim_spacing)
centers <- t(t(centers) + start_loc)
colnames(centers) <- c("x", "y")

#specify font information
font_name <- "Source Serif Pro"

#retrieve true location / bounding info for names
if(!exists("name_info")){
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
    names(name_info) <- sapply(labnames_list, paste0, collapse = " ")
    save(name_info, file = name_info_path)  
  }
}

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

#let's now reprocess these metadata in the context of their true metadata
if(!("bpoly_info" %in% names(name_metadata[[1]]))){
    
  temp_filename <- paste0(tempfile(), ".svg")
  svg(filename = temp_filename, width = dim_plot["w"], height = dim_plot["h"])
  
  par(mar = c(0,0,0,0))
  plot(1,1, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = F,
       xlim = c(0, dim_plot["w"]),
       ylim = c(0, dim_plot["h"])
  )
  
  dna_info <- draw_DNA(dh_bounds = pi/2*c(-1,1), amplitude = 1,
                       rot = 90, strand_thickness = 0.01, 
                       extend_straight = 1, 
                       extend_straight_by = 0, col = 1, 
                       topleft_xy = c(0,0), box_DNA = F, 
                       return_info = T,
                       proportional_distance = prop_dist)
  
  for(i in 1:prod(dim_grid)){
  
    cat(paste0(i, " "))
    if(i <= length(labnames_list)){
      nm <- name_metadata[[i]]
      for(j in 1:length(labnames_list[[i]])){
        bpoly_info <- true_text_params(string = nm$text[j], 
                                       target_center = c(centers[i,1] + nm$dx[j], centers[i,2] + nm$dy[j]), 
                                       cex = nm$cex[j], device = "svg", font_name = font_name)
        name_metadata[[i]]$bpoly_info[j] <- list(bpoly_info)
        
        text(x = centers[i,1] + nm$dx[j], 
             y = centers[i,2] + nm$dy[j], 
             labels = nm$text[j], family = font_name, 
             cex = nm$cex[j], xpd = NA, col = 1)
        
      }
    }
  }
  
  dev.off()

}

#using these as initial values, let's optimize the placement of graphical elements
plot_stuff <- F
b2_h <- max(sapply(dna_info$bp, function(bp) diff(range(bp$y))))
b2_w <- diff(range(dna_info$s1[[1]]$x))
cell_metadata <- lapply(1:prod(dim_grid), function(i) NULL)
for(i in 1:prod(dim_grid)){
  
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
  mean_iar2f <- b2_w / b2_h # inverse aspect ratio for b2 (w / h), corresponding to the DNA strand top nucs
  
  #optimize placement of elements in fully generic sense
  rectangles <- place_rectangles_bfgs(
    prop_dist   = prop_dist,
    total_width = total_width,
    total_height = total_height,
    ar1         = ar1,
    ar3         = ar3,
    start_w1    = 0.5,
    start_w2    = 0.4,
    start_h2    = 0.25,
    big_penalty = 1e2,
    ratio_mean_b1b2_w  = 2.5,
    ratio_sd_b1b2_w    = 1,
    ratio_mean_b1b2_h  = 2.5,
    ratio_sd_b1b2_h    = 0.25,
    ratio_mean_b2b3_w = 1,
    ratio_sd_b2b3_w = 2,
    ratio_mean_b1_b2b3_w = 0.925,
    ratio_sd_b1_b2b3_w = 0.025,
    mean_iar2f = mean_iar2f,
    sd_iar2f = 0.2,
    iar2f_weight = 0.5,
    two_stage_iar2f = T, 
    nrep = 50
  )
  rect_info <- get_rect_info(rectangles)
  ntwists <- round(rect_info$ar2_factor)
  if(plot_stuff){
    plot_rectangles(rectangles)
  }
  
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
  transform_points <- function(x, orig_span, orig_range, relative_range) {
    new_range <- orig_span[1] + relative_range * diff(orig_span)
    mapped_x <- new_range[1] + (x - orig_range[1]) * (diff(new_range) / diff(orig_range))
    return(mapped_x)
  }
  dr <- function(x) diff(range(x))
  drrat <- function(x) dr(x[,1]) / dr(x[,2])
  if(i <= length(labnames_list)){
    
    #write name
    nm_polys <- list()
    nm_polys_glyphs <- list(NULL, NULL)
    nm_polys_bbox <- list()
    for(j in 1:length(labnames_list[[i]])){
      
      nm_poly <- do.call(rbind, nm$bpoly_info[[j]]$bpoly)
      nm_poly_plot <- cbind(nm_poly[,1] - centers[i,"x"] + dim_unit["w"]/2, 
                            nm_poly[,2] - centers[i,"y"] + dim_unit["h"]/2)
      
      if(resc_rect){
        resc_poly <- nm_poly_plot
        if(j==1){
          resc_poly[,1] <- transform_points(x = resc_poly[,1],
                                            orig_span = range(inner_bbox[,1]),
                                            orig_range = range(resc_poly[,1]),
                                            relative_range = c((1 - rect_info$total_width) / 2,
                                                               (1 - rect_info$total_width) / 2 + rectangles$w1))
          resc_poly[,2] <- transform_points(x = resc_poly[,2],
                                            orig_span = range(inner_bbox[,2]),
                                            orig_range = range(resc_poly[,2]),
                                            relative_range = c((total_height - (total_height - rect_info$total_height) / 2 - rectangles$h1 + rectangles$h2 * prop_dist / 2) / total_height,
                                                               (total_height - (total_height - rect_info$total_height) / 2 + rectangles$h2 * prop_dist / 2) / total_height
                                                                ))
        } else {
          resc_poly[,1] <- transform_points(x = resc_poly[,1],
                                            orig_span = range(inner_bbox[,1]),
                                            orig_range = range(resc_poly[,1]),
                                            relative_range = c(
                                                               1 - (1 - rect_info$total_width) / 2 - rectangles$w3,
                                                               1 - (1 - rect_info$total_width) / 2
                                                               )
                                            )
          diff(c((1 - rect_info$total_width) / 2,
            (1 - rect_info$total_width) / 2 + rectangles$w3))
          resc_poly[,2] <- transform_points(x = resc_poly[,2],
                                            orig_span = range(inner_bbox[,2]),
                                            orig_range = range(resc_poly[,2]),
                                            relative_range = c(
                                                               ((total_height - rect_info$total_height) / 2 + rectangles$h2 * prop_dist / 2) / total_height,
                                                               ((total_height - rect_info$total_height) / 2 + rectangles$h3 + rectangles$h2 * prop_dist / 2) / total_height
                                                               )
                                            )
        }
        #record the overall polygon
        nm_poly_plot <- nm_polys[[j]] <- resc_poly
        
        #split the polygon into individual glyphs (using the original poly indices)
        poly2glyph <- function(concat_poly, orig_poly){
          glyph_inds <- rep(1:length(orig_poly), sapply(orig_poly, nrow))
          return(split(data.frame(concat_poly), glyph_inds))
        }
        nm_polys_glyphs[[j]]$orig_poly <- poly2glyph(concat_poly = resc_poly, 
                                                     orig_poly = nm$bpoly_info[[j]]$bpoly)
        
        #convert these to valid polygons for plotting
        fix_polys <- function(glyphs){
          lapply(glyphs, function(glyph){
            glyph_poly <- sf::st_polygon(list(as.matrix(rbind(glyph, glyph[1,]))))
            glyph_poly <- st_make_valid(glyph_poly)
            glyph_sf <- st_sfc(glyph_poly)
            return(glyph_sf)
          })
        }
        nm_polys_glyphs[[j]]$fixed_poly <- fix_polys(glyphs = nm_polys_glyphs[[j]]$orig_poly)
      }
      
      #plot name
      if(plot_stuff){
        lines(nm_poly_plot)
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
      lines(dna_bbox, lty = 2, col = "grey20")
    }
    
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
  
}

#### plot ####

#now do the plotting
svg(filename = "~/lab_names_grid.svg", width = dim_plot["w"], height = dim_plot["h"])

#initiate plotting window
par(mar = c(0,0,0,0))
plot(1,1, col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = F,
     xlim = c(0, dim_plot["w"]),
     ylim = c(0, dim_plot["h"])
)

plot_outer_rect <- T
plot_inner_rect <- F
plot_name_rect <- F
plot_DNA_rect <- F
for(i in 1:64){
# for(i in 63:63){
    
  cat(paste0(i, " "))
  
  #retrieve name info
  nm <- name_metadata[[i]]
  cm <- cell_metadata[[i]]
  disp_x <- centers[i,"x"] - dim_unit["w"]/2
  disp_y <- centers[i,"y"] - dim_unit["h"]/2
  
  #plot rects for troubleshooting
  if(plot_outer_rect){
    lines(cm$outer_bbox[,1] + disp_x,
          cm$outer_bbox[,2] + disp_y, col = "grey20")
  }
  
  if(plot_inner_rect){
    lines(cm$inner_bbox[,1] + disp_x,
          cm$inner_bbox[,2] + disp_y, lty = 2, col = "grey20")
  }
  

  #plot names
  if(i <= length(labnames_list)){
    
    #write name
    for(j in 1:length(labnames_list[[i]])){
      for(gi in 1:length(cm$nm_polys_glyphs[[j]]$fixed_poly)){
        plot(cm$nm_polys_glyphs[[j]]$fixed_poly[[gi]] + c(disp_x, disp_y), col = 1, border = NA,
             add = T)
      }
      # text(x = nm$bpoly_info[[j]]$text_center[["x"]], 
      #      y = nm$bpoly_info[[j]]$text_center[["y"]], 
      #      labels = nm$text[j], family = font_name, 
      #      cex = nm$cex[j], xpd = NA, col = 1)
      
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
      locs <- find_gaps(df = surname_poly, 
                        y = min(cm$dna_bbox[,2]) + disp_y, 
                        widen_by = 0.02)
      gaps <- list(dir = "x",
                   locs = locs)  
    } else {
      gaps <- NULL
    }
    
    
    draw_DNA(dh_bounds = cm$dna_params$dh_bounds, amplitude = 1,
             rot = 90, strand_thickness = 0.01, 
             extend_straight = cm$dna_params$extend_straight, 
             extend_straight_by = cm$dna_params$extend_straight_by, col = 1, 
             topleft_xy = cm$dna_params$tl_xy_DNA + c(disp_x, disp_y), 
             forced_box = cm$dna_params$box_dim, box_DNA = T,
             straight_extension_real_units = T, return_info = F,
             proportional_distance = prop_dist, gaps = gaps)
    
  
  }
  
}

dev.off()

