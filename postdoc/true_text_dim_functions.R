ifelse2 <- function(test, yes, no){if(test){return(yes)}else{no}}

is_showtext_enabled <- function() {
  hooks <- getHook("plot.new")
  any(vapply(hooks, inherits, logical(1), "showtext_hook"))
}


text_distortion <- function(string,
                            font_name = "sans",
                            device = c("svg", "png")[1]){
  
  if(device == "png"){
    using_png <- T
  } else {
    using_png <- F
  }
  
  #turn on showtext to embed glyphs as polygons vs tags
  showtext_enabled <- is_showtext_enabled()
  if(!showtext_enabled){
    showtext::showtext_auto(T)
  }
  
  #add font for use with showtext if it has not been added already
  if(!(font_name %in% sysfonts::font_families())){
    font_styles <- names(formals(sysfonts::font_add))
    font_styles <- setdiff(font_styles, c("family", "symbol"))
    all_fonts <- systemfonts::system_fonts()
    matching_fonts <- subset(all_fonts, family == font_name)
    if(nrow(matching_fonts) == 0){stop("font not found on system")}
    found_styles <- gsub(" ", "", tolower(matching_fonts$style))
    matching_font_paths <- setNames(matching_fonts$path[found_styles %in% font_styles], 
                                    found_styles[found_styles %in% font_styles])
    sysfonts::font_add(family = font_name, 
                       regular = ifelse2(!is.na(matching_font_paths["regular"]), 
                                         matching_font_paths["regular"], NULL), 
                       bold = ifelse2(!is.na(matching_font_paths["bold"]), 
                                      matching_font_paths["bold"], NULL), 
                       italic = ifelse2(!is.na(matching_font_paths["italic"]), 
                                        matching_font_paths["italic"], NULL), 
                       bolditalic = ifelse2(!is.na(matching_font_paths["bolditalic"]), 
                                            matching_font_paths["bolditalic"], NULL))
  }
  
  #find reasonable plotting parameters
  str_rat <- c(w = strwidth(cex = 1, string, units = "inches", family = font_name),
               h = strheight(cex = 1, string, units = "inches", family = font_name))
  str_rat <- str_rat / str_rat[1]
  npix_max <- 1E3
  npix <- round(str_rat / max(str_rat) * npix_max)
  
  #generate plot
  dev_id <- dev.cur()  # Store current device before opening a new one
  temp_filename <- paste0(tempfile(), ifelse(using_png, ".png", ".svg"))
  if(using_png){
    png(filename = temp_filename, units = "px", width = npix["w"], height = npix["h"])
  } else {
    svg(filename = temp_filename, width = npix["w"] / 72, height = npix["h"] / 72)
  }
  
  #open plotting window
  par(mar = c(0,0,0,0))
  plot.new()
  plot.window(xlim = c(0, 1),
              ylim = c(0, 1),
              xaxs="i", yaxs="i")
  
  #initial dimension estimate
  est_dims <- c(w = strwidth(string, cex = 1, units = "user", family = font_name),
                h = strheight(string, cex = 1, units = "user", family = font_name))
  max_cex <- 1 / max(est_dims)
  cex2use <- max_cex * 0.5
  est_dims <- cex2use * est_dims
  
  #plot the text
  text_center <- c(x = 0.5, y = 0.5)
  text(x = text_center["x"],
       y = text_center["y"],
       labels = string, family = font_name,
       cex = cex2use, xpd = NA, col = 1)
  
  #close the window and device
  dev.off()
  if (dev_id %in% dev.list()) {
    dev.set(dev_id)
  }
  
  #process these files to recover the true text info
  if(using_png){
    
    #extract pixel indices
    arr <- png::readPNG(temp_filename)
    mat <- round(arr[,,1])
    coords <- which(mat==0, arr.ind = T)
    coords[,1] <- (dim(mat)[1] - coords[,1]) / npix[2]
    coords[,2] <- coords[,2] / npix[1]
    
    #find bounding polygon
    elem_idx <- 1
    yr <- range(coords[,1])
    xr <- range(coords[,2])
    true_dims <- c(w = diff(xr), h = diff(yr))
    dim_scale <- true_dims / est_dims
    true_loc <- c(x = mean(xr), y = mean(yr))
    loc_scale_disp <- (true_loc - 0.5) / true_dims
    bounding_poly <- concaveman::concaveman(as.matrix(coords), concavity = 2)[,2:1]
    
  } else {
    
    #process svg file
    svgdf <- svgparser::read_svg(temp_filename, obj_type = 'data.frame')
    elem_idx <- svgdf$elem_idx
    coords <- data.frame(x = svgdf$x, y = svgdf$y)
    bbxr <- range(svgdf$x[1:4])  # or range of x if we invert indexing
    bbyr <- range(svgdf$y[1:4])
    coords$x <- (coords$x - bbxr[1]) / diff(bbxr)
    coords$y <- 1 - (coords$y - bbyr[1]) / diff(bbyr)
    
    #remove bounding frame from coordinates
    coords <- coords[svgdf$elem_idx != 1,]
    
    #find plotting distortion effects
    xr <- range(coords$x)
    yr <- range(coords$y)
    true_dims <- c(w = diff(xr), h = diff(yr))
    dim_scale <- true_dims / est_dims * 6 / 5
    true_loc <- c(x = mean(xr), y = mean(yr))
    loc_scale_disp <- (true_loc - 0.5) / true_dims
    bounding_poly <- split(coords, svgdf$elem_idx[svgdf$elem_idx != 1])
    
  }
  
  #turn off showtext to embed glyphs as tags
  if(!showtext_enabled){
    showtext::showtext_auto(showtext_enabled)
  }
  
  return(list(
    string = string,
    font_name = font_name,
    device = device,
    dim_scale = dim_scale,
    loc_scale_disp = loc_scale_disp,
    coords = coords,
    bounding_poly = bounding_poly,
    true_loc = true_loc,
    true_dims = true_dims,
    elem_idx = elem_idx,
    svgdf = svgdf
  ))
}

true_text_params <- function(string, target_center, cex = 1, td_res = NULL,
                             font_name = "sans", device = c("svg", "png")[1]){
  
  if(device == "png"){
    using_png <- T
  } else {
    using_png <- F
  }
  
  if(is.null(td_res)){
    td_res <- text_distortion(string, font_name, device)
  }
  
  est_wh <- c(w = strwidth(string, cex = cex,
                           units = "user", family = font_name),
              h = strheight(string, cex = cex,
                            units = "user", family = font_name))
  true_wh <- est_wh * td_res$dim_scale
  disp_center_by <- td_res$loc_scale_disp * true_wh
  text_center <- target_center - disp_center_by
  
  #draw bounding box / polygon
  if(using_png){
    bpoly <- t(td_res$bounding_poly) - td_res$true_loc
    bpoly <- t(bpoly / td_res$true_dims * true_wh + target_center)
    bbox <- bpoly2bbox(bpoly)
    
  } else {
    bpoly <- t(td_res$coords) - td_res$true_loc
    bpoly <- data.frame(t(bpoly / td_res$true_dims * true_wh + target_center))
    bpoly <- split(bpoly, td_res$elem_idx[td_res$elem_idx != 1])
    cbpoly <- do.call(rbind, bpoly)
    bbox <- bpoly2bbox(cbpoly)
  }

  return(list(bpoly = bpoly,
              bbox = bbox,
              text_center = text_center,
              disp_center_by = disp_center_by))
}

bpoly2bbox <- function(bpoly){
  bbox <- expand.grid(range(bpoly[,1]), range(bpoly[,2]))
  bbox <- bbox[c(1, 3, 4, 2, 1),] #reorder as vertices
  return(bbox)
}
