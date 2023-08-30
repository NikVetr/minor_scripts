# functions
text2 <- function(x, y, pos = NULL, cex = 1, labels = NULL, drect = F, ...){
  adj_x <- x + ifelse(any(pos %in% c(2,4)), ifelse(any(pos == 2), -1, 1) * strwidth(labels, cex = cex) / 2, 0)
  adj_y <- y + ifelse(any(pos %in% c(1,3)), ifelse(any(pos == 1), -1, 1) * strheight(labels, cex = cex) / 2, 0)
  text(x = adj_x, y = adj_y, labels = labels, pos = NULL, cex = cex, ...)
  if(drect){
    rect(xleft = adj_x - strwidth(labels, cex = cex) / 2, 
         xright = adj_x + strwidth(labels, cex = cex) / 2, 
         ybottom = adj_y - strheight(labels, cex = cex) / 2, 
         ytop = adj_y + strheight(labels, cex = cex) / 2)
    # abline(h = adj_y - strheight(labels, cex = cex) / 2, lwd = 0.5)
  }
}

extract_table <- function(rlines){
  rlines <- do.call(rbind, lapply(lapply(strsplit(rlines, "\\|"), as.data.frame), t))
  rlines <- trimws(rlines)
  colnames(rlines) <- gsub(" ", "_", rlines[1,])
  rownames(rlines) <- NULL
  rlines <- rlines[,-1]
  rlines <- rlines[-1,]
  rlines <- as.data.frame(rlines)
}

#make folders
wanderings_dir <- "~/wanderings/"
wanderings_fonts_dir <- paste0(wanderings_dir, "fonts/")
wanderings_controlnet_words <- paste0(wanderings_dir, "controlnet_words/")
dloaded_fonts_dir <- paste0(wanderings_dir, "dloaded_fonts/")
clean_fonts_dir <- paste0(wanderings_dir, "clean_fonts/")

if(!dir.exists(dloaded_fonts_dir)){dir.create(dloaded_fonts_dir)}
if(!dir.exists(wanderings_dir)){dir.create(wanderings_dir)}
if(!dir.exists(wanderings_fonts_dir)){dir.create(wanderings_fonts_dir)}
if(!dir.exists(wanderings_controlnet_words)){dir.create(wanderings_controlnet_words)}
if(!dir.exists(clean_fonts_dir)){dir.create(clean_fonts_dir)}

#read in initial descriptions
txt <- readLines("~/Documents/Documents - nikolai/wanderings.txt", warn = F)
txt <- sapply(sapply(strsplit(txt, "\\."), tail, 1), trimws)
txt <- lapply(txt, strsplit, ":")
txt <- data.frame(name = sapply(lapply(txt, unlist), head, 1), 
           desc = sapply(lapply(txt, unlist), tail, 1))
txt$desc <- trimws(txt$desc)
common <- "by greg rutkowski, by caspar david friedrich, detailed, 4k, award-winning"
txt$word <- tolower(gsub("Wanderer above the Sea of ", "", txt$name))
rownames(txt) <- NULL

#also extra descriptions
edescr <- tolower(readLines("~/Documents/Documents - nikolai/wanderings_extra_descr", warn = F))
edescr <- edescr[edescr != ""]
edescr <- strsplit(edescr, ":")
edescr <- lapply(lapply(edescr, as.matrix), t)
edescr <- as.data.frame(do.call(rbind, edescr))
colnames(edescr) <- c("word", "description")
edescr$word <- trimws(gsub(".+\\.", "", edescr$word))
edescr$description <- trimws(gsub("\\.", "", edescr$description))
edescr <- edescr[match(txt$word, edescr$word),]

#extract job
cat(paste0(edescr$word, ": ", 
       sapply(strsplit(edescr$description, ","), function(x) trimws(x[3])), "\n"))

# objects
things <- tolower(readLines("~/Documents/Documents - nikolai/wanderings_held_objects.txt", warn = F))
things <- extract_table(things)
things <- things[match(txt$word, things$term),]

#person descriptions
colpers <- readLines("~/Documents/Documents - nikolai/wanderings_colors_person", warn = F)
colpers <- extract_table(colpers)
colpers$Dominant_Colors <- tolower(colpers$Dominant_Colors)
colpers$Theme <- tolower(colpers$Theme)
colpers <- colpers[match(txt$word, colpers$Theme),]

#construct input to SD (+ prompt)
txt$pos_descr <- paste0(txt$name, ", ", 
                        edescr$description, ", ", 
                        paste0("(", tolower(txt$word), ")"), ", ", 
                        colpers$Dominant_Colors, ", ", 
                        colpers$Figure, ", ", 
                        paste0("figure holding long ", things$object), ", ",
                        common)
txt <- txt[match(unique(txt$name), txt$name),]

#merge the two
txt <- cbind(txt, fonts[,-which(colnames(fonts) == "Unique_Word")])

#dload and match fonts
# fonts <- readLines("~/Documents/Documents - nikolai/wanderings_fonts_matched.txt", warn = F)
# fonts <- extract_table(fonts)
# fonts <- readLines("~/Documents/Documents - nikolai/wanderings_realfonts.txt", warn = F)
font_info <- data.table::fread("~/dafont_top1000.txt")
fonts <- readLines("~/Documents/Documents - nikolai/wanderings_font_multimatch.txt", 
                   warn = F)
fonts <- extract_table(fonts)
fonts <- fonts[trimws(tolower(fonts$Font_Name)) %in% trimws(tolower(font_info$font_names)),]
fonts$n <- as.numeric(table(fonts$Font_Name)[fonts$Font_Name])

#find text similarity to rank fonts

#first read in or generate latent embeddings
library(text)
if(!file.exists("~/wanderings-keyword-embeddings.txt") || 
   any(rownames(read.csv("~/wanderings-keyword-embeddings.txt", row.names = 1)) != txt$word)){
  txt_descr_embeddings <- lapply(1:nrow(txt), function(i){
    system(sprintf('echo "%s "', paste0(i, collapse="")))
    textEmbed(paste0(txt$word[i], ", ", txt$desc[i]))$x
  })
  txt_descr_embeddings <- as.data.frame(do.call(rbind, txt_descr_embeddings))
  rownames(txt_descr_embeddings) <- txt$word
  write.csv(txt_descr_embeddings, "~/wanderings-keyword-embeddings.txt", row.names = T)
} else {
  txt_descr_embeddings <- read.csv("~/wanderings-keyword-embeddings.txt", row.names = 1)
  txt_descr_embeddings_tb <- tibble::tibble(txt_descr_embeddings)
}
font_embeddings <- read.csv("~/dafont_top1000_name-embeddings.txt", row.names = 1)
font_embeddings_tb <- tibble::tibble(font_embeddings)

#then find similarities... can do dim reduction if too slow
similarities <- parallel::mclapply(1:nrow(txt_descr_embeddings_tb), function(i){
  system(sprintf('echo "%s "', paste0(i, collapse="")))
  text::textSimilarity(txt_descr_embeddings_tb[rep(i, nrow(font_embeddings_tb)),], font_embeddings_tb)
}, mc.cores = 12)
similarities <- do.call(rbind, similarities)
rownames(similarities) <- rownames(txt_descr_embeddings)
colnames(similarities) <- rownames(font_embeddings)

#top 3 fonts for each term
lapply(setNames(rownames(similarities), rownames(similarities)), function(i){
  tail(sort(similarities[i,]), 3)
})

#assign fonts
fonts <- split(fonts, fonts$Term)
fonts <- do.call(rbind, lapply(fonts, function(x){
  scores <- as.numeric(x$Score)
  top_scores <- tail(sort(unique(scores)), 2)
  xtop <- x[scores %in% top_scores,]
  xtop <- xtop[order(xtop$n),]
  xtop$rank <- rank(100-as.numeric(xtop$Score), ties.method = "first")
  # xtop[which.min(xtop$Score),]
  xtop[!duplicated(xtop),]
}))
fonts$Term <- tolower(trimws(gsub(".+\\.", "", fonts$Term)))

# fonts <- fonts[match(txt$word, tolower(fonts$Term)),]
fonts$font_id <- setNames(font_info$font_id, font_info$font_names)[fonts$Font_Name]
fonts$font_id <- gsub("-", "_", fonts$font_id)
fonts$font_url <- setNames(font_info$dload_link, font_info$font_names)[fonts$Font_Name]
fonts$font_url <- gsub("-", "_", fonts$font_url)
fonts_to_dload <- !(fonts$font_id %in% gsub(".zip", "", list.files(wanderings_fonts_dir)))
cat(paste0("wget --content-disposition ", 
           paste0("\"", unique(fonts$font_url[fonts_to_dload]), "\"", collapse = " ")))

font_filepaths <- list.files(wanderings_fonts_dir, full.names = T)
unzips <- lapply(setNames(font_filepaths, gsub(".zip", "", basename(font_filepaths))), function(x){
  out <- try(unzip(x, exdir = dloaded_fonts_dir))
  setNames(list(out), basename(x))
})
unzips <- lapply(unzips, function(i) data.frame(font = rep(gsub(".zip", "", names(i)), length(i[[1]])), file = as.character(i[[1]])))
unzips <- do.call(rbind, unzips)
unzips[startsWith(unzips$file, "Error"),]
unzips$file[unzips$font == "birds_of_paradise"] <- "/Users/nikgvetr/wanderings/dloaded_fonts//Birds of Paradise И PERSONAL USE ONLY.ttf"
unzips$file[unzips$font == "moms_typewriter"] <- "/Users/nikgvetr/wanderings/dloaded_fonts//Mom差___.ttf"
unzips$ext <- sapply(strsplit(unzips$file, "\\."), tail, 1)

unique(unzips$ext)
rownames(unzips) <- NULL
unzips <- unzips[unzips$ext %in% c("ttf", "otf", "woff", "woff2", "TTF", "rtf"),]
unzips$file[unzips$ext %in% c("woff", "woff2")] <- gsub("woff|woff2", "ttf", unzips$file[unzips$ext %in% c("woff", "woff2")])
unzips <- unzips[match(unique(unzips$file), unzips$file),]
file.copy(unzips$file, clean_fonts_dir, overwrite = T)

unzips$fontnames <- sapply(unzips$file, function(x){
  sysout <- system(paste0("mdls ", "\"", x, "\""), intern = TRUE)
  fullname_line <- sysout[grep("com_apple_ats_name_full", sysout) + 1]
  fullname <- gsub(".*\" ", "", fullname_line)
  fullname <- gsub("\"", "", fullname)
  trimws(fullname)
})
fonts$family <- setNames(unzips$fontnames, unzips$font)[fonts$font_id]


#original approach to doing all this
dload_fonts <- F
if(dload_fonts){
  
  
  # download.file(txt$Estimated_Download_Link[1], 
  #               destfile = paste0(wanderings_fonts_dir, 
  #                                 gsub(".+?f=", "", txt$Estimated_Download_Link[1]),
  #                                 ".zip"))
  
  # wget --content-disposition
  # cat(txt$Estimated_Download_Link)
  # scp -r nikgvetr@smsh11dsu-srcf-d15-38.scg.stanford.edu:/oak/stanford/groups/smontgom/nikgvetr/fonts/ ~/wanderings/fonts/
  
  # check which fonts correctly downloaded
  dloaded_fonts <- list.files(wanderings_fonts_dir)
  failed <- dloaded_fonts[grepl("index.html", dloaded_fonts)]
  failed_fonts <- txt$Suggested_Font[match(gsub(".+f=", "", failed), 
                                           gsub(".+f=", "", txt$Estimated_Download_Link))]
  cat(paste0(failed_fonts, sep = ", "))
  
  #snag alternates
  alt_fonts <- readLines("~/Documents/Documents - nikolai/wanderings_fonts_matched_alts.txt", 
                         warn = F)
  alt_fonts <- do.call(rbind, lapply(lapply(strsplit(alt_fonts, "\\|"), as.data.frame), t))
  alt_fonts <- trimws(alt_fonts)
  colnames(alt_fonts) <- gsub(" ", "_", alt_fonts[1,])
  rownames(alt_fonts) <- NULL
  alt_fonts <- alt_fonts[,-1]
  alt_fonts <- alt_fonts[-1,]
  alt_fonts <- as.data.frame(alt_fonts)
  # cat(paste0("https://dl.dafont.com/dl/?f=", unlist(alt_fonts[,-1])))
  
  #check which ones made it through
  dloaded_altfonts <- list.files(wanderings_altfonts_dir)
  failed_altfonts <- dloaded_altfonts[grepl("index.html", dloaded_altfonts)]
  success_altfonts <- setdiff(dloaded_altfonts, failed_altfonts)
  failed_altfonts <- gsub(".+f=", "", failed_altfonts)
  alt_found <- apply(alt_fonts, 2, function(x) 
    paste0(x, ".zip") %in% success_altfonts)
  font_ind_found <- apply(alt_found, 1, function(x) min(which(x)))
  font_ind_found[font_ind_found == Inf] <- NA
  alt_fonts$success <- sapply(1:nrow(alt_fonts), function(i){
    alt_fonts[i, font_ind_found[i]]
  })
  alt_fonts$success[sapply(alt_fonts$success, is.null)] <- NA
  alt_fonts$success <- unlist(alt_fonts$success)
  
  #subset to only those that made it through
  txt <- txt[!(txt$Suggested_Font %in% alt_fonts$Original_Font[is.na(alt_fonts$success)]),]
  orig_success <- gsub(".zip", "", dloaded_fonts[grepl(".zip", dloaded_fonts)])
  
  #match up fonts to words
  txt$font_id <- gsub(".+f=", "", txt$Estimated_Download_Link)
  txt$font_to_use <- txt$font_id
  txt$font_to_use[!(txt$font_id %in% orig_success)] <- NA
  txt$font_to_use[is.na(txt$font_to_use)] <- setNames(alt_fonts$success, alt_fonts$Original_Font)[txt$Suggested_Font[is.na(txt$font_to_use)]]
  
  #open and install fonts
  all_files <- list.files(c(wanderings_altfonts_dir, wanderings_fonts_dir), full.names = T)
  font_paths <- sapply(paste0(txt$font_to_use, ".zip"), function(fi)
    all_files[grepl(fi, all_files)][1]
  )
  file.copy(font_paths, dloaded_fonts_dir, overwrite = T)
  sapply(list.files(dloaded_fonts_dir, full.names = T), unzip, exdir = dloaded_fonts_dir)
  font_filepaths <- list.files(dloaded_fonts_dir, recursive = T, full.names = T)
  font_files <- basename(font_filepaths)
  font_endings <- gsub(".+\\.", "", basename(list.files(dloaded_fonts_dir, recursive = T)))
  unique(font_endings)
  
  
  file.copy(font_filepaths[font_endings %in% c("ttf", "TTF", "otf", "rtf")], 
            clean_fonts_dir, overwrite = T)
}

#### making the pictures ####

library(extrafont)
library(extrafontdb)

#import fonts
# font_import()
# loadfonts()

#do the plotting
dims <- c(900, 1152)
font_height <- 0.1
extra_height_prop <- 0.15
full_dims <- round(dims * c(1, 1 + extra_height_prop))

for(i in 1:nrow(txt)){
  
  #strheight and strwidth are unreliable, so we plot the image first
  #then read it back in and adjust as needed
  cat(paste0(i, ", "))
  
  #read in corresponding fonts
  word <- txt$word[i]
  word_fonts <- unique(fonts$family[fonts$Term == word])
  
  for(j in 1:length(word_fonts)){
    
    #retrieve font
    font <- word_fonts[j]
    
    #make base plot
    imgpath <- paste0(wanderings_controlnet_words, 
                      paste0(rep(0, 5 - nchar(i)), collapse = ""), 
                      i, ".", j, "-", txt$word[i], ".png")
    png(filename = imgpath, width = full_dims[1], height = full_dims[2])
    par(mar = c(0,0,0,0))
    plot(1,1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", 
         frame = F, col = "white", xlim = c(0,1), ylim = c(0,1))
    
    #calc word params
    strd <- c(w = strwidth(word, units = "user", family = font, cex = 1), 
              h = strheight(word, units = "user", family = font, cex = 1))
    space_available <- c(w = 1 * 0.9, h = font_height * 1.08)
    cex <- min(space_available / strd) / 2
    
    #plot word
    text(labels = word, family = font, x = 0.5, y = 0.5, cex = cex, xpd = NA)
    
    dev.off()
    
    #read word back in to replot
    img <- png::readPNG(imgpath)
    img_pix <- which(img[,,1] < 0.5, arr.ind = T)
    
    #delete if font failed
    if(nrow(img_pix) < 50){
      file.remove(imgpath)
      next()
    }
    
    # Find the size of the bounding box for the text
    xdim <- range(img_pix[,2]) / ncol(img)
    ydim <- range(img_pix[,1]) / nrow(img)
    txts <- c(w = diff(xdim), h = diff(ydim))
    
    #find new params
    cex_err <- min(space_available / txts)
    cex <- cex_err * cex
    xloc <- 0.5 + c(0.5 - mean(xdim)) * cex_err
    yloc <- 0.5 + c(0.5 - mean(ydim)) * cex_err + 0.52 - 
      txts["h"] * cex_err / 2 - (extra_height_prop - font_height) / 2
    
    png(filename = imgpath, width = full_dims[1], height = full_dims[2])
    par(mar = c(0,0,0,0))
    plot(1,1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", 
         frame = F, col = "white", xlim = c(0,1), ylim = c(0,1))
    
    #plot the final image
    text(labels = word, family = font, x = xloc, y = yloc, cex = cex, xpd = NA)  
    
    dev.off()
  }
  
}

#evaluate successful completions
completed_txt <- txt[txt$word %in% gsub(".+-|.png.*", "", list.files(wanderings_controlnet_words)),]
writeLines(completed_txt$pos_descr, "~/wanderings/wanderings_positive_prompts.txt")
zip(zipfile = path.expand("~/wanderings/controlnet_words.zip"), 
    files = path.expand("~/wanderings/controlnet_words/"))

zip(zipfile = ,
    files = )
