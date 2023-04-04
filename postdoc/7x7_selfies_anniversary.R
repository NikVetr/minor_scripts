#libraries
library(magick)

#functions
img_grid <- function(imgs, nr, nc){
  nimg <- length(imgs)
  imgmat <- matrix(1:nimg, nrow = nr, nc = nc, byrow = T)
  imgrid <- image_join(lapply(1:nr, function(ri){
    image_append(imgs[imgmat[ri,]])
  }))
  imgrid <- image_append(imgrid, stack = T)
  return(imgrid)
}

#specify img locations
pdir <- "~/Pictures/7x7_selfies/"
picdir <- paste0(pdir, "pics/")


#specify filenames
all_files <- paste0(picdir, list.files(picdir, pattern = "\\."))
not_HEIC <- all_files[!grepl("HEIC", all_files)]

#convert to common format
not_HEICs <- sapply(not_HEIC, image_read)
sapply(1:length(not_HEICs), function(i){
  image_write(not_HEICs[[i]], gsub("\\..*", "\\.HEIC", not_HEIC), format = "HEIC")
})

#read in all HEICs
files <- paste0(picdir, list.files(picdir, pattern = "\\.HEIC"))

#randomize order, but specify some locations
files <- sample(files)
#2769 peepeeska looking down in 1,1
#2616 peepeeska looking straight in 7,7
#3296 us on boat in 4,4
swap_inds <- cbind(c(2769, 2616, 3296), 
                   c(1, 49, 25))
specific_files <- sapply(swap_inds[,1], function(x) which(grepl(x, files)))
for(i in 1:length(specific_files)){
  foo <- files[swap_inds[i,2]]
  files[swap_inds[i,2]] <- files[specific_files[i]]
  files[specific_files[i]] <- foo
}

#load in images
imgs <- sapply(files, image_read)

#rescale images
imgs <- sapply(imgs, image_scale, geometry = "1000")
imgs <- image_join(imgs)

# #this one colors the images in a weird way? no wait is it because of the jpg / heic conflict? no?
# imgrid <- image_montage(imgs, tile = '7x7',
#                         geometry = "200x200+0+0",
#                         bg = "white", shadow = F)

imgrid <- img_grid(imgs, 7, 7)
  
image_write(imgrid, paste0(pdir, "imgrid.jpg"), format = "jpg")


