# Load the necessary library
library(magick)
library(imager)
library(plotwidgets)
library(clue)

# Function to create the image grid
img_grid <- function(imgs, nr, nc){
  nimg <- length(imgs)
  if(nimg < nr * nc){
    img_inds <- c(1:nimg, rep(NA, nr * nc - nimg))
  } else {
    img_inds <- 1:(nr * nc)
  }
  imgmat <- matrix(img_inds, nrow = nr, ncol = nc, byrow = TRUE)
  imgrid <- lapply(1:nr, function(ri){
    img_rinds <- imgmat[ri, !is.na(imgmat[ri,])]
    
    if(length(img_rinds) == 0){
      return(NA)
    } else {
      row_imgs <- imgs[img_rinds]
      return(image_append(image_join(row_imgs)))
    }
    
  })
  imgrid <- imgrid[!is.na(imgrid)]
  imgrid <- image_join(imgrid)
  imgrid <- image_append(imgrid, stack = TRUE)
  return(imgrid)
}

# Function to process and create the image grid from a directory
create_image_grid <- function(image_dir, grid_dims, output_file, npix = "1000", imgs = NA) {
  
  # Specify the number of rows and columns
  nr <- grid_dims[1]
  nc <- grid_dims[2]
  
  # Read in all valid image files
  if(all(is.na(imgs))){
    
    # Get all image files from the directory
    all_files <- list.files(image_dir, full.names = TRUE)
    valid_imgs <- all_files[grepl("\\.(jpg|jpeg|png|heic|gif|tiff|bmp)$", tolower(all_files))]
    imgs <- lapply(valid_imgs, image_read)  
    
    # Ensure imgs list is not empty
    if (length(imgs) == 0) {
      stop("No valid imgs to process in the specified directory.")
    }
    
    # Randomize the order of imgs
    imgs <- sample(imgs)
    
    # Rescale imgs
    imgs <- lapply(imgs, image_scale, geometry = npix)
    
  }
  
  # Create the image grid
  imgrid <- img_grid(imgs, nr, nc)
  
  # Save the grid image
  image_write(imgrid, output_file, format = "jpg")
}

average_color <- function(img) {
  img <- as.data.frame(img)
  avg_color <- colMeans(img[,1:3])
  return(avg_color)
}



#### for testing, generate random color images ####
# Load necessary library
if (!require("png")) {
  install.packages("png")
  library(png)
}

# Create the directory if it doesn't exist
output_dir <- "~/random_cols/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to generate a single color image and save it
generate_single_color_image <- function(index) {
  
  # Generate random RGB values
  color <- gtools::rdirichlet(1, alpha = c(0.5,0.5,0.5))
  
  # Create a 10x10 image with the single color
  img <- array(NA, dim = c(10, 10, 3))
  img[,,1] <- color[1]
  img[,,2] <- color[2]
  img[,,3] <- color[3]
  
  # File path
  file_path <- paste0(output_dir, index, ".png")
  
  # Save the image
  writePNG(img, target = file_path)
}

# Generate and save 1000 images
# for (i in 1:1000) {
#   generate_single_color_image(i)
# }

image_dir <- "~/Pictures/7x7_selfies/pics/"
image_dir <- "/Users/nikgvetr/pictures/band_pics/bp3"
image_dir <- "~/random_cols/"
grid_dims <- c(7, 7)  # User-specified dimensions
output_file <- "~/Pictures/7x7.jpg"
npix = "200"

# retrieve images
if(!exists("imgs") || is.null(imgs[[1]]) || !exists("retrieved_dir") || retrieved_dir != image_dir){
  
  all_files <- list.files(image_dir, full.names = TRUE)
  valid_imgs <- all_files[grepl("\\.(jpg|jpeg|png|heic|gif|tiff|bmp)$", tolower(all_files))]
  imgs <- lapply(valid_imgs, image_read)  

  # Randomize the order of imgs
  imgs <- sample(imgs)
  
  # Rescale imgs
  imgs <- lapply(imgs, image_scale, geometry = npix)
  
  # record dir for later reloading
  retrieved_dir <- image_dir
  
}

#draw grid
create_image_grid(image_dir, grid_dims, output_file, imgs = imgs)



#### now try to color match to specifix pixels in a "super" image ####

#compute color averages
super_image_path <- "/Users/nikgvetr/Downloads/rgb_squares.png"
super_image_path <- "~/Pictures/7x7_selfies/pics/IMG_1003.HEIC"
super_image <- image_read(super_image_path)
super_image_med <- image_scale(super_image, "200") # Scale to grid size
super_image_med

#fix grid dims
super_res <- dim(image_data(super_image_med))[-1]
super_AR <- super_res[2] / super_res[1]
grid_dims <- sqrt(length(imgs) / super_AR)
grid_dims <- floor(c(grid_dims * super_AR, grid_dims))
  
#resize super image to appropriate grid size
super_image_smol <- image_scale(super_image, paste0(grid_dims[2], "x", grid_dims[1])) # Scale to grid size
super_img_data <- as.numeric(image_data(super_image_smol, channels = "rgba"))

# Reshape image data to matrix
super_img_matrix <- apply(super_img_data, c(1,2), identity, simplify = F)
grid_map <- expand.grid(1:grid_dims[1], 1:grid_dims[2])
super_img_data <- do.call(rbind, super_img_matrix[t(t(grid_map))])

# Calculate average colors for sub-images
sub_img_data <- t(sapply(imgs, function(img) {
  img_rgb <- as.numeric(image_data(img, channels = "rgba"))
  return(apply(img_rgb, 3, mean))
}))

#convert to hsl
super_img_hsl <- t(rgb2hsl(t(super_img_data)))
sub_img_hsl <- t(rgb2hsl(t(sub_img_data)))

#find distances
distance_matrix <- outer(
  1:nrow(super_img_hsl), 
  1:nrow(sub_img_hsl), 
  Vectorize(function(i, j) sum((super_img_hsl[i, ] - sub_img_hsl[j, ])^2))
)

#solve linear sum assignment problem
assignment <- solve_LSAP(distance_matrix, maximum = FALSE)
assigned_pairs <- cbind(super_pix = 1:nrow(super_img_hsl), 
                        sub_img = assignment)

super_img_inds <- matrix(assigned_pairs[,2], byrow = T, nrow = grid_dims[1], ncol = grid_dims[2])
# super_img_inds[grid_dims[1]:1,]
imgs <- imgs[c(super_img_inds)]
create_image_grid(image_dir, grid_dims, output_file, imgs = imgs)
