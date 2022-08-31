#packages
library(foreach)
library(doParallel)
library(parallel)
library(png)
library(EBImage)

#functions
recenter <- function(img_array, zoom_ratio = 2, middle = F, bottom = F, topright = F, background_col = "red", transparent_background = F){
  img <- img_array
  if(dim(img)[3] == 3){
    img <- abind(img, matrix(1,dim(img)[1],dim(img)[2])) #add in alpha channel if none exists
  }
  
  new_img <- array(1, c(dim(img)[1:2] * zoom_ratio, 4))
  bcol <- col2rgb(background_col) / 255
  new_img[,,1] <- bcol[1]; new_img[,,2] <- bcol[2]; new_img[,,3] <- bcol[3]
  
  if(transparent_background){
    new_img[,,4] <- 0
  }
  
  if(middle){
    new_img[(dim(img)[2] / 2 + 1):(dim(img)[1] / 2 + dim(img)[1]), 
            (dim(img)[2] / 2 + 1):(dim(img)[2] / 2 + dim(img)[2]),] <- img  
  }
  
  if(bottom){
    new_img[(dim(img)[2] + 1):(dim(new_img)[2]), (dim(img)[2] / 2 + 1):(dim(img)[1] / 2 + dim(img)[1]),] <- img  
  }
  
  if(topright){
    new_img[(1):(dim(new_img)[2] - dim(img)[2]), (dim(img)[2] + 1):(dim(new_img)[1]),] <- img  
  }
  
  return(new_img)
}

#### center particular image in larger frame for masking ####
# img <- readPNG("~/Downloads/dall-e/DALL·E 2022-07-20 14.01.49 - scientific aging clock.png")
# new_img <- recenter(img, middle=T, transparent_background = T)
# writePNG(new_img, target = "~/Downloads/dall-e/test.png")

#### create a video that zooms out of a set of images ####

#load files in
base_dir <- "~/Downloads/dall-e/zoom/shambala/"
base_dir <- "~/Downloads/dall-e/zoom/ruins/"
# base_dir <- "~/Downloads/dall-e/zoom/stanford/"
frame_dir <- paste0(base_dir, "frames/")
recenter_dir <- paste0(base_dir, "recenter/")

if(!dir.exists(recenter_dir)){dir.create(recenter_dir)}


files <- paste0(base_dir, sort(list.files(base_dir, pattern = ".png")))
#remove alt frame
files <- files[!grepl("7a", files)]
# files <- files[!grepl("7.png", files)]
# files <- files[grepl("2048", files)]

d <- lapply(files, readPNG)

#write files in scaled frames
write_scaled_frames <- T
if(write_scaled_frames){
  for(i in 1:length(d)){
    for(j in 1:3){
      cat(paste0(i, " "))
      filename <- paste0(recenter_dir, gsub("\\.", paste0("_recenter\\.", j, "\\."), tail(do.call(cbind, strsplit(files, "/")), 1)[i]))
      if(file.exists(filename)){next()}
      writePNG(recenter(d[[i]], middle=j==1, bottom=j==2, topright=j==3, transparent_background =T), 
               target = filename)
    }
  }
}


#### do the actual rendering ####

if(!dir.exists(frame_dir)){dir.create(frame_dir)}
file.remove(paste0(frame_dir, list.files(path = frame_dir, pattern = ".png")))

n_pics <- length(files)
fps <- 60
ns <- 10
nf <- ns * fps
trim_prop <- 0.2
framethin <- 1
render_video <- T
incr <- 1
zoom_style <- c(1,1,1,1,2,2,1,1,1,1,1,1) #1 is center, 2 is bottom
res_scale <- 2
res <- dim(d[[1]])[1:2] * res_scale

# zoom_rats <- 2^seq(0, (n_pics-1), length.out = nf) #ie linear on log scale (so constant perceived zoom)
zoom_rats <- 2^(pbeta(1:nf / nf, 4, 4) * (n_pics-1)) #alternatively, sigmoid transformed
zoom_rats_resid <- 2^(log2(zoom_rats)%%1)

if(!exists("cl")){
  cl <- makeCluster(4, outfile="")
  registerDoParallel(cl)
  runs <- 0
}
getDoParWorkers()

foreach(fi=seq(1, nf, framethin), .packages = c("png", "EBImage")) %dopar% {
# for(fi in c(15,40)){
  
  cat(paste0(fi, " "))
  
  #accommodate special case of last frame
  if(fi == nf){
    writePNG(d[[length(d)]], target = paste0(frame_dir, paste0(rep(0, 5-nchar(((fi - 1) / framethin) + 1)), collapse = ""), ((fi - 1) / framethin) + 1,".png"))
    next()
  }
  
  
  #get transition params
  curr_trans <- log2(zoom_rats[fi])
  # prop_trans <- curr_trans %% 1
  # prop_trans <- (1 - 1 / (1 + curr_trans %% 1)) * 2
  prop_trans <- 1 - (1 / zoom_rats_resid[fi] - 0.5) * 2
  # prop_trans <- 2^(curr_trans %% 1) - 1
  pis <- c(floor(curr_trans), ceiling(curr_trans)) + 1 #NEED TO ACTUALLY USE 3X LEVELS
  if(pis[1] == pis[2]){pis[2] <- pis[2] + 1}
  pis[3] <- pis[2] + 1
  
  #load in images
  small <- resize(d[[pis[1]]], w = res[1] * (0.5 + ((1 - prop_trans) / 2)), h = res[2] * (0.5 + ((1 - prop_trans) / 2)))
  px2t_s <- round(dim(small)[1:2] * trim_prop)
  small <- small[px2t_s[1]:(dim(small)[1]-px2t_s[1]), px2t_s[2]:(dim(small)[2]-px2t_s[2]), ]
  
  big <- resize(d[[pis[2]]], w = res[1] * (2 - prop_trans), h = res[2] * (2 - prop_trans))
  px2t_b <- round(dim(big)[1:2] * trim_prop)
  orig_dim_big <- dim(big)
  
  if(any((dim(big)[1] - px2t_b) < res) & (pis[3] <= n_pics)){
    big <- big[px2t_b[1]:(dim(big)[1]-px2t_b[1]), px2t_b[2]:(dim(big)[2]-px2t_b[2]), ]
    huge <- resize(d[[pis[3]]], w = res[1] * (2 - prop_trans) * 2, h = res[2] * (2 - prop_trans) * 2)
  }
  
  if(all(dim(big)[1:2] >= res)){
    
    extra_pix <- ((dim(big)[1] - res[1]) / 2)
    miss_pix <- ((res - dim(small)[1:2]) / 2)
    
    if(zoom_style[pis[1]] == 1){
      curr_frame <- big[(extra_pix+1):(dim(big)[1]-extra_pix), (extra_pix+1):(dim(big)[2]-extra_pix),]
      curr_frame[(miss_pix[1]+1):(res[1] - miss_pix[1]), (miss_pix[1]+1):(res[1] - miss_pix[1]), ] <- small
    } else if(zoom_style[pis[1]] == 2){
      curr_frame <- big[(dim(big)[2] - res[2] + 1):(dim(big)[2]),(extra_pix+1):(dim(big)[1]-extra_pix),]
      curr_frame[(res[2] - dim(small)[2] - px2t_s[2] + 1):(res[2] - px2t_s[2]), (miss_pix[1]+1):(res[1] - miss_pix[1]), ] <- small
    }
    
  } else {
    
    extra_pix <- ((dim(huge)[1] - res[1]) / 2)
    miss_pix_b <- ((res - dim(big)[1:2]) / 2)
    miss_pix_s <- ((res - dim(small)[1:2]) / 2)
    
    if(zoom_style[pis[1]] == 1){
      
      if(zoom_style[pis[2]] == 1){
        curr_frame <- huge[(extra_pix+1):(dim(huge)[1]-extra_pix), (extra_pix+1):(dim(huge)[2]-extra_pix),]  
      } else if(zoom_style[pis[2]] == 2){
        curr_frame <- huge[(dim(huge)[2] + 1 + miss_pix_b[2] - px2t_b[2] - res[2]):(dim(huge)[2] + 1 + miss_pix_b[2] - px2t_b[2]),
                           (extra_pix+1):(dim(huge)[1]-extra_pix),]
      }
      curr_frame[(miss_pix_b[1]+1):(res[1] - miss_pix_b[1]), (miss_pix_b[1]+1):(res[1] - miss_pix_b[1]), ] <- big
      curr_frame[(miss_pix_s[1]+1):(res[1] - miss_pix_s[1]), (miss_pix_s[1]+1):(res[1] - miss_pix_s[1]), ] <- small
      
    } else if(zoom_style[pis[1]] == 2){
      
      if(zoom_style[pis[2]] == 2){
        curr_frame <- huge[(dim(huge)[2] - res[2] + 1):(dim(huge)[2]),(extra_pix+1):(dim(huge)[1]-extra_pix),]
      } else if(zoom_style[pis[2]] == 1){
        curr_frame <- huge[(dim(huge)[2]-(dim(huge)[2] - dim(big)[2] - px2t_b[2])/2 - res[2]):(dim(huge)[2]-(dim(huge)[2] - dim(big)[2] - px2t_b[2])/2), 
                           (extra_pix+1):(dim(huge)[1]-extra_pix),]  
        curr_frame <- huge[(dim(huge)[2] - orig_dim_big[2] / 2 - res[2]):(dim(huge)[2] - orig_dim_big[2] / 2), 
                           (extra_pix+1):(dim(huge)[1]-extra_pix),]  
      }
      curr_frame[(res[2] - dim(big)[2] - px2t_b[2] + 1):(res[2] - px2t_b[2]), (miss_pix_b[1]+1):(res[1] - miss_pix_b[1]), ] <- big
      curr_frame[(res[2] - dim(small)[2] - px2t_s[2] + 1):(res[2] - px2t_s[2]), (miss_pix_s[1]+1):(res[1] - miss_pix_s[1]), ] <- small
    }
    
  }
  
  
  curr_frame <- resize(curr_frame, w = res[1] / res_scale, h = res[2] / res_scale)
  writePNG(curr_frame, target = paste0(frame_dir, paste0(rep(0, 5-nchar(((fi - 1) / framethin) + 1)), collapse = ""), ((fi - 1) / framethin) + 1,".png"))
  
  ####
  
}

runs <- runs + incr
if(runs > 5){
  stopCluster(cl = cl)
  rm(cl)
}

if(render_video){
  #stitch frames together
  base_filename <- paste0("dalle_zoom_", res[1] / res_scale, "x" , res[2] / res_scale)
  raw_filename <- paste0("raw_", base_filename, ".mp4")
  final_filename <- paste0("final_", base_filename, ".mp4")
  
  if(file.exists(paste0(base_dir, raw_filename))){file.remove(paste0(base_dir, raw_filename))}
  if(file.exists(paste0(base_dir, "rev_", raw_filename))){file.remove(paste0(base_dir, "rev_", raw_filename))}
  if(file.exists(paste0(base_dir, final_filename))){file.remove(paste0(base_dir, final_filename))}
  if(file.exists(paste0(base_dir, "alt_", final_filename))){file.remove(paste0(base_dir, "alt_", final_filename))}
  if(file.exists(paste0(base_dir, "smaller_alt_", final_filename))){file.remove(paste0(base_dir, "smaller_alt_", final_filename))}
  
  system(paste0("cd ", base_dir,"; ffmpeg -r ", fps / framethin," -f image2 -s ", res[1] / res_scale, "x" , res[2] / res_scale," -i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ", raw_filename))
  system(paste0("cd ", base_dir,"; ", 
                "ffmpeg -i ", raw_filename, " -vf reverse rev_", raw_filename,"; ",
                "touch input.txt;",
                "echo \"file ", raw_filename,"\nfile rev_", raw_filename,"\" > input.txt;",
                "ffmpeg -f concat -i input.txt -codec copy ", final_filename))
  system(paste0("cd ", base_dir,"; ", 
                "touch input.txt;",
                "echo \"file rev_", raw_filename,"\nfile ", raw_filename,"\" > input.txt;",
                "ffmpeg -f concat -i input.txt -codec copy alt_", final_filename))
  system(paste0("cd ", base_dir,"; ",
                "ffmpeg -i alt_", final_filename, " -vf scale=", res[1] / res_scale / 2,":", res[2] / res_scale / 2," -preset slow -crf 24 smaller_alt_", final_filename))

  # system(paste0("cd ", base_dir,"; ", 
  #               "ffmpeg -i alt_", final_filename, " -vcodec libx265 -crf 28 x265_alt_", final_filename))
  # system(paste0("cd ", base_dir,"; ", 
  #               "ffmpeg -i alt_", final_filename, " -vcodec libx264 -crf 35 x264_alt_", final_filename))
  
  
  
}
