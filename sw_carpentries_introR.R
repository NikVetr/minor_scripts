#sw carpentries intro to R

#check your current working directory
getwd()

#change working directory to somewhere appropriate, with path specified relative to home directory
setwd("~/SWC_Workshop") #relative to home directory
setwd(dir = "/Volumes/macOS/Users/nikolai/SWC_Workshop") #absolute path

#create a new directory to hold the data for this tutorial, and change the current working directory to it
dir.create("r-novice-inflammation")
setwd("r-novice-inflammation/")

#download the data to your new folder, naming it "r-novice-inflammation-data"
download.file(url = "https://swcarpentry.github.io/r-novice-inflammation/data/r-novice-inflammation-data.zip", destfile = "r-novice-inflammation-data")

#check to make sure it is there
list.files()

#the downloaded file is zipped, so unzip it
unzip(zipfile = "r-novice-inflammation-data")

#check to make sure that the files have downloaded
list.files("data")

#make sure the first one can be read in
read.csv(file = "data/inflammation-01.csv", header = FALSE)

#...start the lesson here

## typo ##
# Note that the X is prepended just a number would not be a valid variable name. 

weight_kg <- 55
weight_kg <- 57.5
# weight in kilograms is now
weight_kg
weight_lb <- 2.2 * weight_kg

dat <- read.csv(file = "data/inflammation-01.csv", header = FALSE)
head(dat)
(c(1,2,3)) %*% t(c(2,3,4))
1 %in% 1:10
mean(c("a", "b"))

apply(X = dat[1:5,], MARGIN = 1, mean)
apply(X = dat[,1:10], MARGIN = 2, mean)
apply(X = dat[,(1:(length(dat[1,])/2)*2)], MARGIN = 2, mean)
