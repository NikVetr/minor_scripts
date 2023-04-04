library(magick)

mat2img <- function(x){
  image_read(1-abind::abind(x, x, x, along = 3))
}

x <- matrix(F, 100, 150)
x[4:10, 4:10] <- T
x[84:90, 4:10] <- T
x[84:90, 94:100] <- T
x[0:10, 94:100] <- T


x[50:60, 0:6] <- T
x[50:60, 95:100] <- T

x[cbind(sample(1:nrow(x), 5000, T), sample(1:ncol(x), 5000, T))] <- T

ximg <- mat2img(x)
inds <- which(x, T)

xr <- c(1, ncol(x))
yr <- c(1, nrow(x))
plot(inds, xlim = xr, ylim = yr)
abline(h = yr)
abline(v = xr)

dimk <- 2
ckern <- matrix(1,dimk,dimk)
ximg2 <- image_convolve((ximg), kernel = ckern, bias = '50%')
foo <- image_data(ximg2)
foo <- Reduce("+", lapply(1:3, function(i) matrix(as.integer(foo[i,,]), nrow(foo[i,,]), ncol(foo[i,,]))))
foo <- t(scale2(abs(foo - median(foo))) >= 0.5)
foo_inds <- which(foo, T)
points(foo_inds, pch = 19)

#now try deconvolving
dimkr <- c(ceiling(-dimk/2), floor(dimk/2)) - c(0, dimk %% 2 == 0)
dimkr_combos <- expand.grid(dimkr[1]:dimkr[2], dimkr[1]:dimkr[2])
truevals <- do.call(rbind, lapply(1:nrow(dimkr_combos), function(dkri) t(t(foo_inds) + unlist(dimkr_combos[dkri,]))))
truevals <- truevals[truevals[,1] >= yr[1] &
                       truevals[,1] <= yr[2] &
                       truevals[,2] >= xr[1] &
                       truevals[,2] <= xr[2], ]
truevals <- truevals[!duplicated(truevals),]
xdc <- foo
xdc[truevals] <- T
foo_inds_dc <- which(xdc, T)
points(foo_inds_dc, pch = 19, col = adjustcolor(2, 0.2))

#now put it in a function
deconvolve_1kern <- function(boolmat, dimk){
  boolmat_inds <- which(boolmat, T)
  dimkr <- c(ceiling(-dimk/2), floor(dimk/2)) - c(0, dimk %% 2 == 0) #our window around each convolved point
  dimkr_combos <- expand.grid(dimkr[1]:dimkr[2], dimkr[1]:dimkr[2]) #vals to add
  truevals <- do.call(rbind, lapply(1:nrow(dimkr_combos), function(dkri) 
    t(t(boolmat_inds) + unlist(dimkr_combos[dkri,])))) #do the addition
  
  #filter out impossible or duplicated values
  xr <- c(1, ncol(boolmat))
  yr <- c(1, nrow(boolmat))
  truevals <- truevals[truevals[,1] >= yr[1] &
                         truevals[,1] <= yr[2] &
                         truevals[,2] >= xr[1] &
                         truevals[,2] <= xr[2], ]
  
  #this is costlier than just doing the redundant assignment
  # truevals <- truevals[!duplicated(truevals),]
  
  #create and return new matrix of bools
  xdc <- boolmat
  xdc[truevals] <- T
  xdc
}

points(which(deconvolve_1kern(foo, dimk), T), col = adjustcolor(2, 0.2), pch = 19)
