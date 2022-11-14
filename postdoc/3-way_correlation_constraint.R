library(corpcor)
ifelse2 <- function(bool, opt1, opt2){if(bool){return(opt1)}else{return(opt2)}}
gradleg <- function(loc, cols, labels, pos = 4, main = "", rasterize = F, border_in = NA, border_out = T, direc = "v", ...){
  
  #get metadata
  n <- length(cols)
  
  #get rect positions
  if(direc == "v"){
    locs <- data.frame(x1 = rep(loc[1], n),
                       x2 = rep(loc[2], n),
                       y1 = seq(loc[3], loc[4], length.out = n+1)[-(n+1)],
                       y2 = seq(loc[3], loc[4], length.out = n+1)[-1]
    )  
  } else {
    locs <- data.frame(x1 = seq(loc[1], loc[2], length.out = n+1)[-(length(n)+1)],
                       x2 = seq(loc[1], loc[2], length.out = n+1)[-1],
                       y1 = rep(loc[3], n),
                       y2 = rep(loc[4], n)
    )
  }
  
  #draw rects
  rect(locs$x1, locs$y1, locs$x2, locs$y2,
       col = cols, border = border_in)
  if(border_out){
    rect(loc[1], loc[3], loc[2], loc[4])
  }
  
  #draw text
  nl <- length(labels)
  if(nl != n){
    if(direc == "v"){
      locs <- data.frame(x1 = rep(loc[1], nl),
                         x2 = rep(loc[2], nl),
                         y1 = seq(loc[3], loc[4], length.out = nl+1)[-(nl+1)],
                         y2 = seq(loc[3], loc[4], length.out = nl+1)[-1]
      )  
    } else {
      locs <- data.frame(x1 = seq(loc[1], loc[2], length.out = nl+1)[-(length(nl)+1)],
                         x2 = seq(loc[1], loc[2], length.out = nl+1)[-1],
                         y1 = rep(loc[3], nl),
                         y2 = rep(loc[4], nl)
      )
    }
  }
  if(direc == "v") text(x = ifelse2(pos == 4, locs$x2, locs$x1) + 0.03, y = (locs$y1 + locs$y2) / 2, pos = pos, labels = labels)
  if(direc == "h") text(y = ifelse2(pos == 1, locs$x1, locs$x2), x = (locs$x1 + locs$x2) / 2, pos = pos, labels = labels)
  
  #draw title
  text(x = mean(loc[1:2]), y = loc[4], labels = main, pos = 3)
  
}

cmat_3 <- function(a, b, c){
  R <- diag(3)
  R[1,2] <- R[2,1] <- a
  R[1,3] <- R[3,1] <- b
  R[2,3] <- R[3,2] <- c
  R
}

check_psd <- function(R) all(eigen(R)$values > -1E-6)

try_bs <- function(a, b, c, n = 0, max_n = 20, tol = 1E-3){
  if(n <= max_n){
    works <- check_psd(cmat_3(a, b, c))
    new_n <- n + !works
    new_b <- b - 1/(2^n) * ifelse(works, 1, -1/2)
    if(abs(b-new_b) < tol){
      new_n <- max_n + 1
    }
    best_b <- try_bs(a, new_b, c, n = new_n, max_n = max_n)
  } else{
    best_b <- b - 1/(2^(n-1)) * ifelse(works, 1, -1)
  }
  return(best_b)
}
try_bs(0.5, 0.4, 0.6)

eigen(cmat_3(0.5, try_bs(0.5, 0.4, 0.6), 0.6))

try_bs(0.5, 0.4, 0.6)
0.5*0.6 - sqrt(1-0.5^2)*sqrt(1-0.6^2)

a <- 0.5
b <- 0.2
eigen(cmat_3(a, a*b + sqrt(1-a^2)*sqrt(1-b^2), b))
eigen(cmat_3(a, a*b - sqrt(1-a^2)*sqrt(1-b^2), b))

#whoops lol no numerical approx necessary!



as <- cs <- setNames(1:10/10, 1:10/10)
bs <- -10:10/10

as <- cs <- setNames(0:200/200, 0:200/200)
da <- dc <- as.numeric(diff(as)[1])
bs <- -100:100/100

min_bs <- parallel::mclapply(as, function(ai){
  sapply(cs, function(ci){
    ai*ci - sqrt(1-ai^2)*sqrt(1-ci^2)
  })  
}, mc.cores = 12)

min_bs <- do.call(rbind, min_bs)
cols <- viridis::viridis(length(bs))
par(mar = c(5,5,2,4))
plot(as, cs, col = "white",
     xlab = latex2exp::TeX("correlation between X$_1$ & X$_2$"),
     ylab = latex2exp::TeX("correlation between X$_2$ & X$_3$"))
for(ai in seq_along(as)){
  for(ci in seq_along(cs)){
    rect(xleft = as[ai] - da/2,
         xright = as[ai] + da/2,
         ybottom = cs[ci] - dc/2,
         ytop =  cs[ci] + dc/2,
         col = cols[which.min(abs(min_bs[ai, ci] - bs))],
         border = NA)
  }
}
par(xpd = NA)
gradleg(loc = c(1.125, 1.075, 0.2, 1), labels = -5:5/5, cols = cols, main = latex2exp::TeX("r$_{1,3}$"))
text(x = 1.10, y = 1.1, labels = "min.")
