library(plotrix)
library(mvtnorm)

# make animation panes
cumu <- c(0,0,0,0)
nrep <- 10000
for(i in 1:nrep){
  print(i)
  # sample center of projectile, with sigma reflecting accuracy @ hitting the origin, currently symmetrically distributed according to two univariate normals
  centerProjectile <- rmvnorm(1, sigma = diag(c(5,5)))
  
  # sample orientation of projectile
  # popular trick to sample uniformly from surface of m-dim hypersphere by rescaling two univ iid standard normals by vector length
  # when m = 2 sphere is a circle lol
  # sample non-uniformly by drawing from different distributions
  widthProjectile <- 8
  unscaledRadius <- c(rnorm(n = 1, mean = 0, sd = 2), rnorm(n = 1, mean = 0, sd = 1))
  scaledRadius <- unscaledRadius / (sum(unscaledRadius^2)^0.5) * widthProjectile / 2
  discretizeProjectileLength <- seq(from = 0, to = 2, by = 0.01)
  hypotheticalProjectileCoordinates <- sapply(1:length(discretizeProjectileLength), function (x) scaledRadius + -scaledRadius * discretizeProjectileLength[x])
  projectileCoordinates <- sapply(1:length(discretizeProjectileLength), function (x) hypotheticalProjectileCoordinates[,x] + t(centerProjectile))
  coordinatesDistanceFromCenter <- sapply(1:length(discretizeProjectileLength), function (x) sum(projectileCoordinates[,x]^2)^.5)
  
  # define radii of target circles
  firstRadius <- 1
  secondRadius <- 2
  thirdRadius <- 4
  
  hitFirstCircle <- coordinatesDistanceFromCenter < firstRadius
  hitSecondCircle <- firstRadius < coordinatesDistanceFromCenter & coordinatesDistanceFromCenter < secondRadius
  hitThirdCircle <- secondRadius < coordinatesDistanceFromCenter & coordinatesDistanceFromCenter < thirdRadius
  hitNothing <- coordinatesDistanceFromCenter > thirdRadius
  
  discountinuity2 <- diff(which(hitSecondCircle))
  discountinuity3 <- diff(which(hitThirdCircle))
  discountinuity4 <- diff(which(hitNothing))
  
  hitInner <- sum(any(hitFirstCircle))
  hitMiddle <- sum(any(hitSecondCircle))
  hitOuter <- sum(any(hitThirdCircle))
  missed <- sum(any(hitNothing))
  
  allHits <- c(hitInner, hitMiddle, hitOuter, missed)
  cumu[min(which(allHits == 1))] <- cumu[min(which(allHits == 1))] + 1
  
  if(i < 50 | (i < 500 & i %% 10 == 0) | (i < 5000 & i %% 100 == 0) | (i < 50000 & i %% 1000 == 0)){
    png(paste0("C:\\Users\\Nikolai\\Desktop\\Dartboard\\Width", widthProjectile, "Rep", i, ".png"), width = 600, height = 600)
    plot(0,type='n',axes=F,ann=FALSE, frame.plot = T, xlim = c(-5, 5), ylim = c(-5, 5), asp = 1)
    draw.circle(x = 0, y = 0, radius = c(firstRadius, secondRadius, thirdRadius), col = c(rgb(0,0,0,.1),rgb(0,0,0,.05),rgb(0,0,0,.025)), lwd = 2)
    
    title(main = paste0("Throwing Projectiles at a Dartboard Simulation\n", "Projectile Width = ", widthProjectile, ", n = ", i))
    
    # firstCircle
    lines(x = projectileCoordinates[1,hitFirstCircle], y = projectileCoordinates[2,hitFirstCircle], lwd = 5, col = rgb(1,0,0,1))
    # secondCircle
    if(!any(discountinuity2 > 1)){
      lines(x = projectileCoordinates[1,hitSecondCircle], y = projectileCoordinates[2,hitSecondCircle], lwd = 4, col = rgb(1,0,0,.66))
    } else {
      lines(x = projectileCoordinates[1,which(hitSecondCircle)[1:which.max(discountinuity2)]], y = projectileCoordinates[2,which(hitSecondCircle)[1:which.max(discountinuity2)]], lwd = 4, col = rgb(1,0,0,.66))
      lines(x = projectileCoordinates[1,which(hitSecondCircle)[(which.max(discountinuity2)+1):length(which(hitSecondCircle))]], y = projectileCoordinates[2,which(hitSecondCircle)[(which.max(discountinuity2)+1):length(which(hitSecondCircle))]], lwd = 4, col = rgb(1,0,0,.66))
    }
    # thirdCircle
    if(!any(discountinuity3 > 1)){
      lines(x = projectileCoordinates[1,hitThirdCircle], y = projectileCoordinates[2,hitThirdCircle], lwd = 4, col = rgb(1,0,0,.33))
    } else {
      lines(x = projectileCoordinates[1,which(hitThirdCircle)[1:which.max(discountinuity3)]], y = projectileCoordinates[2,which(hitThirdCircle)[1:which.max(discountinuity3)]], lwd = 4, col = rgb(1,0,0,.33))
      lines(x = projectileCoordinates[1,which(hitThirdCircle)[(which.max(discountinuity3)+1):length(which(hitThirdCircle))]], y = projectileCoordinates[2,which(hitThirdCircle)[(which.max(discountinuity3)+1):length(which(hitThirdCircle))]], lwd = 4, col = rgb(1,0,0,.33))
    }
    # missed
    if(!any(discountinuity4 > 1)){
      lines(x = projectileCoordinates[1,hitNothing], y = projectileCoordinates[2,hitNothing], lwd = 4, col = rgb(0,1,0,.66))
    } else {
      lines(x = projectileCoordinates[1,which(hitNothing)[1:which.max(discountinuity4)]], y = projectileCoordinates[2,which(hitNothing)[1:which.max(discountinuity4)]], lwd = 4, col = rgb(0,1,0,.66))
      lines(x = projectileCoordinates[1,which(hitNothing)[(which.max(discountinuity4)+1):length(which(hitNothing))]], y = projectileCoordinates[2,which(hitNothing)[(which.max(discountinuity4)+1):length(which(hitNothing))]], lwd = 4, col = rgb(0,1,0,.66))
    }
    
    text(x = 4.5, y = 5.1, labels = "BEST HIT", cex = 1.5, col = rgb(1,0,0,1))
    text(x = 4.5, y = 4.6, labels = paste0("inner = ", round(cumu[1]/i*100, digits = 1), "%"), cex = 1.25)
    text(x = 4.5, y = 4.1, labels = paste0("middle = ", round(cumu[2]/i*100, digits = 1), "%"), cex = 1.25)
    text(x = 4.5, y = 3.6, labels = paste0("outer = ", round(cumu[3]/i*100, digits = 1), "%"), cex = 1.25)
    text(x = 4.5, y = 3.1, labels = paste0("missed = ", round(cumu[4]/i*100, digits = 1), "%"), cex = 1.25)
    
    dev.off()
  }
}



# examine behavior quantitatively

widths <- c(.1,.5,1,2,3,4,8,16)
cumu <- matrix(nrow = length(widths), ncol = 4, data = 0)
nrep <- 15000
for(j in 1:length(widths))
  for(i in 1:nrep){
    print(c(j,i))
    # sample center of projectile, with sigma reflecting accuracy @ hitting the origin, currently symmetrically distributed according to two univariate normals
    centerProjectile <- rmvnorm(1, sigma = diag(c(5,5)))
    
    # sample orientation of projectile
    # popular trick to sample uniformly from surface of m-dim hypersphere by rescaling two univ iid standard normals by vector length
    # when m = 2 sphere is a circle lol
    # sample non-uniformly by drawing from different distributions
    widthProjectile <- widths[j]
    unscaledRadius <- c(rnorm(n = 1, mean = 0, sd = 2), rnorm(n = 1, mean = 0, sd = 1))
    scaledRadius <- unscaledRadius / (sum(unscaledRadius^2)^0.5) * widthProjectile / 2
    discretizeProjectileLength <- seq(from = 0, to = 2, by = 0.01)
    hypotheticalProjectileCoordinates <- sapply(1:length(discretizeProjectileLength), function (x) scaledRadius + -scaledRadius * discretizeProjectileLength[x])
    projectileCoordinates <- sapply(1:length(discretizeProjectileLength), function (x) hypotheticalProjectileCoordinates[,x] + t(centerProjectile))
    coordinatesDistanceFromCenter <- sapply(1:length(discretizeProjectileLength), function (x) sum(projectileCoordinates[,x]^2)^.5)
    
    # define radii of target circles
    firstRadius <- 1
    secondRadius <- 2
    thirdRadius <- 4
    
    hitFirstCircle <- coordinatesDistanceFromCenter < firstRadius
    hitSecondCircle <- firstRadius < coordinatesDistanceFromCenter & coordinatesDistanceFromCenter < secondRadius
    hitThirdCircle <- secondRadius < coordinatesDistanceFromCenter & coordinatesDistanceFromCenter < thirdRadius
    hitNothing <- coordinatesDistanceFromCenter > thirdRadius
    
    discountinuity2 <- diff(which(hitSecondCircle))
    discountinuity3 <- diff(which(hitThirdCircle))
    discountinuity4 <- diff(which(hitNothing))
    
    hitInner <- sum(any(hitFirstCircle))
    hitMiddle <- sum(any(hitSecondCircle))
    hitOuter <- sum(any(hitThirdCircle))
    missed <- sum(any(hitNothing))
    
    allHits <- c(hitInner, hitMiddle, hitOuter, missed)
    cumu[j, min(which(allHits == 1))] <- cumu[j, min(which(allHits == 1))] + 1
  }

cumu <- cumu/nrep*100
plot(0,type='n',axes=T,ann=T, frame.plot = T, xlim = c(0, 16), xlab = "Width of Projectile", ylim = c(0, 50), ylab = "Percent of Time Circle Was the Maximum Circle Hit", main = "Accuracy at Hitting Target with Projectiles of Different Widths")
lines(x = widths, y = cumu[,1], lwd = 4, col = rgb(1,0,0,1))
lines(x = widths, y = cumu[,2], lwd = 2, col = rgb(1,0,0,.66))
lines(x = widths, y = cumu[,3], lwd = 2, col = rgb(1,0,0,.33))
lines(x = widths, y = cumu[,4], lwd = 2, col = rgb(0,1,0,.66))
legend('topright', c("Inner Circle", "Middle Donut", "Outer Donut", "Complete Miss"), col = c(rgb(1,0,0,1), rgb(1,0,0,0.66), rgb(1,0,0,0.33), rgb(0,1,0,.66)), lwd = c(4,2,2,2))
