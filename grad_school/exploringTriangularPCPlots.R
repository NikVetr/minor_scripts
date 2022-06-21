#generate an n x n within-group covariance matrix and rescale variances to be less... variable; n is the # of traits
n <- 10
cov <- matrix(rnorm(n^2), n, n); cov <- cov %*% t(cov)
cor <- cov2cor(cov)
sds <- runif(n, 1.5, 2)  
cov <- diag(sds) %*% cor %*% diag(sds)

#sample three populations from multivariate normals whose means are the vertices of a triangle in n-dimensional space
pop1size <- 100
pop1mean <- runif(n = n, 0, 100)
pop1 <- t(chol(cov)) %*% matrix(rnorm(pop1size*n), n, pop1size) + pop1mean
pop2mean <- runif(n = n, 0, 100)
pop2size <- 150
pop2 <- t(chol(cov)) %*% matrix(rnorm(pop2size*n), n, pop2size) + pop2mean
pop3size <- 125
pop3mean <- runif(n = n, 0, 100)
pop3 <- t(chol(cov)) %*% matrix(rnorm(pop3size*n), n, pop3size) + pop3mean

#sample additional populations whose means lie along an edge of the triangle formed above
npops <- 50
smallPopSize <- 30
pop4means <- sapply(1:10, function(x) seq(pop1mean[x], pop3mean[x], length.out = npops))
pop4s <- lapply(1:npops, function(x) t(chol(cov)) %*% matrix(rnorm(smallPopSize*n), n, smallPopSize) + pop4means[x,])
#concatenate
pop4 <- pop4s[[1]]
for(i in 2:npops){
  pop4 <- cbind(pop4, pop4s[[i]])
}

#do PCA on all samples concatenated
alldata <- cbind(pop1, pop2, pop3, pop4)
allCov <- cov(t(alldata))
V <- eigen(allCov)$vectors
L <- eigen(allCov)$values
VarExpl <- round(L / sum(L), 4) * 100
PC1_Scores <- t(V[,1]) %*% alldata / L[1]^0.5 - mean(t(V[,1]) %*% alldata / L[1]^0.5)
PC2_Scores <- t(V[,2]) %*% alldata / L[2]^0.5 - mean( t(V[,2]) %*% alldata / L[2]^0.5 )
#note -- mean centering just makes the average PC-score c(0,0,...); does not affect shape of curve
palette <- colorRampPalette(colors=c("#A50F15", "#006D2C"))
cols <- c(rep("#A50F15", pop1size), rep("#08519C", pop2size), rep("#006D2C", pop3size), rep(palette(npops), each = smallPopSize))

#plot first two PCs
par(mfrow = c(2,1))
plot(PC1_Scores, PC2_Scores, xlab = paste0("PC1 Score, ", VarExpl[1], "% Residual Variance Reduced"),
     ylab = paste0("PC2 Score, ", VarExpl[2], "% Residual Variance Reduced"), col = cols,
     main = "Large Sample of Vertex Pops")

#subsample populations 1-3
samplingFractions <- runif(3, 0.1, 0.5)
someData <- cbind(pop1[,sample(1:pop1size, round(samplingFractions[1]*pop1size), replace = F)], 
                  pop2[,sample(1:pop2size, round(samplingFractions[2]*pop2size), replace = F)], 
                  pop3[,sample(1:pop3size, round(samplingFractions[3]*pop3size), replace = F)], 
                  pop4)
allCov <- cov(t(someData))
V <- eigen(allCov)$vectors
L <- eigen(allCov)$values
VarExpl <- round(L / sum(L), 4) * 100
PC1_Scores <- t(V[,1]) %*% someData / L[1]^0.5 - mean(t(V[,1]) %*% someData / L[1]^0.5)
PC2_Scores <- t(V[,2]) %*% someData / L[2]^0.5 - mean(t(V[,2]) %*% someData / L[2]^0.5)
library(grDevices)
palette <- colorRampPalette(colors=c("#A50F15", "#006D2C"))
cols <- c(rep("#A50F15", round(samplingFractions[1]*pop1size)), rep("#08519C", round(samplingFractions[2]*pop2size)), rep("#006D2C", round(samplingFractions[3]*pop3size)), rep(palette(npops), each = smallPopSize))
plot(PC1_Scores, PC2_Scores, xlab = paste0("PC1 Score, ", VarExpl[1], "% Residual Variance Reduced"), 
     ylab = paste0("PC2 Score, ", VarExpl[2], "% Residual Variance Reduced"), col = cols,
     main = "Small Sample of Vertex Pops")
