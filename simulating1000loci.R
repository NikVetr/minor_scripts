nLoci <- 1000
indiv <- cbind(rbinom(n = nLoci, size = 1, prob = 0.9), rbinom(n = nLoci, size = 1, prob = 0.5))
nIndiv <- 100
indivs <- lapply(1:nIndiv, function(x) indiv <- cbind(rbinom(n = 1000, size = 1, prob = 0.5), rbinom(n = 1000, size = 1, prob = 0.5)))
popMean <- mean(sapply(1:length(indivs), function(x) sum(indivs[[x]])))
popMean
nGen <- 100
popMeans <- rep(0, nGen)
popMeans[1] <- popMean
for(i in 2:nGen){
  indivs <- lapply(1:nIndiv, function(x) cbind(indivs[[sample(1:nIndiv, size = 1)]][,sample(1:2,1)],indivs[[sample(1:nIndiv, size = 1)]][,sample(1:2,1)]))
  popMean <- mean(sapply(1:length(indivs), function(x) sum(indivs[[x]])))
  popMeans[i] <- popMean
  if(i %% 100 == 0) {plot(popMeans[1:i], type = "l")}
}
plot(popMeans[1:i], type = "l")

# indiv <- cbind(rbinom(n = nLoci, size = 1, prob = 0.999), rbinom(n = nLoci, size = 1, prob = 0.5))
# popMean <- mean(sapply(1:length(indivs), function(x) sum(indivs[[x]])))
# popMeansDisps <- rep(0, nGen)
nGen <- 10000
for(i in 1:nGen){
  indivsDisp <- lapply(1:nIndiv, function(x) cbind(indivs[[sample(1:nIndiv, size = 1)]][,sample(1:2,1)],indivs[[sample(1:nIndiv, size = 1)]][,sample(1:2,1)]))
  popMeansDisps[i] <- mean(sapply(1:length(indivsDisp), function(x) sum(indivsDisp[[x]]))) - popMean
  if(i %% 100 == 0) {print(i);hist(popMeansDisps[1:i], breaks = 50)}
}
mean(popMeansDisps)
ks.test(popMeansDisps, y = pnorm)
