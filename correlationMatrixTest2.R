library(mvtnorm)
library(MCMCpack)

meandiff <- vector()
for(i in 1:1000){
  print(i)
  simD <- rmvnorm(3,sigma = cov)
  meandiff[i] <- mean(cov(simD) - cov)
}
hist(meandiff)


c <- matrix(c(1,.5,.7,.5,1,.6,.7,.6,1), nrow = 3)


niter <- 1e5
df <- 300
wishes <- list()
for(i in 1:niter){
  print(i)
  wishes[[i]] <- rwish(df, diag(3))
}

wishez <- list()
for(i in 1:niter){
  print(i)
  simdata <- rmvnorm(n = df, sigma = diag(3))
  wishez[[i]] <- cov(simdata) * df
}


var1 <- sapply(1:length(wishes), function (x) wishes[[x]][1,1])
var2 <- sapply(1:length(wishez), function (x) wishez[[x]][1,1])

hist(var1, breaks = 100)
hist(var2, breaks = 100)

corres <- list()
for(i in 1:niter){
  corres[[i]] <- cov2cor(wishes[[i]])
}

var1 <- sapply(1:length(wishes), function (x) wishes[[x]][1,1])
cor1 <- sapply(1:length(corres), function (x) corres[[x]][1,2])

plot(cor1, log(var1))

#Wishart samples
rwishart <- function(r,R)
{
  X <- rmvnorm(r,sig=R)
  t(X)%*%X
}


#create an MCMC sampler for unknown off diag components

c <- matrix(c(1,.5,.7,.5,1,.6,.7,.6,1), nrow = 3)

data <- rmvnorm(n = 100, sigma = c)

sub <- c[1:2, 1:2]
iter <- 1e4
matrices <- list()
matricesWish <- list()
likelihoods <- vector()
posteriors <- vector()
priors <- vector()
proposalTuner <- 10000
matricesWish[[1]] <- rwish(v = proposalTuner, S = diag(3))
matrices[[1]] <- cov2cor(matricesWish[[1]])
subDF <- 10
subP <- log(dwish(W = matrices[[1]][1:2, 1:2], v = subDF, S = sub * subDF))
wholeDF <- 3
cP <- log(dwish(W = matrices[[1]], v = wholeDF, S = wholeDF * diag(3)))
priors[1] <- subP + cP
likelihoods[1] <- sum(dmvnorm(x = data, sigma = matrices[[1]], log = T))
posteriors[1] <- priors[[1]] + likelihoods[[1]]

for(i in 2:iter){
  if(i %% 100 == 0) {print(i)}
  current <- matrices[[i-1]]
  currentWish <- matricesWish[[i-1]]
  proposedWish <- rwish(v = proposalTuner, S = current)
  proposed <- cov2cor(proposedWish)
  proposedPrior <- log(dwish(W = proposed[1:2, 1:2], v = subDF, S = sub * subDF)) + 
    log(dwish(W = proposed, v = wholeDF, S = wholeDF * diag(3)))
  proposedLikelihood <- sum(dmvnorm(x = data, sigma = proposed, log = T))
  proposedPosterior <- proposedLikelihood + proposedPrior
  acceptanceProbLog <- proposedPosterior - posteriors[[i-1]] + 
    log(dwish(W = proposedWish, v = proposalTuner, S = current)) - 
          log(dwish(W = currentWish, v = proposalTuner, S = proposed))
  acceptanceNumber <- log(runif(n = 1, min = 0, max = 1))
  if(acceptanceNumber > acceptanceProbLog){
    posteriors[i] <- proposedPosterior
    matrices[[i]] <- proposed
    matricesWish[[i]] <- proposedWish
    likelihoods[i] <- proposedLikelihood
  } else {
    posteriors[i] <- posteriors[i-1]
    matrices[[i]] <- matrices[[i-1]]
    likelihoods[i] <- likelihoods[i-1]
    matricesWish[[i]] <- matricesWish[[i-1]]
    
  }
}

noBurnLLs <- likelihoods[iter/10:iter]
noBurnLLs <- noBurnLLs[seq(1, length(noBurnLLs), by =100)]
plot(1:length(noBurnLLs), noBurnLLs, type = "l")



#how far apart are the chimp and human covariance matrices?


d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\chimp.csv")
d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\PanSorted.csv")
d <- d.orig
men <- d[d$Sex == "M",]
d <- men
str(d)


head(as.matrix(d))
unique(d$Taxon)[1]
sp <- sapply(1:length(unique(d$Taxon)), function (x) as.matrix((t(d[d$Taxon == unique(d$Taxon)[x],-(1:2)]))))
ntraits <- 27
nspecies <- 4

#get traits
traits <- list()
for(i in 1:ntraits){ 
  print(i)
  trait <- list()
  trait[[1]] <- sp[[1]][i,]
  for (j in 2:nspecies){
    trait[[j]] <- sp[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = ntraits, ncol = ntraits)
for (i in 1:ntraits){
  print(i)
  for(j in 1:ntraits){
    trait1 <- traits[[i]]
    trait2 <- traits[[j]]
    trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
    trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
    trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
    trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
    deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
    sumALL <- 0
    for (k in (1:length(deviationsProducts))){
      sumALL <- sumALL + sum(deviationsProducts[[k]])
    }
    indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
    covariance <- sumALL/(indiv - length(deviationsProducts))
    cov[i,j] <- covariance
  }
}

offDiagCorr <- c(cov2cor(cov))
offDiagCorr <- offDiagCorr[offDiagCorr < 1] 
par(mfrow = c(1,2))
hist(offDiagCorr, main = "")
title("Correlations Between Traits")
hist(diag(cov), main = "")
title("Variances of Traits")

npop <- length(sp)
#compute matrix of Species means
traitsCHIMP <- matrix(nrow=npop, ncol=ntraits)
for(i in 1:npop){
  for(j in 1:ntraits){
    traitsCHIMP[i,j] <- mean(sp[[i]][j,])
  }
}
rownames(traitsCHIMP) <- unique(d$Taxon)
colnames(traitsCHIMP) <- colnames(d)[-(1:2)]

rownames(cov) <- colnames(cov) <- colnames(d)[-(1:2)]
covChimp <- cov


d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\Howell.csv")
d <- d.orig
# str(d)
# d[d == 0] <- NA
# d <- d[complete.cases(d),]
men <- d[d$Sex == "M",]
#men <- d
#head(as.matrix(men))
#unique(men$Population)[1]
pops <- sapply(1:length(unique(men$Population)), function (x) as.matrix((t(men[men$Population == unique(men$Population)[x],-(1:4)]))))

#get traits
traits <- list()
for(i in 1:82){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:30){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 82, ncol = 82)
for (i in 1:82){
  print(i)
  for(j in 1:82){
    trait1 <- traits[[i]]
    trait2 <- traits[[j]]
    trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
    trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
    trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
    trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
    deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
    sumALL <- 0
    for (k in (1:length(deviationsProducts))){
      sumALL <- sumALL + sum(deviationsProducts[[k]])
    }
    indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
    covariance <- sumALL/(indiv - length(deviationsProducts))
    cov[i,j] <- covariance
  }
}

offDiagCorr <- c(cov2cor(cov))
offDiagCorr <- offDiagCorr[offDiagCorr < 1] 
par(mfrow = c(1,2))
hist(offDiagCorr, main = "")
title("Correlations Between Traits")
hist(diag(cov), main = "")
title("Variances of Traits")

npop <- length(pops)
#compute matrix of population means
traitsHOWELL <- matrix(nrow=npop, ncol=82)
for(i in 1:npop){
  for(j in 1:82){
    traitsHOWELL[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELL) <- unique(men$Population)
colnames(traitsHOWELL) <- colnames(men)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(men)[-(1:4)]
covHuman <- cov



populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU", "BERG", "ZALAVAR", "AINU", "NORSE", "EASTER I", "MOKAPU")
#populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU")
numTraits <- 57 #for the linear measurements 
#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHuman <- traitsHOWELL[,linMeasNames]
linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]


##################################################
###### Finding Intersect of Human and Chimp ######
##################################################

chimpHumanTraits <- intersect(rownames(covChimp), rownames(covHuman))
covC <- sharedChimpCov <- covChimp[chimpHumanTraits, chimpHumanTraits]
covH <- sharedHumanCov <- covHuman[chimpHumanTraits, chimpHumanTraits]

# how far is the chimp from the human?

rwish(v = 27, S = covH) / 27
