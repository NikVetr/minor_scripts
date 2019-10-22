# Maximize Quality*Time With Partner
ageAtFirstSample <- 12 #when we first start evaluating potential matches
ageAtDeath <- 80 #when we die
evalDaysPerPerson <- 30 #how many days we sample each partner for
ageSeq <- seq(from = ageAtFirstSample, to = ageAtDeath, by = 1) #range of ages across which we can begin committing
noCommitmentReductionFactor <-.8 #how good casual dating is compared to a committed relationship
ageWisdomBonus <- 0 #how much to shift the mean of our distribution at age 80 relative to age 12?
perPartnerWisdomBonus <- ageWisdomBonus / ((ageAtDeath - ageAtFirstSample) * 365 / evalDaysPerPerson)
perceptionUncertaintySDMax <- 0
perPartnerUncertaintySD <- perceptionUncertaintySDMax / ((ageAtDeath - ageAtFirstSample) * 365 / evalDaysPerPerson)

CDFval <- vector(length = length(ageSeq))
samples <- 1e3 #how many simulations from each age do we want to run?
EQT <- vector(length = length(ageSeq))
ageAtCommitment <- vector(length = length(ageSeq))
for(j in 2:length(ageSeq)){ #make a for loop took look at age at which you might begin committing
  jEQT <- vector(length = samples) #initialize an empty vector to hold data from our simulations
  jAAC <- vector(length = samples) #initialize an empty vector to hold data from our simulations
  jCDF <- vector(length = samples) #initialize an empty vector to hold data from our simulations
  for(i in 1:samples){
    if (i %% 50 == 0) {print(c(j, i))} #print a counter to tell us where in the sims we are
    search <- rnorm((ageSeq[j] - ageAtFirstSample) * 365 / evalDaysPerPerson, mean = 5, sd = 1) #simulate partners during initial sampling
    search <- search + seq(length.out=length(search), from = 0, by = perPartnerWisdomBonus)
#     errorSearch <- sapply(seq(from = perceptionUncertaintySDMax, by = -(perPartnerUncertaintySD), length.out = length(search)), 
#                           function(x) rnorm(1,mean=0,sd=x))
#     searchF <- search + errorSearch
    life <- rnorm((ageAtDeath - ageSeq[j]) * 365 / evalDaysPerPerson, mean = 5, sd = 1) #simulate partners during the rest of our lives
    life <- life + seq(length.out = length(life), from = length(search) * perPartnerWisdomBonus, by = perPartnerWisdomBonus)
#     errorLife <- sapply(seq(from = perceptionUncertaintySDMax - length(search) * perPartnerUncertaintySD, by = -(perPartnerUncertaintySD), 
#                             length.out = length(life)), function(x) rnorm(1,mean=0,sd=x))
#     lifeF <- life + errorLife
    theOne <- which(life > max(search)) #find which partners during the rest of our lives are 
    #better than the best partner from initial sampling. Change to lifeF and searchF when accomodating uncertainty
    if(!is.na(theOne[1])) #if there aren't any better partners during the rest of our lives, we keep dating forever
      { jEQT[i] <- (length(life) - theOne[1]) * evalDaysPerPerson * life[theOne[1]] + #Quality Time with our chosen partner
      (sum(search) + sum(life[1:(theOne[1]-1)])) * evalDaysPerPerson * noCommitmentReductionFactor;#Quality Time with our casual dates
      jAAC[i] <- ageSeq[j] + theOne[1] * 30 / 365
      jCDF[i] <- pnorm(life[theOne[1]], mean=5, sd=1)
       } else 
      {jEQT[i] <- (sum(search) + sum(life)) * evalDaysPerPerson * noCommitmentReductionFactor} #Quality Time with our lifetime of dating if we never commit
    
  }
  EQT[j] <- sum(jEQT)/samples #find the average quality time across simulations for each age of begin commitment
  ageAtCommitment[j] <- mean(setdiff(jAAC, 0))
  CDFval[j] <- mean(setdiff(jCDF, 0))
}
EQT
plot(ageSeq[-1], EQT[-1], xlab="ageSeq", ylab="EQT")
plot(ageSeq[3:length(ageSeq)-1], EQT[3:length(ageSeq)-1], xlab="ageSeq", ylab="EQT",ylim=c(100000,170000))
ageSeq[(which.max(EQT))]
plot(ageSeq, ageAtCommitment)
abline(a=0, b=1)


samples <- 1e3 #how many simulations from each age do we want to run?
EQT <- vector(length = length(ageSeq))
ageAtCommitment <- vector(length = length(ageSeq))
for(j in 1:length(ageSeq)){ #make a for loop took look at age at which you might begin committing
  jEQT <- vector(length = samples) #initialize an empty vector to hold data from our simulations
  jAAC <- vector(length = samples) #initialize an empty vector to hold data from our simulations
  for(i in 1:samples){
    if (i %% 50 == 0) {print(c(j, i))} #print a counter to tell us where in the sims we are
    search <- rnorm((ageSeq[j] - ageAtFirstSample) * 365 / evalDaysPerPerson, mean = 5, sd = 1) #simulate partners during initial sampling
    life <- rnorm((ageAtDeath - ageSeq[j]) * 365 / evalDaysPerPerson, mean = 5, sd = 1) #simulate partners during the rest of our lives
    theOne <- which(life > max(search)) #find which partners during the rest of our lives are better than the best partner from initial sampling
    if(!is.na(theOne[1])) #if there aren't any better partners during the rest of our lives, we keep dating forever
    { jEQT[i] <- (length(life) - theOne[1]) * evalDaysPerPerson * life[theOne[1]] + #Quality Time with our chosen partner
      (sum(search) + sum(life[1:(theOne[1]-1)])) * evalDaysPerPerson * noCommitmentReductionFactor;#Quality Time with our casual dates
    jAAC[i] <- ageSeq[j] + theOne[1] * 30 / 365
    jCDF[i] <- pnorm(life[theOne[1]], mean=5, sd=1)
    } else 
    {jEQT[i] <- (sum(search) + sum(life)) * evalDaysPerPerson * noCommitmentReductionFactor} #Quality Time with our lifetime of dating if we never commit
    
  }
  EQT[j] <- sum(jEQT)/samples #find the average quality time across simulations for each age of begin commitment
  ageAtCommitment[j] <- mean(setdiff(jAAC, 0))
}


betterCDF <- vector(length = samples)
wait <- vector(length = samples)
for(i in 1:samples){
  if (i %% 50 == 0) {print(i)}
  search <- rnorm((ageAtBeginCommitment - ageAtFirstSample) * 365 / evalDaysPerPerson, mean = 5, sd = 1)
  life <- rnorm((ageAtDeath - ageAtBeginCommitment) * 365 / evalDaysPerPerson, mean = 5, sd = 1)
  theOne <- which(life > max(search))
  betterCDF[i] <- pnorm(life[theOne[1]], mean=5, sd = 1)
  wait[i] <- (theOne[2] - theOne[1]) * evalDaysPerPerson / 365
}

##testing multivariate -> univariate
quantile(betterCDF, c(.1, .9), na.rm = T)
quantile(wait, c(.1, .9), na.rm = T)
median(wait, na.rm=T)


sig <- matrix(c(1,2,3,
                2,5,7,
                3,7,15), ncol=3)

sig <- matrix(c(10,2,3,1,4,
                2,5,7,1,2,
                3,7,15,3,1,
                1,1,3,11,3,
                4,2,1,3,14), ncol=5)

test <- rmvnorm(n=1e5, sigma = sig)
test2 <- sqrt(test[,1]^2 + test[,2]^2 + test[,3]^2 + test[,4]^2 + test[,5]^2)
dens((test2))

#uniform distribution?

ageAtFirstSample <- 12
ageAtBeginCommitment <- 20
ageAtDeath <- 80
evalDaysPerPerson <- 1000
ageSeq <- seq(from = ageAtFirstSample, to = ageAtDeath, by = 1)


samples <- 1e3
EQT <- vector(length = length(ageSeq))
for(j in 1:length(ageSeq)){
  jEQT <- vector(length = samples)
  for(i in 1:samples){
    if (i %% 50 == 0) {print(c(j, i))}
    search <- runif((ageSeq[j] - ageAtFirstSample) * 365 / evalDaysPerPerson, min = 0, max = 10)
    life <- runif((ageAtDeath - ageSeq[j]) * 365 / evalDaysPerPerson, min = 0, max = 10)
    theOne <- which(life > max(search))
    jEQT[i] <- (length(life) - theOne[1]) * evalDaysPerPerson * life[theOne[1]] 
  }
  jEQT[is.na(jEQT)] <- 0
  EQT[j] <- sum(jEQT)/samples
}
EQT
plot(ageSeq, EQT)