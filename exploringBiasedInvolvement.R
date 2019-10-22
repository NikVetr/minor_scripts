#quick refresh on sojourns of a poisson process
# niter <- 1e5
# rate <- 0.01
# time <- rep(0, niter)
# for(i in 1:niter){
#   if(i %% 1e3 == 0){print(i)}
#   x <- rexp(1e4, .1)
#   interval <- 0.01
#   d <- which(x < interval)[1]
#   time[i] <- interval*d+x[d]
# }
# nbreaks = 50
# hist(rexp(niter, .1), breaks = nbreaks, col = rgb(0,1,0,0.3), xlim = c(0,80))
# hist(time, breaks = nbreaks, add = T, col = rgb(1,0,0,0.3))
#yep these are equivalent

#simulate discrete-time, bounded Brownian motion of "relative recruitment source" over 2,000 days
ndays <- 2000
nrecs <- 10 #number of recruitment sources
foo <- sapply(1:nrecs, function(x) abs(cumsum(rnorm(n = ndays, mean = rnorm(1,0,0.025), sd = 5))))
foo <- t(t(foo) + sample(x = 10:200, size = nrecs, replace = T)) #add variable starting weights
weights <- t(sapply(1:length(foo[,1]), function(x) foo[x,]/sum(foo[x,])))
#visualize relative recruitment source weights through time
plot(weights[,1], type ="l",xlab="day", ylab=paste0("proportion of recruitment from each source, # sources = ", nrecs), 
                 xlim=c(1, ndays), ylim=c(0, 0.3));for(i in 2:nrecs){lines(weights[,i], col = RColorBrewer::brewer.pal(10, "Set3")[i])}

#simulate absolute and and source-specific recruiting effort through time with trended, bounded BM
start <- 100 #starting recruitment effort
nindiv <- abs(50 + round(cumsum(rnorm(ndays, mean = 0.1, sd = 5)))) #number of individuals recruited each day, where the first day recruits 50
recruits <- t(sapply(1:ndays, function(x) rmultinom(n = 1, size = nindiv[x], prob = weights[x,])))
expLatentEAness <- sample(5:30, size = nrecs, replace = T) #sample a latent "EAness" variable for each recruitment source that will influence expected source-specific involvement and source-specific extinction rate
expInvolvementPerSrc <- expLatentEAness*rnorm(nrecs, mean = 1, sd = 0.05) + rnorm(nrecs, 0, sd = 0.2)
dailyExpExtinctionRatePerSrc <- MCMCpack::rdirichlet(n=1, alpha = 100/expLatentEAness)/20

allExtantEAs <- matrix(0, ncol = 3, nrow = 1E6) #initialize large matrix to avoid rbind/cbind overhead slowdown
colnames(allExtantEAs) <- c("recruitmentSource", "tenureLength", "involvement")
nSurvivorsTotal <- 0 #initialize starting number of survivors (day 0)
for(i in 1:ndays){
  if(i%%100==0){print(paste0("day ", i))}
  for(j in 1:nrecs){
    day <- i
    source <- j
    if(recruits[day,source] != 0){ #if new members were recruited that day
      RecruitsNoise <- abs(rnorm(recruits[day,source], 1, 0.1)) #recruit-specific noise that influences both involvement and prob(extinction)
      RecruitsInvolvement <- RecruitsNoise * rnorm(recruits[day,source], mean = expInvolvementPerSrc[source], sd = 3) #simulate involvement metrics for each recruit
      RecruitsExtinctionRate <- (1/RecruitsNoise) * dailyExpExtinctionRatePerSrc[source] #simulate extinction rates for each recruit 
      extinctionSojourn <- rexp(recruits[day,source], RecruitsExtinctionRate) #model extinction of each recruit as a poisson process, with exponentially distributed waiting times between events (i.e. until extinction)
      survivalToPresent <- extinctionSojourn > (ndays - day) #ask if each of the given recruits survived to be sampled in the present
      if(any(survivalToPresent)){
        nSurvivors <- sum(survivalToPresent)
        survivingInvolvementStart <- RecruitsInvolvement[survivalToPresent]
        involvementLength <- ndays - day
        survivingInvolvement <- survivingInvolvementStart + rnorm(n = nSurvivors, mean = involvementLength/200, sd = involvementLength/200) #let's say those who remain also get more involved, typically
        allExtantEAs[((nSurvivorsTotal+1):(nSurvivorsTotal+nSurvivors)),] <- cbind(source, involvementLength, survivingInvolvement)
        nSurvivorsTotal <- nSurvivorsTotal+nSurvivors #update the number of survivors (for indexing purposes)
      }
    }
  }
}
allExtantEAs <- allExtantEAs[allExtantEAs[,1] != 0,] #trim all the zeros
#probability of responding to survey proportional to involvement; survey samples 5000 individuals
probSurveyed <- allExtantEAs[,3]/sum(allExtantEAs[,3]); probSurveyed[probSurveyed < 0] <- 0; probSurveyed <- probSurveyed / sum(probSurveyed) 
surveySample <- sample(x = 1:length(allExtantEAs[,3]), size = 5000, replace = F, prob = probSurveyed)
surveySample <- as.data.frame(allExtantEAs[surveySample,])
library(ggplot2); library(ggthemes) #do some quick visualization
surveySample$recruitmentSource <- as.factor(surveySample$recruitmentSource)
ggplot(surveySample, aes(x=recruitmentSource, y=involvement)) + geom_violin() + theme_minimal()
plot(surveySample$recruitmentSource, surveySample$involvement, xlab = "recruitment source", ylab = "involvement")
cor(surveySample[,2], surveySample[,3])
plot(surveySample$tenureLength, surveySample$involvement, xlab = "length (in days) of involvement", ylab = "involvement", col = rgb(0,0,0,0.5))

#fit three basic linear models
library(rethinking)
#ok, now let's have m1, but with patient-specific and protein-specific effects, but have them be deviations from some mean effect (i.e. an intercept)
d <- surveySample
d$recruitmentSource <- as.integer(d$recruitmentSource)
d$tenureLengthStd <- (d$tenureLength - mean(d$tenureLength)) / sd(d$tenureLength)
d$involvementStd <- (d$involvement - mean(d$involvement)) / sd(d$involvement)

m0 <- map2stan(
  alist(
    involvementStd ~ dnorm(mu, sigma),
    mu <- a + c[recruitmentSource], 
    a ~ dnorm(0,1),
    c[recruitmentSource] ~ dnorm(0, sigC),
    c(sigC, sigma) ~ dcauchy(0,1)
  ) ,
  data= d,
  iter = 3000, warmup = 3000, chains = 1)
precis(m0, depth=2)
m0s <- extract.samples(m0)
plot(expInvolvementPerSrc, apply(m0s$c, 2, mean), main = "no time predictor",
     xlab = "true source-specific involvement effect", ylab = "inferred standardized source-specific involvement effect (posterior mean)")

#add linear predictor for time
m1 <- map2stan(
  alist(
    involvementStd ~ dnorm(mu, sigma),
    mu <- a + b*tenureLengthStd + c[recruitmentSource], 
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    c[recruitmentSource] ~ dnorm(0, sigC),
    c(sigC, sigma) ~ dcauchy(0,1)
  ) ,
  data= d,
  iter = 3000, warmup = 3000, chains = 1)
precis(m1, depth=2)
m1s <- extract.samples(m1)
plot(expInvolvementPerSrc, apply(m1s$c, 2, mean), main = "time predictor",
     xlab = "true source-specific involvement effect", ylab = "inferred standardized source-specific involvement effect (posterior mean)")
# postcheck(m1, window = 500)

m2 <- map2stan(
  alist(
    involvementStd ~ dnorm(mu, sigma),
    mu <- a + b*tenureLengthStd + c[recruitmentSource], 
    a ~ dnorm(0,1),
    b <- b_base + b_interaction[recruitmentSource],
    b_base ~ dnorm(0,1),
    b_interaction[recruitmentSource] ~ dnorm(0,sigB),
    c[recruitmentSource] ~ dnorm(0, sigC),
    c(sigB, sigC, sigma) ~ dcauchy(0,1)
  ) ,
  data= d,
  iter = 3000, warmup = 3000, chains = 1)
precis(m2, depth=2)
m2s <- extract.samples(m2)
plot(expInvolvementPerSrc, apply(m2s$c, 2, mean), xlab = "true source-specific involvement effect", main = "time predictor (with interactions)",
     ylab = "inferred standardized source-specific involvement effect (posterior mean)")

compare(m0,m1,m2)
