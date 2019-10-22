earliestHomo <- -2.8
sediba <- -2.0
obsDiff <- sediba - earliestHomo
temporalRange <- 2
nrep <- 1e4

#########################################################################################################
# sample sediba temporal intervals, assuming the observed sediba is a uniform sample from that interval #
#########################################################################################################
sedibaStart <- sediba - runif(nrep) * temporalRange
sedibaEnd <- sedibaStart + temporalRange

#sample fictitious sediba fossil dates
sampleSediba <- runif(n = nrep, min = sedibaStart, max = sedibaEnd)

#compute 1-tailed p-value assuming the observed earliest Homo is truly the earliest Homo
sum((sampleSediba - earliestHomo > obsDiff) & (sedibaStart < earliestHomo)) / nrep

################################################################################################################################
# what if we also sample Homo? assuming the observed earliest Homo is just a single sample from the total temporal range of Homo
################################################################################################################################

#We can constrain our temporal range with the same sliding mechanism as we did sediba #
HomoStart <- earliestHomo - runif(nrep) * temporalRange
HomoEnd <- HomoStart + temporalRange

#sample fictitious Homo fossil dates
sampleHomo <- runif(n = nrep, min = HomoStart, max = HomoEnd)

#compute 1-tailed p-value assuming the observed earliest Homo is just a single sample from the total temporal range of Homo
sum((sampleSediba - sampleHomo > obsDiff) & (sedibaStart < HomoStart)) / nrep

###############################################################################################################################
# what if we also sample Homo? but sample multiple Homo, and then see if the oldest one is at least obsDiff younger than sediba
###############################################################################################################################

#sample fictitious Homo fossil dates
nHomo <- 10
sampleHomo <- sapply(1:nrep, function(x) min(runif(n = nHomo, min = HomoStart[x], max = HomoEnd[x])))

#compute 1-tailed p-value assuming the observed earliest Homo is just a single sample from the total temporal range of Homo
sum((sampleSediba - sampleHomo > obsDiff) & (sedibaStart < sampleHomo)) / nrep

################################################################################################################################
# what if we also sample Homo? but sample multiple Homo, and then set the currently observed Homo to be the youngest sample Homo
################################################################################################################################

#sample fictitious Homo fossil dates
nHomo <- 10
sampleHomo <- sapply(1:nrep, function(x) min(runif(n = nHomo, min = HomoStart[x], max = HomoEnd[x])))
HomoStart2 <- HomoStart + (earliestHomo - sampleHomo)
HomoEnd2 <- HomoStart + temporalRange
sampleHomo <- sapply(1:nrep, function(x) min(runif(n = nHomo, min = HomoStart2[x], max = HomoEnd2[x])))


#compute 1-tailed p-value assuming the observed earliest Homo is just a single sample from the total temporal range of Homo
sum((sampleSediba - sampleHomo > obsDiff) & (sedibaStart < sampleHomo)) / nrep
