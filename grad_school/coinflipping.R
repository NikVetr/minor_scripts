#You have 9 fair and 1 double headed coin
# You pick one randomly. You flip it 5 times. You observe 5 heads. What is probability to get head in a sixth flip?
nsim <- 1000000
coins <- rbinom(n = nsim, size = 1, prob = .9)
flips <- t(sapply(1:nsim, function (x) if(coins[x] == 1){rbinom(n = 6, size = 1, prob = .5)} else {rbinom(n = 6, size = 1, prob = 1)}))
sixth <- flips[which(sapply(1:nsim, function(x) all(flips[x,1:5] == c(1,1,1,1,1)))),][,6]
sum(sixth)/length(sixth)
1 - (1 - (.1 / ((.9 * .5^5) + .1)) ) / 2
