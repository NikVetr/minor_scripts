#generate fictional "empirical" dataset to test

#sample 'A's from uniform in interval 0 to 1 (can rescale to any interval by multiplying by total length)
nA <- 25
A <- rbeta(n = nA, shape1 = 1, shape2 = 1)

#sample 'B's from beta with expectations drawn from As
abratio <- function(expectation){return((1-expectation)/expectation)}
a <- 50
b <- a * abratio(A)
nB <- 37
B <- rbeta(n = nB, shape1 = a, shape2 = sample(b, nB, T))

#calculate average minimum distance of each B from A
avgMinDist <- mean(sapply(1:length(B), function(x) min(abs(A - B[x]))))

#get sampling distribution from null model where A, B are iid ~unif(0,1)
nSamp <- 1E4
sampDistr <- sapply(1:nSamp, function(x) {A <- runif(n = nA, 0, 1); 
                                          B <- runif(n = nB, 0, 1); 
                                          mean(sapply(1:length(B), function(x) min(abs(A - B[x]))))})
p.value <- sum(avgMinDist > sampDistr)/nSamp
p.value
plot(density(sampDistr)); abline(v = avgMinDist, col = 2, lwd = 2)


