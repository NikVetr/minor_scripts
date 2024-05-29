mu <- 3
sig <- 5
lambda <- 0.7
n <- 1E6
x <- rnorm(n, mu, sig)
sds <- rexp(n, rate = lambda)
y <- rnorm(n, x, sds)

varexp <- 1/lambda^2
expexp <- 1/lambda

var(y)
# hmm var of diff betw two normals is 2*sig^2
#if drawing a single constant sds, var of mixture would be sig^2 + sds^2
#but if we have a scale mixture, we need to add 2*var(sds), so 2 * 1/lambda^2
#but we also need a term for the mean in there somewhere? though as n -> inf this converges
var(y)
sig^2 + 2*varexp
sig^2 + varexp + expexp^2 #or equivalently (see derivation)

#let's try it out for other scale distributions
siglog <- 0.7
mulog <- 2
sds <- (rlnorm(n, mulog, siglog))
y <- rnorm(n, x, sds)


varlogn <- (exp(siglog^2)-1)*exp(2*mulog+siglog^2)
explogn <- exp(mulog + siglog^2/2)

mean(sds)
explogn

var(sds)
varlogn

#wait, so if the lognormal W is the distribution of sds
#then its contribution to the mixture variance is E[W^2]?
#bc for independent random variable X & Y, E[XY] = E[X]E[Y]
#the scale mixture is the same as taking our initial normal for the means
#and adding (a standard normal X) * (a scale RV sqrt(W)), where W is the dist
#of the variance.

#if W is the dist of the sds, we instead need to add Var(W) + E[W]^2
#see written note for full derivation, 
#but relies on identity E[XY] = Cov(X,Y) + E[X]E[Y]

var(y)
sig^2 + varlogn + explogn^2

#let's confirm with one more distribution, a gamma
gshape = 1
gscale = 1/0.7
sds <- rgamma(n, shape = gshape, scale = gscale)
y <- rnorm(n, x, sds)
vargam <- gshape*gscale^2
expgam <- gshape*gscale

var(y)
sig^2 + vargam + expgam^2

#woot yes it is correct

#repeating above but sampling variances and not SDs

#gamma mixture
sds <- sqrt(rgamma(n, shape = gshape, scale = gscale))
y <- rnorm(n, x, sds)
vargam <- gshape*gscale^2
expgam <- gshape*gscale
var(y)
sig^2 + expgam

#lnormal mixture
sds <- sqrt(rlnorm(n, mulog, siglog))
y <- rnorm(n, x, sds)
varlogn <- (exp(siglog^2)-1)*exp(2*mulog+siglog^2)
explogn <- exp(mulog + siglog^2/2)
var(y)
sig^2 + explogn

#exp mixture
sds <- sqrt(rexp(n, rate = lambda))
y <- rnorm(n, x, sds)
varexp <- 1/lambda^2
expexp <- 1/lambda
var(y)
sig^2 + expexp

#so for empirical Bayes, is posteriors are roughly all normal w/ difft means and variances
#we don't care about the means bc they are in the total variance already
#but we need to compute the distribution of their sds and then 
#subtract out the (sample variance of that dist) and the (sample expectation of that dist)^2
#or compute their variances and subtract out the sample expectation of those variances