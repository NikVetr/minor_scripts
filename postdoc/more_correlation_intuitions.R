
#sample two random vectors in n-d
n <- 3
u <- rnorm(n)
v <- rnorm(n)

#center and scale
x <- u - mean(u) 
x <- x / norm(x, "2")
y <- v - mean(v)
y <- y / norm(y, "2")

#check correlations
cor(x, y)
sum(x * y)
# sum(x * y) / norm(x, "2") / norm(y, "2")
norm(x, "2") 
norm(y, "2")
cor_xy <- c(cor(x, y))

#introduce a new vector, z, that has some correlation theta with y
r_yz <- 0.75
theta <- acos(r_yz)
vecs2theta <- function(a,b) acos( sum(a*b) / ( sqrt(sum(a^2)) * sqrt(sum(b^2)) ) )
rotate_vec <- function(y, z, theta) {
  dot_product <- sum(y * z)
  cross_product <- c(y[2] * z[3] - y[3] * z[2],
                     y[3] * z[1] - y[1] * z[3],
                     y[1] * z[2] - y[2] * z[1])
  z * cos(theta) +
    cross_product * sin(theta) +
    y * dot_product * (1 - cos(theta))
}

# Generate a random vector orthogonal to y
o <- rnorm(n-1)
o <- c(o, -sum(y[1:(n-1)] * o) / y[n])
o <- o / sqrt(sum(o^2))
o <- o / sd(o) / sqrt(2)
# cov(cbind(y, o))
#o <- c(-y[2], y[1], 0) / sqrt(sum(y[1:2]^2)) #alternatively

# compute cholesky factor and transform to get z0
R_yz <- diag(2) + r_yz - diag(2) * r_yz
L_yz <- t(chol(R_yz))
yz0 <- t(L_yz %*% t(cbind(y, o)))
z0 <- yz0[,2]
z0 <- z0 - mean(z0) 
z0 <- z0 / norm(z0, "2")

#confirm properties
vecs2theta(z0,y)
theta
cor(z0,y)
r_yz

#try rotating vector around y?
z1 <- rotate_vec(y, z0, 0.5)
sum(z1^2)
vecs2theta(z1,y)
theta
#still on the unit hypersphere
cor(z1,y)
r_yz
#but no longer has the same correlation
mean(z1)
#also no longer has mean 0

#ok, so we need to get off the unit hypersphere
#other points on the hypersphere
ts <- seq(0,2*pi,length.out=100)
zrs <- sapply(ts, function(t) cor(y, rotate_vec(y, z0, t)))
plot(ts, zrs, type = "l")
zts <- sapply(ts, function(t) vecs2theta(y, rotate_vec(y, z0, t)))
plot(ts, zts, type = "l")

#so only two points (antipodes) in the rotation that satisfy this property
#also means that all of the other points in the cone do not have the desired correlation

#what does the *correct* shape look like on the unit hypersphere? let's sample and find out
R_yz <- diag(2) + r_yz - diag(2) * r_yz
L_yz <- t(chol(R_yz))
zs <- lapply(1:1E3, function(i){
  
  #sample random orthogonal vector
  index <- sample(1:3, 1)
  o <- rep(0,3)
  o[-index] <- rnorm(2)
  o[index] <- -sum(y[-index] * o[-index]) / y[index]
  o <- o / sqrt(sum(o^2))
  o <- o / sd(o) / sqrt(2)
  o
  
  # compute cholesky factor and transform to get z0
  yz0 <- t(L_yz %*% t(cbind(y, o)))
  z0 <- yz0[,2]
  z0 <- z0 - mean(z0) 
  z0 <- z0 / norm(z0, "2")
  z0
  
})
zs <- do.call(rbind, zs)
zs
#hm wait *are* there only two vectors that satisfy this property? no there can't be
#there need to be lots, pointing in different directions too, to allow for different
#correlations with x

#let's try again
zs <- lapply(1:1E3, function(i){
  
  #sample random orthogonal vector
  rz <- matrix(rnorm(6), 3, 2)
  rz <- t(solve(t(chol(cov(rz)))) %*% t(rz))
  
  # compute cholesky factor and transform to get z0
  # rz <- rz - (rz[,1] - y)
  # z0 <- rz[,2]
  # z0 <- z0 - mean(z0) 
  # z0 <- z0 / norm(z0, "2")
  # cor(z0, y)
  
})
zs <- do.call(rbind, zs)
zs

# test correlation matrix / cholesky factor relationship?
p <- 5
r <- 0.5
R <- diag(p) + r - diag(p) * r
L <- t(chol(R))
L
vecs2theta(L[,1], L[,2])
vecs2theta(L[,3], L[,2])

#construct set of mean = 0 unit vectors with desired properties
p <- 5 #number of distinct features (aka number of vectors)
n <- 6 #number of dimensions to the space (aka obs per vector)
vecs <- matrix(0, nrow = n, ncol = p)
vecs[1:2,1] <- c(-sqrt(1/2),sqrt(1/2)) #ensures mean 0 and unit length constraint

#next dimension
i<-2
l2 <- R[1,2]
# w2 <- sin(acos(R[1,2])) * c(rep(0,i),1,rep(0,n-i-1))
w2 <- sqrt(1-R[1,2]^2)
w2 <- c(0,0,sqrt(w2^2/2) * c(-1,1),0,0)
vecs[,i] <- l2 * vecs[,i-1] + w2
# norm(vecs[,i], "2")
# mean(vecs[,i])
# sum(vecs[,i] * vecs[,i-1])

#ok, the trick is to represent each new unit vector as a combination of the previous
#ones plus an orthogonal component, w
#so v3 = l1 * v1 + l2 * v2 + w3
#we can substitute back in to the angle constraint
#v3 dot v1 = R1,3
#v3 dot v2 = R2,3
#to get a linear system R[1:(i-1),1:(i-1)] * l_d = R[1:(1-i),i]
#which we can solve in the usual way
i <- 3
l <- solve(R[1:(i-1),1:(i-1)]) %*% R[i,1:(i-1)]
#then we solve for the orthogonal component, w
#so if v is a summation of the previous vectors weighted by l
lvs <- vecs[,1:(i-1)] %*% l
#w needs to be added to this to satisfy v1 * w = 0, v2 * w = 0, ... v_{i-1} * w = 0
#and also |lvs + w| = 1, and since all the earlier vs have norm 1, we get
w_ss <- 1 - sum(l^2)
#... we seem to run out of dimensions, so maybe the initial approach was off

#let's try again
p <- 5 #number of distinct features (aka number of vectors)
n <- 8 #number of dimensions to the space (aka obs per vector)
r <- 0.5
R <- diag(p) + r - diag(p) * r
# R <- cor(matrix(rnorm(n*p), n, p))
vecs <- matrix(0, nrow = n, ncol = p)
vecs[1:2,1] <- c(-sqrt(1/2),sqrt(1/2))
vecs[,1] <- rnorm(n) #or can initialize randomly
vecs[,1] <- (vecs[,1] - mean(vecs[,1])) / norm(vecs[,1] - mean(vecs[,1]), "2")

null_space <- function(mat) {
  qr_decomp <- qr(mat)
  ns <- qr.Q(qr_decomp, complete = TRUE)[, (ncol(mat) + 1):nrow(mat), drop = F]
  return(ns)
}

generate_null_vector <- function(prev_vectors) {
  n <- nrow(prev_vectors)
  
  # Start with a random vector
  w <- rnorm(n)
  
  # Project out components in the span of previous vectors
  if (ncol(prev_vectors) > 0) {
    proj <- prev_vectors %*% solve(t(prev_vectors) %*% prev_vectors) %*% t(prev_vectors) %*% w
    w <- w - proj
  }
  
  # Normalize to unit norm
  w <- w / sqrt(sum(w^2))
  
  return(w)
}


for(i in 2:p){
  l <- solve(R[1:(i-1),1:(i-1)]) %*% R[i,1:(i-1)]
  lv <- vecs[,1:(i-1)] %*% l
  w_ss <- 1 - sum(lv^2)
  # w_ns <- null_space(vecs[,1:(i-1), drop = F])[,1] #use first null space vector
  w_ns <- generate_null_vector(vecs[,1:(i-1), drop = F])
  w <- w_ns - mean(w_ns) 
  w <- w / norm(w, "2") * sqrt(w_ss)
  vecs[,i] <- lv + w  
}

apply(vecs, 2, mean) #all have mean 0
apply(vecs, 2, norm, "2") #all have unit norm
t(vecs) %*% vecs - R #appropriately yields correlation matrix (all vectors have right angles)

#also, the volume enclosed by the transpose of this vectors, relative to that
#enclosed by orthogonal vectors, is the determinant of the correlation matrix

#fals when p > n because we run out of degrees of freedom! 
#the last vector has n degrees of freedom initially
#but it loses one to the unit length constraint,
#one to the mean 0 constraint
#and n-1 to the angle constraints
#but that requires n+1 degrees of freedom


#the relationship to determinants
#the volume captured by these vectors is given by det(G) / nrow(G)!, where G=t(vecs) %*% vecs
#eg in 2D, if we have 2 vectors x1 = 0.5, sqrt(0.75), x2 = sqrt(0.5), sqrt(0.5) + the origin
x3d <- cbind(c(0.5, sqrt(0.75)),
             c(sqrt(0.5), sqrt(0.5)),
             c(0,0))
G <- x3d %*% t(x3d)
sqrt(det(G)) / factorial(nrow(G))
#herons formula for area is sqrt(s(s-a)(s-b)(s-c))), where s=(a+b+c)/2
sides <- c(1, 1, sqrt(sum((x3d[,1] - x3d[,2])^2)))
s <- sum(sides)/2
sqrt(s*prod(s-sides))

x3d_i <- cbind(c(0,1),
             c(1,0),
             c(0,0))
G_i <- x3d_i %*% t(x3d_i)
sqrt(det(G_i)) / factorial(nrow(G_i))
sides <- c(1, 1, sqrt(sum((x3d_i[,1] - x3d_i[,2])^2)))
s <- sum(sides)/2
sqrt(s*prod(s-sides))

#ok, so that demonstrates the hypervolume idea

#project to the lower dimensional subspace
qr_decomp <- qr(vecs)
Q <- qr.Q(qr_decomp)[, 1:(nrow(vecs)-2)]  # Only first (p-2) columns
vecs_projected <- t(Q) %*% vecs
as.matrix(dist(t(m_vecs))) - as.matrix(dist(t(vecs))) #distances are preserved
orig_angles <- t(vecs) %*% vecs
new_angles <- (t(m_vecs) %*% m_vecs) / outer(sqrt(colSums(m_vecs^2)), sqrt(colSums(m_vecs^2)), "*") #as are angular relationships
round(orig_angles - new_angles, 3)


G <- m_vecs %*% t(m_vecs)
sqrt(det(G)) / factorial(nrow(G))
det(R)

t(vecs) #so this is the transformation matrix from the original space to the lower-d subspace
det(rbind(t(vecs), 0)) #because it projects to a subspace, it squishes things to hypervolume 0
#The determinant of the Gram matrix of the projected basis vectors tells us 
#the fraction of the maximum area occupied relative to the full basis
# Mean 0 unit vectors represent basis in lower dimensional subspace, and the relative volume (to an orthogonal basis) of the paralleliped enclosed in that subspace can be obtained from the determinant of the gram matrix.
# This is exactly the determinant of the original correlation matrix
# The angles between our new basis vectors and the cholesky factor vectors are the same too

#can use these to transform uncorrelated data
nsamp <- 1E3
y <- matrix(rnorm(nsamp*(n)), ncol=n, nrow=nsamp)
yt <- t(t(vecs) %*% solve(chol(cov(y))) %*% t(y))
cor(yt)
R
y <- matrix(rnorm(nsamp*(p)), ncol=p, nrow=nsamp)
yt <- t(t(chol(R)) %*% solve(chol(cov(y))) %*% t(y))
cor(yt)
R


cov(t(solve(chol(cov(y))) %*% t(y))) #gotta be applying a bessel-y correction here
