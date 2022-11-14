library(RMTstat)

#top level parameters
incl_dir_samps <- F
n_obs <- 2E2
k <- 2E1
wvar <- 1

#test I remember stuff (later blocks of lines should be equal)
x <- matrix(rnorm(n_obs * k), n_obs, k)
x_m0 <- t(t(x) - apply(x, 2, mean))
x_m0_sd1 <- t(t(x_m0) / apply(x_m0, 2, sd) * sqrt(wvar))

sqrt(eigen(cov(x))$values)
sqrt(eigen(t(x_m0) %*% x_m0 / (n_obs-1))$values)
svd(x_m0 / sqrt(n_obs-1))$d

sqrt(eigen(cor(x))$values)
sqrt(eigen(t(x_m0_sd1) %*% x_m0_sd1 / (n_obs-1))$values)
svd(x_m0_sd1 / sqrt(n_obs-1))$d

#simulate SVs
n_rep <- 1E3
svals_samp <- t(replicate(n_rep, {
  x_rep <- matrix(rnorm(n_obs * k, sd = sqrt(wvar)), n_obs, k)
  sqrt(eigen(cov(x_rep))$values)
  # svd(t(t(x_rep) - apply(x_rep, 2, mean)) / sqrt(n_obs-1))$d #equivalently
}))
kdens <- density(svals_samp^2, kernel = "gaussian") #can use change of variables to respect bounds
kdens_by_sv <- lapply(1:k, function(i) density(svals_samp[,i]^2, kernel = "gaussian"))


#find distribution analytically
bounds <- qmp(c(0, 1), ndf = n_obs, pdim = k, var = wvar)
svals <- seq(bounds[1], bounds[2], length.out = 101)
mpd <- dmp(svals, ndf = n_obs, pdim = k, var = wvar)

#get specified decision threshold?
lambda <- k / n_obs 
thresh <- wvar * (1 + sqrt(lambda))^2

#plot to compare
par(mfrow = c(2,1))
plot(NA, NA, type = "l", xlim = range(c(kdens$x, svals)), ylim = range(c(kdens$y, mpd)), 
     ylab = "density", xlab = latex2exp::TeX("singular values ($^{2}$)"),
     main = latex2exp::TeX(paste0("n$_{df}$ = ", n_obs, ", p$_{dim}$ = ", k, ", $\\sigma^2$ = ", wvar)))
polygon(svals, mpd, col = adjustcolor(2, 0.5))
polygon(kdens$x, kdens$y, col = adjustcolor(4, 0.5))
text(range(svals_samp^2), rep(par("usr")[3], 2), labels = paste0("\u2191"), cex = 3, pos = 1, xpd = NA)
text(range(svals_samp^2), rep(par("usr")[3]-strheight("\u2191", cex = 3) * 1.4, 2), 
     labels = paste0("obs. ", c("min", "max"), "."), cex = 1.5, font = 2, pos = 1, xpd = NA)
text(thresh, kdens$y[which.min(abs(kdens$x - thresh))] + diff(par("usr")[3:4])/ 20, labels = paste0("\u2193"), cex = 3, pos = 3, xpd = NA)
text(thresh, kdens$y[which.min(abs(kdens$x - thresh))] + diff(par("usr")[3:4])/ 20 + strheight("\u2193", cex = 3) * 1.1, 
     labels = "decision\nthreshold", cex = 1, pos = 3, xpd = NA)



legend("topright", legend = latex2exp::TeX(c("Marchenko-Pastur Distribution", "Monte-Carlo Sampled SV$^2$s KDE")),
       pt.bg = adjustcolor(c(2, 4), 0.5), pch = 22, pt.cex = 2, bty="n")

#draw kdes for each individual SV^2
cols <- viridis::viridis(k)
plot(NA, NA, type = "l", 
     xlim = range(unlist(do.call(rbind, kdens_by_sv)[,"x"])),
     ylim = range(unlist(do.call(rbind, kdens_by_sv)[,"y"])) / k,
     ylab = "density", xlab = latex2exp::TeX("singular values ($^{2}$)"))
for(i in 1:k){
  polygon(kdens_by_sv[[i]]$x, kdens_by_sv[[i]]$y / k, col = adjustcolor(cols[i], 0.5))
}
gradloc <- list(xleft = par("usr")[2] - diff(par("usr")[1:2]) / 50, 
                ybottom = par("usr")[4] - diff(par("usr")[3:4]) / 1.5, 
                xright = par("usr")[2] - diff(par("usr")[1:2]) / 100, 
                ytop = par("usr")[4] - diff(par("usr")[3:4]) / 15)
plotrix::gradient.rect(gradloc$xleft, gradloc$ybottom, gradloc$xright, gradloc$ytop, 
                       col = viridis::viridis(100), gradient = "y")
text(x = rep(gradloc$xleft, 2), y = gradloc[c("ytop", "ybottom")], 
     labels = latex2exp::TeX(paste0(c("smallest", "largest"), " SV$^2$")), pos = 2)

#sample directly from distribution and add to plot
if(incl_dir_samps){
  mp_samps <- rmp(n_rep, ndf = n_obs, pdim = k, var = wvar)
  hist(mp_samps, probability = T, add = T, col = adjustcolor(7, 0.5)) 
}


#follow-up eigenvector question
random_orthonormal_matrix <- function(k){
  v <- matrix(0, k, k)
  v1 <- rnorm(k)
  v1 <- v1 / sqrt(sum(v1^2))
  v[,1] <- v1
  v_sub <- t(t(v[,1]))
  
  for(i in 1:(k-1)){
    v_temp <- rnorm(k-i)
    if(i != (k-1)){
      z <- t(t(v_temp) %*% t(t(v[1:(k-i),1:i])))
    } else {
      z <- t(t(v_temp) %*% t(v[1:(k-i),1:i]))
    }
    v_temp <- c(v_temp, solve(t(v_sub[(k-i+1):k,])) %*% -z)
    v_temp <- v_temp / sqrt(sum(v_temp^2))
    v[,i+1] <- v_temp
    v_sub <- v[,1:(i+1)]
  }
  v
}

v <- random_orthonormal_matrix(k)
v %*% t(v)
det(v)
sapply(1:k, function(i) det(v[,i] %*% t(v[,i])))
#or alternatively
# correlation matrix determinant shrinks when correlations get bigger bc the determinant is the product of the eigenvalues and the geometric mean is always smaller than the arithmetic mean https://en.wikipedia.org/wiki/Inequality_of_arithmetic_and_geometric_means (since their sum is constrained, greater inequality of eigenvalues results in a smaller product)


v %*% t(v)
l <- sort(gtools::rdirichlet(1, rep(1,k)), T) * k
nR <- (v) %*% diag(l) %*% t(v)

# find eigenvectors given this basis, if any exist
specific_corr_mat <- T
if(specific_corr_mat){
  rij <- 0.9
  R <- diag(k) + rij - diag(k) * 0.9
  v <- eigen(R)$vectors
}
solve(sapply(1:k, function(i) diag(v[,i] %*% t(v[,i])))) %*% t(t(rep(1, k)))
apply(solve(sapply(1:k, function(i) diag(v[,i] %*% t(v[,i])))), 1, sum)

t(v) %*% R %*% v

#is generating bases in this way related to an lkj(1) any? eg looking at diagonal elements


#rotation matrices in n-dimensions
d2rot <- function(theta) matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
len <- function(x) sqrt(sum(x^2))

pts <- t(t(c(0.5, 0.9)))
plot(t(pts), pch = 19, cex = 3, xlim = c(-len(pts), len(pts)), ylim = c(-len(pts), len(pts))); abline(h=0); abline(v=0)
points(t(d2rot(pi/4) %*% x), pch = 19, cex = 3)

# k <- 50
dkrot <- function(k, thetas) {
  angles <- data.frame(cbind(t(combn(1:k, 2)), theta = thetas))
  rotmats <- lapply(angles$theta, d2rot)
  big_rotmats <- lapply(1:choose(k, 2), function(i){
    bigrm <- diag(k)
    bigrm[angles$V1[i], angles$V1[i]] <- rotmats[[i]][1,1]
    bigrm[angles$V1[i], angles$V2[i]] <- rotmats[[i]][1,2]
    bigrm[angles$V2[i], angles$V1[i]] <- rotmats[[i]][2,1]
    bigrm[angles$V2[i], angles$V2[i]] <- rotmats[[i]][2,2]
    bigrm
  })
  rotmat <- Reduce("%*%", big_rotmats)
  rotmat
}

k <- 20
hist(replicate(1000, dkrot(k, runif(choose(k, 2), 0, 2*pi))), breaks = 100)

k <- 20
rotmat <- dkrot(k, runif(choose(k, 2), 0, 2*pi))
# rotmat <- random_orthonormal_matrix(k)
rotv <- eigen(rotmat)$vectors
rotl <- eigen(rotmat)$values
sum(rotv[,1]^2)
Re((rotv) %*% diag(rotl) %*% t(rotv)) - rotmat
round(Re(t(rotv) %*% rotv))
round(Re((rotv) %*% t(rotv)))
round(Im(t(rotv) %*% rotv), 3)

eigen(cor(x))$vectors %*% diag(eigen(cor(x))$values) %*% t(eigen(cor(x))$vectors) - cor(x)
t(eigen(cor(x))$vectors) %*% diag(eigen(cor(x))$values) %*% (eigen(cor(x))$vectors) - cor(x)
eigen(cor(x))$vectors %*% t(eigen(cor(x))$vectors) 
t(eigen(cor(x))$vectors) %*% eigen(cor(x))$vectors

rv <- rotmat %*% v
rv %*% t(rv)

a <- (eigen(eigen(R)$vectors)$vectors)
round(Re(t(a) %*% a))
round(Im(t(a) %*% a))

#find eigenvectors given some eigenvalues by optimizing elements of the rotation matrix
evectors_from_evalues <- function(data, par) {
  
  #undo params
  thetas <- par

  #undo data
  k <- data$k
  l <- data$l
  
  #find rotmat
  mat <- dkrot(k, thetas)
  
  #compose corrmat
  cmat <- mat %*% l %*% t(mat)
  
  #return log diag prod (optim minimizes by default, so take negative)
  return(-sum(log(diag(cmat))))
  
}

k <- 5
l <- diag(sort(gtools::rdirichlet(1, rep(1,k)), T) * k)
l <- diag(c(k-1, sort(gtools::rdirichlet(1, rep(1,k-1)), T)))
nrep <- 10
optim_outs <- t(replicate(nrep, {
  optim_out <- optimx::optimx(par = runif(choose(k, 2), 0, pi), 
                 fn = evectors_from_evalues, 
                 data = list(k = k, l = l), method = "nlminb", 
                 lower = rep(0, choose(k,2)), 
                 upper = rep(pi, choose(k,2)), #not 2*pi cos of identifiability
                 control = list(maxit = 1E3, trace = 0, dowarn = F))
  as.numeric(optim_out[1,paste0("p", 1:choose(k,2))])
}))

est_vs <- lapply(1:nrep, function(i) dkrot(k, optim_outs[i,]))
est_Rs <- lapply(1:nrep, function(i) est_vs[[i]] %*% l %*% t(est_vs[[i]]))
est_Rcs <- lapply(1:nrep, function(i) cov2cor(est_Rs[[i]]))
est_Rcs_ij <- sapply(1:nrep, function(i) est_Rcs[[i]][upper.tri(est_Rcs[[i]])])


sapply(est_Rcs, function(x) det(x)) #as expected
apply(est_Rcs_ij, 2, function(x) mean(abs(x)))

#confirm evs are correct
diag(l)
t(sapply(1:nrep, function(i) eigen(est_Rcs[[i]])$values))

#now do the inverse question -- can we solve for eigenvalues of an arbitrary orthonormal matrix? 

#first check rotation matrices and random orthonormal matrix functions are the same
k <- 20
nrep <- 1E2
hist(c(replicate(nrep, random_orthonormal_matrix(k))), breaks = 100, freq = F)
# hist(c(replicate(nrep, random_orthonormal_matrix(k)[lower.tri(diag(k))])), breaks = 100, freq = F)
# hist(c(replicate(nrep, random_orthonormal_matrix(k)[upper.tri(diag(k))])), breaks = 100, freq = F)
# hist(c(replicate(nrep, diag(random_orthonormal_matrix(k)))), breaks = 100, freq = F)
hist(c(replicate(nrep, dkrot(k, runif(choose(k, 2), 0, 2 * pi)))), breaks = 100, freq = F)
#interestingly, no!

evalues_from_evectors <- function(data, par) {
  
  #undo data
  k <- data$k
  v <- data$v
  
  #undo params
  l <- rev(cumsum(c(1, par)))
  l <- diag(l / sum(l) * k)
  
  #compose corrmat
  cmat <- v %*% l %*% t(v)
  
  #return log diag prod (optim minimizes by default, so take negative)
  return(-sum(log(diag(cmat))))
  
}

k <- 5
v <- random_orthonormal_matrix(k)
v <- dkrot(k, runif(choose(k, 2), 0, 2 * pi))
v <- eigen(cor(x)[1:k, 1:k])$vectors
nrep <- 10
optim_outs <- t(replicate(nrep, {
  optim_out <- optimx::optimx(par = runif(k-1, 0, 1), 
                              fn = evalues_from_evectors, 
                              data = list(k = k, v = v), method = "nlminb", 
                              lower = rep(0, k-1), 
                              upper = rep(1, k-1), 
                              control = list(maxit = 1E3, trace = 0, dowarn = F))
  as.numeric(optim_out[1,paste0("p", 1:(k-1))])
}))

est_ls <- lapply(1:nrep, function(i) {
  l <- rev(cumsum(c(1, optim_outs[i,])))
  l <- diag(l / sum(l) * k)
})
est_Rs <- lapply(1:nrep, function(i) round((v %*% est_ls[[i]] %*% t(v)) * 1000) / 1000)
est_Rs

#solve for a below
sapply(2:nrep, function(i) (diag(est_ls[[1]]) - 1) / (diag(est_ls[[i]]) - 1))


est_Rcs <- lapply(1:nrep, function(i) cov2cor(est_Rs[[i]]))
est_Rcs_ij <- sapply(1:nrep, function(i) est_Rcs[[i]][upper.tri(est_Rcs[[i]])])
apply(est_Rcs_ij, 2, function(x) mean(abs(x)))




#check row reshuffling question
R <- cor(x)
# rij <- 0.9
# R <- diag(k) + rij - diag(k) * 0.9
L <- eigen(R)$values
V <- eigen(R)$vectors

V %*% diag(L) %*% t(V) - R
V %*% diag(k) %*% t(V)

a <- 1 / (1-min(L))
a <- 1 / (1-max(L)) #weight for original corr
b <- 1 - a #weight for identity
nR <- R * a + diag(k) * b
nR

plot(eigen(nR)$vectors, V)
abline(0,1)

plot(eigen(nR)$vectors, V[,k:1])
abline(0,1)


V %*% ((diag(L) * a + diag(k) * b)) %*% t(V) 


par(mfrow = c(3,1))
hist(nR[upper.tri(nR)])
hist(R[upper.tri(R)])
eigen(nR)$VaLues


V %*% ((diag(L) * a + diag(k) * b)) %*% t(V) 
cor(x) * a + diag(k) * b


diag(cor(x) * a + diag(k) * b)
eigen(cor(x) * a + diag(k) * b)

Reduce("+", lapply(1:k, function(i) v[,i] %*% t(v[,i])))
max(abs(Reduce("+", lapply(1:k, function(i) v[,i] %*% t(v[,i]) * l[i])) - cor(x)))

reord <- sample(1:k)
v1 <- v[reord,]
((v1) %*% diag(l) %*% t(v1))[order(reord), order(reord)] - cor(x)

P <- model.matrix(~0+as.factor(sample(1:k)))
v2 <- P %*% v 

(v2) %*% diag(l) %*% t(v2) - P %*% cor(x) %*% t(P)


(v2) %*% diag(l) %*% t(v2)

cor(x) %*% t(t(v[,1])) - t(t(l[1] * v[,1]))


plot((cor(x) %*% v[,1]), v[,1] * l[1]
)
