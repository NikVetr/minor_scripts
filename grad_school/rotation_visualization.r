start <- c(1,1)
rot <- matrix(data = c(0.078,1.033,-1.01,-0.11), nrow = 2, byrow = F); det(rot)
rot <- matrix(data = c(0,1,-.93,-0.08)+matrix(rnorm(4, 0, 0.1), ncol = 2), nrow = 2, byrow = F); det(rot)
incr_size <- 0.01
# rot <- matrix(data = c((1-incr_size^2)^0.5,incr_size,-incr_size,(1-incr_size^2)^0.5), nrow = 2, byrow = F); det(rot)
post <- rot %*% start 
nrep <- 600
posts <- matrix(rep(post, nrep), ncol = 2, byrow = T)
for(i in 2:nrep){
posts[i,] <- rot %*% posts[i-1,]
}
plot(posts, type = "l")
plot(posts, ty = "s")
lines(posts)
