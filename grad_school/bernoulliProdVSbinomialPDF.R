nrep <- 1e3
log <- T
l <- 1e3
p <- 0.3
p1 <- p2 <- rep(0, nrep)
ps <- 1:nrep/nrep
for(i in 1:nrep){
  d <- sample(c(1,0), size = l, replace = T, prob = c(p, 1-p))
  p1[i] <- dbinom(prob = ps[i], sum(d), l, log = log)
  p2[i] <- -sum(sapply(1:l, function(x) dbinom(prob = ps[i], d[x], 1, log = T)))
  if(!log){
    p2[i] <- exp(p2[i])
  }
}
dev.off()
plot(p1, p2, xlab = paste0(ifelse(log, "log ", ""), "binomial probability mass"), 
             ylab = paste0(ifelse(log, "log ", ""), "product of bernoulli probability masses"),
             main = paste0("number of draws = ", l, ", probability = ", p))
