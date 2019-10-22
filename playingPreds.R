#do some eda with the graph
preds <- matrix(nrow = 7, ncol = 3, data = c(.5,.6,.7,.8,.9,.95,.99,8,22,13,18,25,9,4,5,14,8,16,24,8,4))
par(mfrow = c(1,1))
plot(preds[,1], preds[,3]/preds[,2], xlim = c(.5, 1), ylim = c(0,1), type = "l", col = 2, lwd = 4, 
     xlab = "true probability", ylab = "realized frequency", main = "A Plot of Possible Worlds")
for(i in 1:2000){
  lines(preds[,1],sapply(1:7, function(x) rbinom(n = 1, size = preds[x,2], prob = preds[x,1]))/preds[,2], type = "l", col = rgb(0,0,0,.01))
}       
lines(preds[,1], preds[,3]/preds[,2], xlim = c(.5, 1), ylim = c(.5,1), type = "l", col = 2, lwd = 4)

#monte carlo flavored NHST, weighted
test <- sapply(1:1e5, function(c) sum((sapply(1:7, function(x) rbinom(n = 1, size = preds[x,2], prob = preds[x,1]))/preds[,2] - preds[,1])^2)*(preds[,2]/sum(preds[,2])))
obs <- sum(((preds[,1] - preds[,3]/preds[,2])^2)*(preds[,2]/sum(preds[,2])))
p_val <- sum(obs < test)/length(test)
hist(test, breaks = 50, main = "Average Weighted Mean Squared Error Sampling Distribution", xlab = "weighted test statistic")
abline(v = obs, col= 2, lwd = 2)
text(x = 0.13, y = 40000, labels = paste0("1-tailed p-value = ", round(p_val, 3)))

#monte carlo flavored NHST, unweighted
test <- sapply(1:1e5, function(c) mean((sapply(1:7, function(x) rbinom(n = 1, size = preds[x,2], prob = preds[x,1]))/preds[,2] - preds[,1])^2))
obs <- mean(((preds[,1] - preds[,3]/preds[,2])^2))
p_val <- sum(obs < test)/length(test)
hist(test, breaks = 20, main = "Average Unweighted Mean Squared Error Sampling Distribution", xlab = "unweighted test statistic")
abline(v = obs, col= 2, lwd = 4)
text(x = 0.08, y = 7800, labels = paste0("1-tailed p-value = ", round(p_val, 3)))

#testing deviation
#monte carlo flavored NHST, weighted
test <- sapply(1:1e5, function(c) sum((sapply(1:7, function(x) rbinom(n = 1, size = preds[x,2], prob = preds[x,1]))/preds[,2] - preds[,1])^2)*(preds[,2]/sum(preds[,2])))
obs <- sum(((preds[,1] - preds[,2]/preds[,2])^2)*(preds[,2]/sum(preds[,2])))
p_val <- sum(obs < test)/length(test)
p_val
hist(test, breaks = 50, main = "Average Weighted Mean Squared Error Sampling Distribution", xlab = "weighted test statistic")
abline(v = obs, col= 2, lwd = 2)
text(x = 0.13, y = 40000, labels = paste0("1-tailed p-value = ", round(p_val, 3)))

#monte carlo flavored NHST, unweighted
test <- sapply(1:1e5, function(c) mean((sapply(1:7, function(x) rbinom(n = 1, size = preds[x,2], prob = preds[x,1]))/preds[,2] - preds[,1])^2))
obs <- mean(((preds[,1] - preds[,2]/preds[,2])^2))
p_val <- sum(obs < test)/length(test)
p_val
hist(test, breaks = 20, main = "Average Unweighted Mean Squared Error Sampling Distribution", xlab = "unweighted test statistic")
abline(v = obs, col= 2, lwd = 4)
text(x = 0.08, y = 7800, labels = paste0("1-tailed p-value = ", round(p_val, 3)))

sum(((preds[,1] - preds[,2]/preds[,2])^2)*(preds[,2]/sum(preds[,2])))

p <- c(rep(c(.5,.7), 5), .5)
d <- rbinom(n = 11, size = 1, prob = p)
mean((d-p)^2)

(.5 - sum(d[p == .5])/length(d[p == .5]))^2*6/11 + (.7 - sum(d[p == .7])/length(d[p == .7]))^2*5/11

sum()