logit <- function(p) log(p / (1-p))
invlogit <- function(x) exp(x) / (1+exp(x))

n <- 15^2
opts <- 0:10
conf_scores <- sample(opts, n, T)
lies <- rbinom(n, 1, conf_scores/10)
d <- data.frame(conf = conf_scores, lies = lies)
hist(d$conf[d$lies == 0])
hist(d$conf[d$lies == 1])

prop_lies <- sapply(opts, function(opt) mean(d$lies[d$conf == opt]))
par(mar = c(4,4,1,6), mfrow = c(1,2))
plot(NA, xlim = range(opts) + c(-1,1), ylim = c(0,1), xlab = "confidence in lie", ylab = "proportion")
for(opt in opts){
  rect(xleft = opt-1/2, xright = opt+1/2,
       ybottom = 0, ytop = prop_lies[opt+1], col = "#5d76cb")
  rect(xleft = opt-1/2, xright = opt+1/2,
       ybottom = prop_lies[opt+1], ytop = 1, col = "darkorange")
}
legend(x = 12, y = 1, legend = c("lies", "truths"), col = c("#5d76cb", "darkorange"), pch = 15, cex = 1, xpd = NA)

num_conf <- lapply(c(0,1), function(opt) table(d$conf[d$lies == opt]))
max_opt <- max(table(d$conf))+10
plot(NA, xlim = range(opts) + c(-1,1), ylim = c(0,max_opt), xlab = "confidence in lie", ylab = "count",
     yaxt = "n")
axis(2, at = seq(max_opt, max_opt - max(num_conf[[1]]), by = -2), 
     labels = seq(0, max(num_conf[[1]]), by = 2))
axis(4, at = seq(0, max(num_conf[[2]]), by = 2), 
     labels = seq(0, max(num_conf[[2]]), by = 2))
for(opt in opts){
  rect(xleft = opt-1/2, xright = opt+1/2,
       ybottom = 0, ytop = num_conf[[2]][as.character(opt)], col = "#5d76cb")
  rect(xleft = opt-1/2, xright = opt+1/2,
       ybottom = max_opt - num_conf[[1]][as.character(opt)], ytop = max_opt, col = "darkorange")
  
}
legend(x = 12, y = max_opt, legend = c("lies", "truths"), col = c("#5d76cb", "darkorange"), pch = 15, cex = 1, xpd = NA)
