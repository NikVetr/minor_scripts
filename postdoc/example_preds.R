set.seed(1)
teams <-  c("Animal", "GHD", "XST", "Survey", "WIT", "IAPS", "Ops", "SP")
cols <- c("#307496", "#963074", "#749630", "#965230", "#E4764B", "#6CE44B", "#303a40", "#bde6f8")
n <- length(teams)
mus <- rexp(n) * (rbinom(n, 1, 0.3) - 1/2) + 0.5
sds <- rlnorm(n, mean = -2, sd = 0.2)
obs <- lapply(1:n, function(i) rnorm(1E2, mus[i], sds[i]) * 100)
rx <- range(c(obs))
rx <- (rx - mean(rx)) * 1.1 + mean(rx)
kdes <- lapply(obs, density, from = rx[1], to = rx[2])
kdesx <- sapply(kdes, function(x) x$x)
kdesy <- sapply(kdes, function(y) y$y)
plot(NULL, NULL, xlim = range(kdesx), ylim = c(0, max(kdesy)), xlab = "Net Budget Surplus or Deficit (Thousands USD)", 
     ylab = "", yaxt = "n", xaxt = "n", frame = F)
axis(1, col = "#337799", lwd = 2)
abline(v = 0, lty = 2, lwd = 2, col = adjustcolor("#337799", 0.3))
for(i in 1:n){
  polygon(x = c(kdesx[,i], rev(kdesx[,i])),
          y = c(kdesy[,i], rep(0, length(kdesx[,i]))), 
                col = adjustcolor(cols[i], 0.4), border = cols[i])
  maxy <- max(kdesy[,i])
  maxx <- kdesx[which.max(kdesy[,i]),i]
  text(x = maxx, y = maxy, pos = 3, col = cols[i], teams[i], xpd = NA)
}

plot(NULL, NULL, xlim = range(kdesx), ylim = c(0, max(kdesy)), xlab = "Net Budget Surplus or Deficit (Thousands USD)", 
     ylab = "", yaxt = "n", xaxt = "n", frame = F)
axis(1, col = "#337799", lwd = 2)
prob_neg <- sapply(obs, function(x) mean(x < 0))
for(i in 1:n){
  negx <- kdesx[,i] < 0
  posx <- !negx
  polygon(x = c(kdesx[negx,i], rev(kdesx[negx,i])),
          y = c(kdesy[negx,i], rep(0, length(kdesx[negx,i]))), 
          col = adjustcolor("red", 0.4), border = "red")
  polygon(x = c(kdesx[!negx,i], rev(kdesx[!negx,i])),
          y = c(kdesy[!negx,i], rep(0, length(kdesx[!negx,i]))), 
          col = adjustcolor("green", 0.4), border = "green")
  maxy <- max(kdesy[,i])
  maxx <- kdesx[which.max(kdesy[,i]),i]
  text(x = ifelse(i == which.max(prob_neg), maxx - 18, maxx), 
       y = maxy, pos = 3, col = c("red", "green")[(prob_neg[i] < 0.5)+1], 
       labels = ifelse(i == which.max(prob_neg), 
                       paste0("Pr(deficit) = ", prob_neg[which.max(prob_neg)]), 
                       prob_neg[i]), 
       xpd = NA)
}
abline(v = 0, lty = 2, lwd = 2, col = adjustcolor("#337799", 1))

