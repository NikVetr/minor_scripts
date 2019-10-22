g <- c(13.96, 14.04, 13.88, 11.51, 10.37, 12.49, 12.57, 9.18, 13.39, 12.33, 14.24, 11.59, 12.9, 12.73, 13.88, 12, 12.41, 11.51, 11.84, 13.06, 13.76, 13.22, 11.59, 10.9, 12.65, 11.76, 7.59, 9.96, 10.12, 11.51, 12.08, 9.88, 8.9, 13.71, 13.47, 11.67, 12.65, 13.55, 10.94, 7.43)
g <- g / 16 * 100 + 10
hist(g, freq = F, xlim = c(40,100), ylim = c(0,0.075), at = 4:20*5, density = 0, ylab = "",
     xlab = "Score (Raw %)", main = "", breaks = 20)
d <- density(g)
lines(d$x, d$y, lwd = 2, col = 2)

abline(v = mean(g), lwd = 3, lty = 2)
text(x = mean(g)-0.75, y = 0.0735, "mean", srt = 90, cex = 1.3)

for(i in 1:10){
  abline(v = quantile(g, 1:10/10)[i], lwd = 1.5, lty = 3)
  text(x = quantile(g, 1:10/10)[i]-0.75, y = 0.075, paste0(i*10, "%"), srt = 90, cex = 1.1)
  
}

title("Curved Scores (+ 10%)", cex.main = 3)