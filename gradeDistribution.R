g <- c("47	63.75	64	55.75	67	38	58.75	38.5	56.25	48	66.85	47.05	40	64.7	67.5	51.25	61.15	56.5	49.75	50	55.75	67.75	60.25	53.5	53	66.25	43.75	52	50.75	63.3	52	58.75	52.5	57.5	59.75	63.75	66.3	36.75	58	44.25	54.75	52.75	49.5	59.75	59	68.25	61.5	59.25	47	57.5	44.5	55.25	63.75	31	58.8	53.5	63.5	59.25	38.5	56.15	63.75	54.5	58	41	65	52	49.5	56.75	57.75	49	63	61.25	52.5	56.25	51	60	57.75	60.5	60.25	55.3	62.75	52.5	53.5	54.4	67.35	55	53.5	60.75	40	64	61.75	48.75	57")
g <- as.numeric(strsplit(g, "\t")[[1]])
g <- g / 68.5 * 100 
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

title("Midterm I Grade Distribution", cex.main = 3)