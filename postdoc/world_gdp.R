d <- read.csv("world-gdp-over-the-last-two-millennia.csv")
colnames(d) <- c("entity", "code", "year", "gdp")
d$log_gdp <- log10(d$gdp)

plot(d$year, d$gdp)
lines(d$year, d$gdp)

n <- 4
par(mfrow = c(n,1))

for(i in 1:n){
  nlogged_vals <- eval(parse(text = paste0(paste0(rep("log10(", i), collapse = ""), "d$gdp", paste0(rep(")", i), collapse = ""))))
  ylab <- latex2exp::TeX(paste0("$", paste0(rep("log_{10}(", i), collapse = ""), "gdp", paste0(rep(")", i), collapse = ""), ""))
  plot(d$year, nlogged_vals, xlab = "year", ylab = ylab, type = "l")
  points(d$year, nlogged_vals,  pch = 19, col = "white")
  points(d$year, nlogged_vals,  pch = 19, col = adjustcolor(1,0.5))
}

plot(d$year^20, d$gdp^(1/20))

par(mar = c(5,7,2,2), mfrow = c(1,1))
nvals_to_plot <- 9
plot(d$year[1:nvals_to_plot], d$gdp[1:nvals_to_plot], xlab = "year", ylab = "", pch = 19, col = adjustcolor(1,0.5), 
     yaxt = "n")
lines(d$year[1:nvals_to_plot], d$gdp[1:nvals_to_plot], xlab = "year", ylab = "world gdp", pch = 19, col = adjustcolor(1,0.5))
points(d$year[1:nvals_to_plot], d$gdp[1:nvals_to_plot],  pch = 19, col = "white")
points(d$year[1:nvals_to_plot], d$gdp[1:nvals_to_plot],  pch = 19, col = adjustcolor(1,0.5))


yrange <- range(d$gdp[1:nvals_to_plot])
yvals <- seq(yrange[1], yrange[2], length.out = 5)
yincr <- floor(diff(yvals)[1] / 10^floor(log10(diff(yvals)[1]))) *10^floor(log10(diff(yvals)[1]))
yvals <- seq(yrange[1], yrange[2], by = yincr)

newyr <- yrange / yincr
yvals <- seq(ceiling(newyr[1]), floor(newyr[2]), by = 1) * yincr
ymag <- floor(log10(yincr))
text(x = par("usr")[1], y = yvals + diff(par("usr")[3:4])/120, labels = latex2exp::TeX(paste0(yvals / 10^ymag, " x ", "$10^{", ymag, "}$")), pos = 2, xpd = NA)
segments(x0 = par("usr")[1], y0 = yvals, y1 = yvals, x1 = par("usr")[1] - diff(par("usr")[1:2])/80, xpd = NA)

text(x = par("usr")[1] - diff(par("usr")[1:2])/6, y = mean(yrange), labels = "world gdp (USD, 2011 equiv.)", srt = 90, xpd = NA)
