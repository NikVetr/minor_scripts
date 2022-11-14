d <- read.table("G:\\Downloads\\usa_00001.dat.gz", header=F)
d$year <- as.numeric(substring(d$V1, 1, 4))
d$age <- as.numeric(substring(d$V1, 5, nchar(d$V1)))

con117 <- readLines("G:\\Downloads\\117_congress.txt", warn = F)
con117 <- as.data.frame(do.call(rbind, strsplit(con117, "\t")))
mean(as.numeric(con117$V5) < 40.1)

d <- split(d$age, d$year)
years <- as.numeric(names(d))
prop <- sapply(d, function(x) sum(x > 70.1) / (sum(x > 29.9) * 100 / 535 + sum(x > 24.9) * 435 / 535))
prop2 <- sapply(d, function(x) sum(x > 70.1) / (sum(x > 40.1)))
plot(years, prop, type = "l", lwd = 2, col = "orange",
  ylim = c(0, max(c(prop, prop2))), xlim = c(1800, max(years)),
  ylab = "proportion of US population > 70", xlab = "year")
lines(years, prop2, lwd = 2, col = "blue")
legend("topleft", col = c("blue", "orange"), lwd = 2, legend = c("out of people > 40", "out of weighted people > 25 & > 30"))
