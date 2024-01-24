n <- "Johnnie Walker Red Label
Johnnie Walker Black Label
Johnnie Walker Double Black
Johnnie Walker Island Green
Johnnie Walker Green Label
Johnnie Walker Gold Label Reserve
Johnnie Walker Aged 18 Years
Johnnie Walker Blue Label
Johnnie Walker Platinum Label
Johnnie Walker King George
Johnnie Walker Swing
Johnnie Walker XR 21
Johnnie Walker Explorers' Club Collection
Johnnie Walker White Walker
Johnnie Walker Odyssey"

n <- strsplit(n, "\n")[[1]]
n <- gsub("Johnnie Walker ", "", n)


p <- "$21 
$31 
$43 
$47 
$57 
$85 
$200 
$312 
$105 
$307 
$60 
$155 
$165 
$46 
$1,725"

p <- strsplit(p, " ")[[1]]
p <- log10(as.numeric(gsub("\\$|\n|,", "", p)))

par(mar = c(9,5,3,2))
plot(diff(sort(p)), type = "l")
plot(sort(p), type = "l", lwd = 2, main = "Johnnie Walkers' Bottle Prices",
     ylab = latex2exp::TeX("log$_{10}$(cost / 750ml) (2022 USD)"), xaxt = "n", xlab = "")
i <- 1:length(p)
axis(1, at = i, labels = rep("", length(i)))
text(i + 0.2, y = par("usr")[3] - diff(par("usr")[3:4])/25, labels = n, xpd = NA, srt = 45, pos = 2)
abline(lm(sort(p) ~ i), col = 2, lwd = 2, lty = 2)
abline(h = 3:6/2, lty = 3, lwd = 0.5)
abline(v = i, lty = 3, lwd = 0.5)
legend(x = "topleft", legend = c("observed prices", "linear fit (ols)"), lwd = 2, lty = c(1,2), col = 1:2)
