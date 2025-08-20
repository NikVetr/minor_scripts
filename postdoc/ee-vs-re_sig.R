#simulation code
n <- 500
x <- matrix(rnorm(n*2), nrow = n)
r <- 0.7
R <- diag(2) + r - diag(2) * r
x <- x %*% chol(R)
invlogit <- function(x) exp(x) / (1 + exp(x))
p <- invlogit(abs(x) * 4 - 8)
xs <- cbind(rbinom(n, 1, p[,1]), rbinom(n, 1, p[,2]))
dx <- x[,1] - x[,2]
dp <- invlogit(abs(dx))
ds <- rbinom(n, 1, dp)

#helpers
maj <- factor(xs[,1] + 2*xs[,2],          # 0 = neither, 1 = x1 only,
              levels = c(1,3,0,2),        #   2 = x2 only, 3 = both
              labels = c("EE only", "both", "neither", "RE only"))

sub <- factor(c("x1","diff","x2"),        # minor rows inside each band
              levels = c("x1","diff","x2"))

#graphical params
band_ht  <- 3                             # three sub-rows per band
y_maj    <- as.numeric(maj)               # 1–4
y_x1     <- (4 - y_maj) * band_ht + 2
y_diff   <- (4 - y_maj) * band_ht + 1
y_x2     <- (4 - y_maj) * band_ht
alphas <- c(0.5, 0.1)
pt_cex <- 1

#colors for classes
base_col <- c("firebrick", "purple3", "grey50", "royalblue")[y_maj]

#opacity according to significance of difference
col_pt <- base_col                   # start with the per-point hue vector
col_pt[ ds == 1 ] <- adjustcolor(col_pt[ ds == 1 ], alpha.f = alphas[1])  # Δ sig
col_pt[ ds == 0 ] <- adjustcolor(col_pt[ ds == 0 ], alpha.f = alphas[2])  # Δ not

#set up plotting window
xr <- range(c(x[,1], x[,2], dx))# + c(-1,1) * diff(xr) / 10
par(mar = c(5,8,2,2))
plot(NA, xlim = xr, ylim = c(0, 12),
     xlab = "value", ylab = "",
     yaxt = "n", bty = "n",
     main = "")
abline(v = 0, lty = 2, lwd = 2, col = adjustcolor(1, 0.25))

#guide lines for bands
siglab_x <- xr[1] - diff(xr)/2.5
tlab_x <- xr[1] - diff(xr)/5
abline(h = c(3,6,9) - 0.5, lty = 3, col = "grey80")
text(x = siglab_x, 1 + 0:3*3, labels = rev(levels(maj)), pos = 4, xpd = NA)
text(x = siglab_x, y = 12, labels = "p < 0.05", pos = 4, font = 2, col = 1, cex = 1, xpd = NA)
text(x = tlab_x, y = 0:11, labels = c("EE log2FC", "diff log2FC", "RE log2FC"), 
     pos = 4, font = 1, xpd = NA, cex = 0.75)


#plot points and connecting lines
segments(x[,1], y_x1, dx, y_diff, col = col_pt, lwd = .8)
segments(dx, y_diff,  x[,2], y_x2,  col = col_pt, lwd = .8)
points(x[,1], y_x1, col = col_pt, pch = 16, cex = pt_cex)
points(dx,    y_diff, col = col_pt, pch = 16, cex = pt_cex)
points(x[,2], y_x2, col = col_pt, pch = 16, cex = pt_cex)

#add a legend
pusr <- par("usr")
legend(x = pusr[2] - diff(pusr[1:2]) / 5, y = pusr[4] + diff(pusr[3:4])/12, inset = .02, bty = "n", title = "Δ significant",
       legend = c("yes","no"), lwd = c(1,1),
       pch = 16, col = sapply(alphas, function(af) adjustcolor("black", alpha.f = af)), xpd = NA)
