# dark orange and dark indigo hex codes
col1 <- "#D55E00"   # orange
col2 <- "#4B0082"   # indigo

# set up a 3-panel layout
par(mfrow = c(1, 3), mar = c(3, 1, 3, 1) + 0.1, xpd = NA)

##### panel 1 – location ######################################################
curve(dnorm(x, -1, 1),
      from = -4, to = 4, ylim = c(0, 0.5),
      lwd = 2, col = col1, axes = FALSE, xlab = "", ylab = "",
      main = "location", cex.main = 2)
curve(dnorm(x,  1, 1),
      add = TRUE, lwd = 2, col = col2)

# dashed verticals at the means
segments(-1, 0, -1, dnorm(-1, -1, 1) * 1.05, lty = 2, col = col1)
segments( 1, 0,  1, dnorm( 1,  1, 1) * 1.05, lty = 2, col = col2)

# mu labels
text(-1, dnorm(-1, -1, 1) *1.05, expression(mu[1]), pos = 3, col = col1, cex = 2)
text( 1, dnorm( 1,  1, 1) *1.05, expression(mu[2]), pos = 3, col = col2, cex = 2)

# double-headed arrow between the means
arrows(-1, dnorm(-1, -1, 1) * 1.25, 1, 
       dnorm(-1, -1, 1) *1.25, code = 3, angle = 90, length = 0.08)

##### panel 2 – scale #########################################################
curve(dnorm(x, -2, 1),
      from = -5, to = 5, ylim = c(0, 0.5),
      lwd = 2, col = col1, axes = FALSE, xlab = "", ylab = "",
      main = "scale", cex.main = 2)
curve(dnorm(x,  1, sqrt(2)),
      add = TRUE, lwd = 2, col = col2)

# sigma arrows and labels
y1 <- dnorm(-1, -1, 1) + 0.02
arrows(- 2 - 1, y1, -2 + 1, y1, code = 3, angle = 45, length = 0.07, col = col1)
text(-2, y1, expression(sigma[1]), pos = 3, col = col1, cex = 2)

y2 <- dnorm(1, 1, sqrt(2)) + 0.02
arrows(1 - sqrt(2), y2, 1 + sqrt(2), y2, code = 3, angle = 45, length = 0.07, col = col2)
text(1, y2, expression(sigma[2]), pos = 3, col = col2, cex = 2)

##### panel 3 – modularity (fixed) ############################################
# grid for bivariate normals
x <- seq(-5, 5, length.out = 120)
y <- seq(-5, 5, length.out = 120)

# helper: unit-variance bivariate normal pdf with correlation rho
bvn <- function(x, y, mx, my, rho) {
  z  <- ( (x - mx)^2 - 2 * rho * (x - mx) * (y - my) + (y - my)^2 )
  exp( -z / (2 * (1 - rho^2)) ) / (2 * pi * sqrt(1 - rho^2))
}

# build z-matrices with outer(): rows = x, cols = y  →  correct for contour()
z1 <- outer(x, y, function(x, y) bvn(x, y, -3,  1, 0.1))  # r1 = 0.1
z2 <- outer(x, y, function(x, y) bvn(x, y,  3, -1, 0.8))  # r2 = 0.8

# draw contours

qs <- c(0.0001, 1:9/10)
contour(x, y, z1,
        levels = quantile(z1, qs),
        drawlabels = FALSE, axes = FALSE, xlab = "", ylab = "",
        main = "modularity", col = col1, lwd = 2, 
        xlim = c(-5,5), ylim = c(-5,5), cex.main = 2)
contour(x, y, z2,
        levels = quantile(z2, qs),
        drawlabels = FALSE, axes = FALSE, add = TRUE,
        col = col2, lwd = 2)

# r labels near the means
text(-3,  1, expression(r[1]), col = col1, cex = 2)
text( 3, -1, expression(r[2]),  col = col2, cex = 2)
