#archimedean spiral
n <- 5E3
max_t <- 20 * pi - pi/2
t <- seq(0, max_t, length.out = n)
x <- sin(t) * t
y <- cos(t) * t
plot(x, y, type = "l", asp = 1)
points(tail(x,1), tail(y,1), col = 1, pch = 17)
box(lwd = 2)

#two circles
n <- 1E3
t <- seq(0, 2 * pi, length.out = n)
x <- sin(t)
y <- cos(t)
plot(x, y, type = "l", asp = 1)
# lines(x/2 + 1/5, y/2 + 2/5)
lines(x/2 + 1/2, y/2 + 2/5)
box(lwd = 2)
