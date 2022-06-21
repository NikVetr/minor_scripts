a_1 = 0
v0_1 = 100
x0_1 = 0

a_2 = -10
v0_2 = 100
x0_2 = 0
t <- seq(0, v0 / (-a) * 2, length.out = 200)
Xt <- function(x0, t, v0, a){
  x0 + v0*t + a/2*t^2
}
plot(Xt(x0_1, t, v0_1, a_1), Xt(x0_2, t, v0_2, a_2), type = "l")

