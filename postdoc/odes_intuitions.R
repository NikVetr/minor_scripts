# Load the required package
library(deSolve)

# Define the system of ODEs
ode_system <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx <- -2 * x + 3 * y
    dy <- -1 * y + 2 * z
    dz <- -0.5 * x + 1 * y - 1.5 * z
    list(c(dx, dy, dz))
  })
}

# Initial conditions
initial_state <- c(x = 4, y = 1, z = -2)

# Time points
time_points <- seq(0, 10, 1)

# Solve numerically
numerical_solution <- data.frame(ode(y = initial_state, times = time_points, func = ode_system, parms = NULL))

# Solve using matrix exponentiation
solve_matrix_exponentiation <- function(initial_state, time_points) {
  A <- matrix(c(
    -2,  3,  0,
    0, -1,  2,
    -0.5, 1, -1.5
  ), nrow = 3, byrow = TRUE)
  
  exp_solutions <- sapply(time_points, function(t) {
    expm::expm(A * t) %*% initial_state
  })
  
  return(data.frame(time = time_points, t(exp_solutions)))
}

# Use the expm package for matrix exponentiation
if (!require(expm)) install.packages("expm")
library(expm)

matrix_solution <- solve_matrix_exponentiation(initial_state, time_points)

# Print matrix solution
print("Matrix solution:")
print(matrix_solution)

# Print numerical solution
print("Numerical solution:")
print(numerical_solution)

apply(matrix_solution - numerical_solution, 2, mean)
plot(matrix_solution$X1, numerical_solution$x)

#now test nonlinear version

# Nonlinear ODE system
ode_system_nonlinear <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx <- -2 * x + 3 * y * z - 1.1 * x^2
    dy <- -1 * y + 2 * z
    dz <- -0.5 * x + 1 * y - 1.5 * z
    list(c(dx, dy, dz))
  })
}

# Initial conditions
initial_state <- c(x = 4, y = 1, z = -2)

# Time points (from t=0 to t=10)
time_points <- seq(0, 10, 1)

# Solve numerically
nonlinear_solution <- ode(
  y = initial_state,
  times = time_points,
  func = ode_system_nonlinear,
  parms = NULL
)

nonlinear_solution_df <- data.frame(nonlinear_solution)
print(nonlinear_solution_df)

# Plot solution
plot(nonlinear_solution_df$time, nonlinear_solution_df$x, type="b", col="red",
     ylim=c(min(nonlinear_solution_df[,2:4]), max(nonlinear_solution_df[,2:4])),
     xlab="time", ylab="values", main="Nonlinear ODE solution")
points(nonlinear_solution_df$time, nonlinear_solution_df$y, type="b", col="blue")
points(nonlinear_solution_df$time, nonlinear_solution_df$z, type="b", col="green")
legend("topright", legend=c("x", "y", "z"), col=c("red","blue","green"), lty=1)
