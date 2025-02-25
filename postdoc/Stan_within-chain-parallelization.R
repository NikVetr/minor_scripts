library(cmdstanr)
library(posterior)
library(microbenchmark)

set.seed(123)

#set benchmarking parameters
niter <- 500
max_nt <- 12
nthreads <- 1:max_nt
n_microb_reps <- 50
max_cores <- 12

#set model parameters
N <- 5E5
true_mu <- 2
true_sigma <- 1.5

#simulate data
y <- rnorm(N, mean = true_mu, sd = true_sigma)

#specify stan model
model_parallel_code <- "
functions {
  real partial_sum_lpdf(array[] real y_subset, int start, int end, real mu, real sigma) {
    return normal_lpdf(y_subset | mu, sigma);
  }
}
data {
  int<lower=1> N;
  int<lower=1> N_cores;
  array[N] real y;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  mu ~ normal(0, 10);
  sigma ~ exponential(1);
  int grainsize = N / N_cores / 4;
  target += reduce_sum(partial_sum_lpdf, y, grainsize, mu, sigma);
}
"

#compile model
mod_parallel <- cmdstan_model(write_stan_file(model_parallel_code), 
                              cpp_options = list(stan_threads = TRUE))

#benchmark
bench_parallel <- data.frame(cbind(do.call(rbind, lapply(nthreads, function(threads) {
  stan_data <- list(N = N, y = y, N_cores = threads)
  parallelize_benchmark <- F
  if(parallelize_benchmark){
    #parallelize?
    n_par_poss <- floor(max_cores / threads)
    n_mr <- ceiling(n_microb_reps / n_par_poss)
    mean_time <- mean(unlist(parallel::mclapply(1:n_mr, function(i){
      out <- microbenchmark(
        mod_parallel$sample(data = stan_data,
                            seed = 123,
                            chains = 1,
                            iter_warmup = niter,
                            iter_sampling = niter,
                            parallel_chains = 1,
                            threads_per_chain = threads),
        times = 1
      )
      return(mean(out$time))  
    }, mc.cores = n_par_poss)))
  } else {
    #or not.
    out <- microbenchmark(
      mod_parallel$sample(data = stan_data,
                          seed = 123,
                          chains = 1,
                          iter_warmup = niter,
                          iter_sampling = niter,
                          parallel_chains = 1,
                          threads_per_chain = threads),
      times = n_microb_reps
    )
    mean_time <- mean(out$time)
  }
  return(c(time = mean_time))
  
})), nthreads = nthreads))
bench_parallel$time <- bench_parallel$time / bench_parallel$time[1] #scale to maximum time (single-core)

##### visualize #####

png(paste0("~/Pictures/Stan_Parallel_Test-n=", N, ".png"), width = 1000, height = 800, pointsize = 24)
par(mar = c(5,5,2,5))
plot(bench_parallel$nthreads, bench_parallel$time, ylim = c(0, max(bench_parallel$time) * 1.05),
     xlab = "number of parallel threads", ylab = "relative time", las = 1, main = paste0("N = ", N))
pusr <- par("usr")
lines(bench_parallel$nthreads, bench_parallel$time)
points(nthreads, bench_parallel$time[1] / nthreads, col = 2, pch = 1)
lines(nthreads, bench_parallel$time[1] / nthreads, lty = 2, col = 2)

relative_number_of_CPU_hours <- bench_parallel$time * nthreads / max(bench_parallel$time * nthreads)
points(nthreads, relative_number_of_CPU_hours, col = 4, pch = 2)
lines(nthreads, relative_number_of_CPU_hours, lty = 2, col = 4)

#you add a thread, and it should decrease the time by some amount, but instead it [WHAT] by [WHAT] amount?
marginal_efficiency <- c(1, diff(bench_parallel$time) / diff(1/nthreads))
scaled_marginal_efficiency <- (marginal_efficiency - min(marginal_efficiency)) / diff(range(marginal_efficiency))
points(nthreads, scaled_marginal_efficiency, col = 6, pch = 3)
lines(nthreads, scaled_marginal_efficiency, lty = 2, col = 6)
axis(4, labels = pretty(marginal_efficiency), 
     at = (pretty(marginal_efficiency) - min(marginal_efficiency)) / diff(range(marginal_efficiency)), 
     las = 1
)
text(label = "relative marginal efficiency", 
     x = pusr[2] + diff(pusr[1:2])/8, 
     y = mean(pusr[3:4]), 
     xpd = NA, srt = 270, col = 6)
relative_sme_bounds <- (0:1 - min(marginal_efficiency)) / diff(range(marginal_efficiency))
abline(h = relative_sme_bounds[1], lty = 3, col = adjustcolor(1, 0.4), lwd = 0.5)
abline(h = relative_sme_bounds[2], lty = 3, col = adjustcolor(1, 0.4), lwd = 0.5)
polygon(x = pusr[c(1,2,2,1,1)],
        y = relative_sme_bounds[c(1,1,2,2,1)], border = NA, col = adjustcolor(6, 0.05))

dev.off()

png(paste0("~/Pictures/Stan_Parallel_Test-legend.png"), width = 600, height = 400)
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1))
legend(x = 0, y = 1, lwd = c(1,2,2,2), pch = c(1,1,2,3), lty = c(1,2,2,2),
       legend = c("observed average relative time",
                  "expected average relative time",
                  "relative # of CPU-hours",
                  "relative marginal efficiency"),
       col = c(adjustcolor(1, 0.5), 2, 5, 6), cex = 1.75)

dev.off()

#arrange in grid
file.paths <- c(paste0("~/Pictures/Stan_Parallel_Test-n=", c(500,5000,5E4,5E5), ".png"),
                "~/Pictures/Stan_Parallel_Test-legend.png", 
                "~/Pictures/Stan_Parallel_Test-model.png")
file.paths <- file.paths[c(1,2,6,3,4,5)]
image_grobs <- lapply(file.paths, function(fp) 
  grid::rasterGrob(png::readPNG(fp), interpolate = TRUE))
png(paste0("~/Pictures/Stan_Parallel_Test-combined.png"), width = 2000, height = 1100)
gridExtra::grid.arrange(grobs = image_grobs, nrow = 2, ncol = 3,
                        top = grid::textGrob(c("2019 Macbook Pro (i9-9980HK, 8C/16T)\n",
                                               "2024 Mac Mini (M4, 10C, 4P/6E)\n")[1], 
                                             gp = grid::gpar(fontsize = 20, fontface = "bold")))
dev.off()
