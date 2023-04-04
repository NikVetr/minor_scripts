running_inds <- orig_running_inds
plot(running_inds, type = "l")
lines(filtered_inds, col = 2)

f <- function(x) (x)^(1/3)
finv <- function(x) (x)^3

#transform running inds
running_inds$x <- f(running_inds$x)
running_inds$y <- f(running_inds$y)

#transform and extract magnitudes
wavelet_filter <- c("haar", "la8", "d4", "d6", "mb4")[5]
n_levels <- max(1, floor(log2(nrow(running_inds)) / 2))

waveout.1 <- waveslim::modwt(running_inds[,1], wf = wavelet_filter, n.levels = n_levels)
waveout.2 <- waveslim::modwt(running_inds[,2], wf = wavelet_filter, n.levels = n_levels)

# Get the detail coefficients from the first level
detail_coeffs_1 <- waveout.1$d1
detail_coeffs_2 <- waveout.2$d1

# Estimate noise_sigma using the Median Absolute Deviation (MAD)
noise_sigma_1 <- sd(waveout.1$d1)
noise_sigma_2 <- sd(waveout.2$d1)

# Set the threshold level using a multiplier and the noise standard deviation
threshold_multiplier <- 3
threshold_1 <- threshold_multiplier * noise_sigma_1
threshold_2 <- threshold_multiplier * noise_sigma_2

# Apply soft thresholding to the detail coefficients
thresholded_coeffs_1 <- waveout.1
thresholded_coeffs_2 <- waveout.2

for (i in 1:n_levels) {
  level_name <- paste0("d", i)
  thresholded_coeffs_1[[level_name]] <- sign(thresholded_coeffs_1[[level_name]]) * 
    (abs(thresholded_coeffs_1[[level_name]]) - threshold_1) * (abs(thresholded_coeffs_1[[level_name]]) > threshold_1)
  thresholded_coeffs_2[[level_name]] <- sign(thresholded_coeffs_2[[level_name]]) * 
    (abs(thresholded_coeffs_2[[level_name]]) - threshold_2) * (abs(thresholded_coeffs_2[[level_name]]) > threshold_2)
}

denoised_signal_1 <- waveslim::imodwt(thresholded_coeffs_1)
denoised_signal_2 <- waveslim::imodwt(thresholded_coeffs_2)

# Combine the denoised signals into a 2D time series
filtered_inds <- data.frame(x = denoised_signal_1, y = denoised_signal_2)
filtered_inds <- data.frame(x = finv(filtered_inds$x), y = finv(filtered_inds$y))
