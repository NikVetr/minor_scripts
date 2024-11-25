# Extract diagnostics
diagnostics <- fit$sampler_diagnostics()

# Convert the diagnostics to a data frame
diagnostics_df <- as.data.frame(as_draws_df(diagnostics))

# Extract only the columns related to divergent transitions
divergences <- diagnostics_df[, grep("divergent__", colnames(diagnostics_df))]

# Sum across chains to mark whether an iteration is divergent (1 or more divergences)
divergent_flag <- divergences > 0
divinds <- which(divergent_flag)
ndivinds <- setdiff(1:nrow(samps), divinds)

# Extract the posterior samps as a data frame
samps <- as.data.frame(as_draws_df(fit$draws()))
mean_divergent_quantiles <- sort(apply(samps, 2, function(x) mean(sapply(divinds, function(i) mean(x[i] > x[ndivinds])))))

head(mean_divergent_quantiles[order(0.5 - abs(mean_divergent_quantiles - 0.5))], 30)

# Add the divergent information to the samps
samps$divergent <- divergent_flag

param1 <- "sigma_base"
param2 <- "sigma_ovc"

param1 <- "mu_mean_grade"
param2 <- "mu_ovc"

param1 <- "lodds_AB_base"
param2 <- "lodds_AB_mean_grade"

param1 <- "clodds_loh_mean_grade"
param2 <- "sd_clodds_loh_grade"

param1 <- "mu_grade_csum[5]"
param2 <- "mu_grade_csum[4]"

param1 <- "sigma_mean_grade"
param2 <- "sigma_ovc"


# Plot the two parameters with divergent points highlighted
plot(samps[[param1]][!samps$divergent], samps[[param2]][!samps$divergent],
     pch = 19, xlab = param1, ylab = param2,
     main = "Scatter Plot of Parameters with Divergent Iterations")
points(samps[[param1]][samps$divergent], samps[[param2]][samps$divergent],
     col = adjustcolor("red", 0.5),
)
# legend("topright", legend = c("No Divergence", "Divergent"),
#        col = c("black", "red"), pch = 19)
