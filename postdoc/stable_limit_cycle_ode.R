library(deSolve)

# ============================================================================
# 1) USER PARAMETERS
# ============================================================================
N <- 7          # Total number of states
K <- 3          # Size of the forced-oscillator subnetwork (must be <= N)

# "Repressilator" parameters (used if K == 3)
alpha  <- 25    # Production in the repressor terms
gamma  <- 1     # Decay for the repressilator
hill_n <- 3     # Hill exponent for the forced oscillator
hill_K <- 1     # Hill saturation constant for the forced oscillator

# If K > 3, we do a "negative ring" among 1..K. 
# We'll store ring_k < 0 for each link in the ring.

# Random adjacency for states (K+1..N):
num_extra_params <- 20   # number of random edges
mean_k_pos <- 1.0        # positive edge magnitude
mean_k_neg <- -1.0       # negative edge magnitude
sd_k <- 0.5

# Production/decay for states (K+1..N):
decay_min <- 0.0
decay_max <- 0.05
prod_min  <- 0.5
prod_max  <- 1.5

# Simulation times:
tmax  <- 1e4
dt    <- 0.1
times <- seq(0, tmax, by = dt)

# ============================================================================
# 2) BUILD FORCED OSCILLATOR in STATES 1..K
# ============================================================================
# We'll define two possible oscillators:
#  (1) K == 3 => Repressilator among states {1,2,3}
#  (2) K > 3  => Negative ring among states {1..K}, i.e. i <- i-1 with negative k

make_forced_oscillator <- function(K, alpha, gamma, hill_n, hill_K) {
  # We'll return a function that computes dY[1..K]/dt 
  # given (Y[1..K]) for the forced oscillator part.
  
  if (K == 3) {
    # Repressilator:
    #   dy1 = alpha/(1 + (y3/K)^n) - gamma*y1
    #   dy2 = alpha/(1 + (y1/K)^n) - gamma*y2
    #   dy3 = alpha/(1 + (y2/K)^n) - gamma*y3
    # We'll store them in a closure:
    force_3state <- function(Yvec) {
      y1 <- Yvec[1]
      y2 <- Yvec[2]
      y3 <- Yvec[3]
      dy <- numeric(3)
      dy[1] <- alpha / (1 + (y3 / hill_K)^hill_n) - gamma * y1
      dy[2] <- alpha / (1 + (y1 / hill_K)^hill_n) - gamma * y2
      dy[3] <- alpha / (1 + (y2 / hill_K)^hill_n) - gamma * y3
      dy
    }
    return(force_3state)
    
  } else {
    # Negative ring: 1 <- 2 <- 3 <- ... <- K <- 1
    # We'll do k = -5 for each link, or you could parameterize it further.
    ring_k <- -5
    force_ring <- function(Yvec) {
      # Yvec has length K, states 1..K
      # for i in 1..K:
      #   dYi = ring_k * Hill( Y_{i-1} ) - gamma * Y_i
      # where i-1 is mod K. 
      # We'll do the same Hill function approach:
      #   ring term = k*( Y_{i-1}^n / (hill_K + Y_{i-1}^n ))
      # and also we might add a small production if needed. 
      # But let's see if ring alone suffices to sustain oscillations.
      # We'll also subtract gamma * Y_i as a linear decay.
      # 
      # Because it's tricky to guarantee all negative feedback rings will 
      # produce stable limit cycles without some parameter tinkering,
      # you might want to increase ring_k in magnitude or add baseline production.
      
      dy <- numeric(K)
      for (i in seq_len(K)) {
        # index of previous node in ring: 
        # if i=1 => prev = K, else prev = i-1
        prev <- if (i==1) K else (i-1)
        prev_val <- Yvec[prev]
        
        # negative Hill
        feedback <- ring_k * (prev_val^hill_n / (hill_K + prev_val^hill_n))
        dy[i] <- feedback - gamma*Yvec[i]
      }
      dy
    }
    return(force_ring)
  }
}

osc_fun <- make_forced_oscillator(K, alpha, gamma, hill_n, hill_K)

# ============================================================================
# 3) BUILD RANDOM NETWORK for STATES (K+1..N)
# ============================================================================
if (K < N) {
  # Indices for the "random block" = (K+1..N)
  block_indices <- (K+1):N
  
  # (A) Self-loops
  base_dependencies <- cbind(
    influenced = block_indices,
    influencer = block_indices
  )
  k_self <- rnorm(length(block_indices), 0.5, 0.3)   # you can adjust
  
  # (B) All pairs in block_indices, excluding self-loops
  all_pairs_block <- as.matrix(expand.grid(block_indices, block_indices))
  colnames(all_pairs_block) <- c("influenced","influencer")
  self_mask <- all_pairs_block[,1] == all_pairs_block[,2]
  possible_dependencies_block <- all_pairs_block[!self_mask,]
  
  if (num_extra_params > nrow(possible_dependencies_block)) {
    stop("Not enough possible edges among states in block.")
  }
  
  # (C) Sample random edges
  extra_dependencies_block <- possible_dependencies_block[
    sample(seq_len(nrow(possible_dependencies_block)), num_extra_params),
  ]
  
  # Force half positive, half negative
  half <- num_extra_params/2
  k_pos <- rnorm(ceiling(half), mean_k_pos, sd_k)
  k_neg <- rnorm(floor(half), mean_k_neg, sd_k)
  k_extra <- sample(c(k_pos, k_neg), num_extra_params)
  
  block_random <- rbind(
    cbind(base_dependencies, k = k_self),
    cbind(extra_dependencies_block, k = k_extra)
  )
  
  # (D) Production & decay in that block
  decay_block <- runif(length(block_indices), decay_min, decay_max)
  prod_block  <- runif(length(block_indices), prod_min, prod_max)
} else {
  # If K == N, then no random block at all.
  block_random <- matrix(nrow=0, ncol=3)
  decay_block <- numeric(0)
  prod_block  <- numeric(0)
}

# ============================================================================
# 4) DEFINE THE ODE
# ============================================================================
ode_system <- function(t, state_vec, parms) {
  N <- parms$N
  K <- parms$K
  
  # (A) forced oscillator for 1..K
  if (K > 0) {
    forced_dy <- parms$forced_fun( state_vec[1:K] )  # get dY[1..K]
  } else {
    forced_dy <- numeric(0)
  }
  
  # (B) adjacency-based for (K+1..N)
  # we have "block_random" with columns: influenced, influencer, k
  # plus separate production & decay for each of states K+1..N
  block_dy <- numeric(N - K)
  
  # 1) baseline production
  for (i in seq_len(N-K)) {
    block_dy[i] <- block_dy[i] + parms$prod_block[i]
  }
  
  # 2) hill interactions
  for (row_i in seq_len(nrow(parms$block_random))) {
    infl <- parms$block_random[row_i,"influenced"]
    by   <- parms$block_random[row_i,"influencer"]
    kval <- parms$block_random[row_i,"k"]
    yby  <- state_vec[by]
    
    # index in block_dy is infl - K
    idx_infl <- infl - K
    
    # Hill function
    # for simplicity, we use the same hill_n, hill_K as the forced oscillator.
    # you can parameterize them differently if you want
    h_n  <- parms$hill_n
    h_K  <- parms$hill_K
    block_dy[idx_infl] <- block_dy[idx_infl] + 
      kval * (yby^h_n / (h_K + yby^h_n))
  }
  
  # 3) decay
  for (i in seq_len(N-K)) {
    block_dy[i] <- block_dy[i] - parms$decay_block[i]*state_vec[K + i]
  }
  
  # (C) combine
  dy <- numeric(N)
  if (K>0) {
    dy[1:K] <- forced_dy
  }
  if (K < N) {
    dy[(K+1):N] <- block_dy
  }
  
  list(dy)
}

# ============================================================================
# 5) PACKAGE EVERYTHING & SOLVE
# ============================================================================
params_list <- list(
  N = N,
  K = K,
  forced_fun = osc_fun,          # subnetwork function for states 1..K
  block_random = block_random,   # adjacency among states K+1..N
  decay_block = decay_block,
  prod_block  = prod_block,
  hill_n      = hill_n,  # used for the random block
  hill_K      = hill_K
)

# initial conditions
y0 <- runif(N, 0.1, 1.0)

tmax <- 1E3
dt <- 0.1
times <- seq(0, tmax, by=dt)

out <- ode(
  y = y0, 
  times = times, 
  func = ode_system, 
  parms = params_list,
  method = "bdf",
  rtol = 1e-3, 
  atol = 1e-4,
  maxsteps = 1e4
)

# Plot
matplot(out[,"time"], out[, -1], type="l", lty=1, xlab="time", ylab="states",
        main=paste(N, "State System, with a Forced Oscillator of size", K))

