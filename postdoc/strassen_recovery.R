# ---------------------------------------
# variable-m, variable-(kL,kR) MCMC (base R)
# ---------------------------------------

# indices: 1:a11, 2:a12, 3:a21, 4:a22, 5:b11, 6:b12, 7:b21, 8:b22
flatten_AB <- function(A, B) {
  v <- numeric(8)
  v[1] <- A[1,1]; v[2] <- A[1,2]; v[3] <- A[2,1]; v[4] <- A[2,2]
  v[5] <- B[1,1]; v[6] <- B[1,2]; v[7] <- B[2,1]; v[8] <- B[2,2]
  v
}

# random state for arbitrary m, kL, kR, q_len
random_state <- function(m = 7, kL = 2, kR = 2, q_len = 6, separable = TRUE) {
  left_pool  <- if (separable) 1:4 else 1:8
  right_pool <- if (separable) 5:8 else 1:8
  list(
    T_idx  = cbind(
      matrix(sample(left_pool,  m*kL,  replace = TRUE), nrow = m, ncol = kL),
      matrix(sample(right_pool, m*kR,  replace = TRUE), nrow = m, ncol = kR)
    ),
    r_sign = matrix(sample(c(-1L,0L,1L), m*(kL+kR), replace = TRUE), nrow = m, ncol = kL+kR),
    q_idx  = matrix(sample(1:m, 4*q_len, replace = TRUE), nrow = 4, ncol = q_len),
    s_sign = matrix(sample(c(-1L,0L,1L), 4*q_len, replace = TRUE), nrow = 4, ncol = q_len),
    kL = kL, kR = kR, separable = separable
  )
}

# exact Strassen seed; works when m == 7
# exact Strassen seed; supports (kL,kR) >= 2 via zero padding
strassen_seed <- function(m = 7, kL = 2, kR = 2, q_len = 6, separable = TRUE) {
  stopifnot(m >= 1, kL >= 1, kR >= 1, q_len >= 1)
  if (m != 7 || kL < 2 || kR < 2) return(NULL)
  if (!separable) warning("For universal correctness, use separable=TRUE (A-left, B-right).")
  
  # zero-fill
  T_idx  <- matrix(0L, nrow = 7, ncol = kL + kR)
  r_sign <- matrix(0L, nrow = 7, ncol = kL + kR)
  q_idx  <- matrix(0L, nrow = 4, ncol = q_len)
  s_sign <- matrix(0L, nrow = 4, ncol = q_len)
  
  # convenience aliases for the 4 slots we actually use per side
  L1 <- 1; L2 <- 2; R1 <- kL + 1; R2 <- kL + 2
  
  # --- Fill M1..M7 (classic Strassen) with direct indexing ---
  # M1 = (a11 + a22) * (b11 + b22)
  T_idx[1, c(L1, L2, R1, R2)]  <- c(1L, 4L, 5L, 8L)
  r_sign[1, c(L1, L2, R1, R2)] <- c(+1L,+1L,+1L,+1L)
  
  # M2 = (a21 + a22) * b11
  T_idx[2, c(L1, L2, R1)]      <- c(3L, 4L, 5L)
  r_sign[2, c(L1, L2, R1)]     <- c(+1L,+1L,+1L)
  
  # M3 = a11 * (b12 - b22)
  T_idx[3, c(L1, R1, R2)]      <- c(1L, 6L, 8L)
  r_sign[3, c(L1, R1, R2)]     <- c(+1L,+1L,-1L)
  
  # M4 = a22 * (b21 - b11)
  T_idx[4, c(L1, R1, R2)]      <- c(4L, 7L, 5L)
  r_sign[4, c(L1, R1, R2)]     <- c(+1L,+1L,-1L)
  
  # M5 = (a11 + a12) * b22
  T_idx[5, c(L1, L2, R1)]      <- c(1L, 2L, 8L)
  r_sign[5, c(L1, L2, R1)]     <- c(+1L,+1L,+1L)
  
  # M6 = (a21 - a11) * (b11 + b12)
  T_idx[6, c(L1, L2, R1, R2)]  <- c(3L, 1L, 5L, 6L)
  r_sign[6, c(L1, L2, R1, R2)] <- c(+1L,-1L,+1L,+1L)
  
  # M7 = (a12 - a22) * (b21 + b22)
  T_idx[7, c(L1, L2, R1, R2)]  <- c(2L, 4L, 7L, 8L)
  r_sign[7, c(L1, L2, R1, R2)] <- c(+1L,-1L,+1L,+1L)
  
  # --- Fill Q’s (pad zeros beyond used columns) ---
  # c11 = M1 + M4 - M5 + M7
  q_idx[1, 1:4]  <- c(1L, 4L, 5L, 7L)
  s_sign[1, 1:4] <- c(+1L,+1L,-1L,+1L)
  
  # c12 = M3 + M5
  q_idx[2, 1:2]  <- c(3L, 5L)
  s_sign[2, 1:2] <- c(+1L,+1L)
  
  # c21 = M2 + M4
  q_idx[3, 1:2]  <- c(2L, 4L)
  s_sign[3, 1:2] <- c(+1L,+1L)
  
  # c22 = M1 - M2 + M3 + M6
  q_idx[4, 1:4]  <- c(1L, 2L, 3L, 6L)
  s_sign[4, 1:4] <- c(+1L,-1L,+1L,+1L)
  
  list(T_idx = T_idx, r_sign = r_sign, q_idx = q_idx, s_sign = s_sign,
       kL = kL, kR = kR, separable = separable)
}

# Evaluate M, Q, loss, log-target for arbitrary m,kL,kR,q_len
# --- FIXED evaluate_state (skip t==0) ---
evaluate_state <- function(state, A, B, C, eps = 1e-12, debug = FALSE) {
  vals <- flatten_AB(A, B)
  m     <- nrow(state$T_idx)
  kL    <- state$kL
  kR    <- state$kR
  q_len <- ncol(state$q_idx)
  
  M <- numeric(m)
  for (i in 1:m) {
    left  <- 0
    right <- 0
    
    # left sum
    for (j in 1:kL) {
      t <- state$T_idx[i, j]
      s <- state$r_sign[i, j]
      if (s != 0L && t >= 1L && t <= length(vals)) {
        left <- left + s * vals[t]
      }
    }
    
    # right sum
    for (j in 1:kR) {
      t <- state$T_idx[i, kL + j]
      s <- state$r_sign[i, kL + j]
      if (s != 0L && t >= 1L && t <= length(vals)) {
        right <- right + s * vals[t]
      }
    }
    
    M[i] <- left * right
    if (debug) cat(sprintf("[eval] M%d: L=%.6g R=%.6g V=%.6g\n", i, left, right, M[i]))
  }
  
  Q <- numeric(4)
  for (k in 1:4) {
    s <- 0
    for (j in 1:q_len) {
      idx <- state$q_idx[k, j]
      sig <- state$s_sign[k, j]
      if (sig != 0L && idx >= 1L && idx <= m) s <- s + sig * M[idx]
    }
    Q[k] <- s
    if (debug) cat(sprintf("[eval] Q%d: %.6g\n", k, Q[k]))
  }
  
  Cvec <- c(C[1,1], C[1,2], C[2,1], C[2,2])
  L <- sum((Q - Cvec)^2)
  logpi <- -log(L + eps)
  list(M = M, Q = Q, loss = L, logpi = logpi)
}

# Symmetric single-parameter proposal respecting (kL,kR) and separability
propose_state <- function(state) {
  new_state <- state
  m  <- nrow(state$T_idx)
  kL <- state$kL; kR <- state$kR
  left_pool  <- if (state$separable) 1:4 else 1:8
  right_pool <- if (state$separable) 5:8 else 1:8
  
  which_block <- sample(c("T", "r", "q", "s"), 1)
  if (which_block == "T") {
    i <- sample(1:m, 1); j <- sample(1:(kL+kR), 1)
    cur <- state$T_idx[i, j]
    pool <- if (j <= kL) left_pool else right_pool
    cand <- sample(setdiff(pool, cur), 1)
    new_state$T_idx[i, j] <- cand
  } else if (which_block == "r") {
    i <- sample(1:m, 1); j <- sample(1:(kL+kR), 1)
    cur <- state$r_sign[i, j]
    cand <- sample(setdiff(c(-1L,0L,1L), cur), 1)
    new_state$r_sign[i, j] <- cand
  } else if (which_block == "q") {
    k <- sample(1:4, 1); j <- sample(1:ncol(state$q_idx), 1)
    cur <- state$q_idx[k, j]
    cand <- sample(setdiff(1:m, cur), 1)
    new_state$q_idx[k, j] <- cand
  } else {
    k <- sample(1:4, 1); j <- sample(1:ncol(state$s_sign), 1)
    cur <- state$s_sign[k, j]
    cand <- sample(setdiff(c(-1L,0L,1L), cur), 1)
    new_state$s_sign[k, j] <- cand
  }
  new_state
}

# helper: count M terms in a state
num_M <- function(state) {
  nrow(state$T_idx)
}

# helper: coerce an init state to have exactly m Ms (truncate or pad)
coerce_state_m <- function(state, m) {
  m0 <- num_M(state)
  kL <- state$kL; kR <- state$kR
  q_len <- ncol(state$q_idx)
  
  if (m0 == m) return(state)
  
  if (m0 > m) {
    # truncate to first m rows
    state$T_idx  <- state$T_idx[1:m, , drop = FALSE]
    state$r_sign <- state$r_sign[1:m, , drop = FALSE]
    
    # any q_idx > m becomes inactive by zeroing its s_sign
    for (k in 1:4) {
      for (j in 1:q_len) {
        if (state$q_idx[k, j] > m) {
          state$s_sign[k, j] <- 0L
        }
      }
    }
  } else {
    # m0 < m: pad with zeros (creates dead Ms that don't affect Q)
    add <- m - m0
    state$T_idx  <- rbind(state$T_idx,  matrix(0L, nrow = add, ncol = kL + kR))
    state$r_sign <- rbind(state$r_sign, matrix(0L, nrow = add, ncol = kL + kR))
    # q_idx can already reference up to m; leave as-is
  }
  state
}

# make a fixed batch of random pairs
make_batch <- function(n, dist = rnorm, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pairs <- vector("list", n)
  for (i in 1:n) {
    A <- matrix(dist(4), 2, 2)
    B <- matrix(dist(4), 2, 2)
    pairs[[i]] <- list(A = A, B = B, C = A %*% B)
  }
  pairs
}

# evaluate a state against a batch (sum of L2 losses)
evaluate_state_batch <- function(state, batch, eps = 1e-12) {
  loss <- rep(NA, length(batch))
  for (i in seq_along(batch)) {
    ab <- batch[[i]]
    ev <- evaluate_state(state, ab$A, ab$B, ab$C, eps = eps, debug = FALSE)
    loss[i] <- ev$loss
  }
  totL <- mean(loss)
  # log-target for MH on the *sum* loss
  list(loss = totL, logpi = -log(totL + eps))
}


# MCMC driver
mcmc_strassen <- function(A, B, m = 7, kL = 2, kR = 2, q_len = 6,
                          n_iter = 20000, init = NULL, separable = TRUE,
                          eps = 1e-12, print_every = 1000, seed = NULL,
                          batch = NULL, pow = 1) {
  if (!is.null(seed)) set.seed(seed)
  
  # init (seed or random), enforce m
  if (is.null(init)) {
    st <- strassen_seed(m = m, kL = kL, kR = kR, q_len = q_len, separable = separable)
    init <- if (is.null(st)) random_state(m, kL, kR, q_len, separable) else st
  } else {
    if (is.null(init$kL)) init$kL <- kL
    if (is.null(init$kR)) init$kR <- kR
    if (is.null(init$separable)) init$separable <- separable
    init <- coerce_state_m(init, m)  # from our previous message
  }
  
  # base matrices if no batch is used
  if (is.null(batch)) {
    C <- A %*% B
    cur <- evaluate_state(init, A, B, C, eps = eps, debug = FALSE)
    get_eval <- function(st) evaluate_state(st, A, B, C, eps = eps, debug = FALSE)
  } else {
    cur <- evaluate_state_batch(init, batch, eps = eps)
    get_eval <- function(st) evaluate_state_batch(st, batch, eps = eps)
  }
  
  state <- init
  best <- list(state = state, loss = cur$loss, logpi = cur$logpi)
  accept <- 0L
  loss_trace <- numeric(n_iter)
  logpi_trace <- numeric(n_iter)
  
  for (t in 1:n_iter) {
    cand_state <- propose_state(state)
    cand <- get_eval(cand_state)
    log_alpha <- (cand$logpi - cur$logpi) * pow
    if (log(runif(1)) < log_alpha) {
      state <- cand_state
      cur <- cand
      accept <- accept + 1L
      if (cur$loss < best$loss) {
        best$state <- state; best$loss <- cur$loss; best$logpi <- cur$logpi
      }
    }
    loss_trace[t]  <- cur$loss
    logpi_trace[t] <- cur$logpi
    if (t %% print_every == 0) {
      cat(sprintf("[iter %d] loss=%.3e  logpi=%.3f  acc_rate=%.3f\n",
                  t, cur$loss, cur$logpi, accept/t))
    }
  }
  
  list(best_state = best$state, best_loss = best$loss, best_logpi = best$logpi,
       accept_rate = accept / n_iter, loss_trace = loss_trace, logpi_trace = logpi_trace,
       final_state = state, final_eval = cur)
}


# Pretty-printer for arbitrary kL,kR
describe_state <- function(state) {
  nm <- function(idx) c("a11","a12","a21","a22","b11","b12","b21","b22")[idx]
  m  <- nrow(state$T_idx); kL <- state$kL; kR <- state$kR
  q_len <- ncol(state$q_idx)
  m_lines <- character(m)
  for (i in 1:m) {
    Lparts <- character(0); Rparts <- character(0)
    for (j in 1:kL) {
      s <- state$r_sign[i, j]; t <- state$T_idx[i, j]
      if (s != 0L) Lparts <- c(Lparts, sprintf("%+d*%s", s, nm(t)))
    }
    for (j in 1:kR) {
      s <- state$r_sign[i, kL + j]; t <- state$T_idx[i, kL + j]
      if (s != 0L) Rparts <- c(Rparts, sprintf("%+d*%s", s, nm(t)))
    }
    Lexpr <- if (length(Lparts)) gsub("^\\+", "", paste(Lparts, collapse=" ")) else "0"
    Rexpr <- if (length(Rparts)) gsub("^\\+", "", paste(Rparts, collapse=" ")) else "0"
    m_lines[i] <- sprintf("M%d = (%s) * (%s)", i, Lexpr, Rexpr)
  }
  q_lines <- character(4)
  for (k in 1:4) {
    parts <- character(0)
    for (j in 1:q_len) {
      sig <- state$s_sign[k, j]; idx <- state$q_idx[k, j]
      if (sig != 0L && idx >= 1 && idx <= m) parts <- c(parts, sprintf("%+d*M%d", sig, idx))
    }
    expr <- if (length(parts)) gsub("^\\+", "", paste(parts, collapse = " ")) else "0"
    q_lines[k] <- sprintf("Q%d = %s", k, expr)
  }
  list(M = m_lines, Q = q_lines)
}

# turn a learned state into a callable function f(A,B) -> 2x2 matrix via Q's
make_algorithm <- function(state, eps = 0) {
  kL <- state$kL; kR <- state$kR
  force(state)
  function(A, B) {
    # compute M_i
    vals <- c(A[1,1],A[1,2],A[2,1],A[2,2], B[1,1],B[1,2],B[2,1],B[2,2])
    m <- nrow(state$T_idx)
    M <- numeric(m)
    for (i in 1:m) {
      left <- 0; right <- 0
      for (j in 1:kL) {
        t <- state$T_idx[i, j]; s <- state$r_sign[i, j]
        if (s != 0L && t >= 1L && t <= 8L) left <- left + s * vals[t]
      }
      for (j in 1:kR) {
        t <- state$T_idx[i, kL + j]; s <- state$r_sign[i, kL + j]
        if (s != 0L && t >= 1L && t <= 8L) right <- right + s * vals[t]
      }
      M[i] <- left * right
    }
    # assemble Q -> C
    q_len <- ncol(state$q_idx)
    Q <- numeric(4)
    for (k in 1:4) {
      s <- 0
      for (j in 1:q_len) {
        idx <- state$q_idx[k, j]; sig <- state$s_sign[k, j]
        if (sig != 0L && idx >= 1L && idx <= m) s <- s + sig * M[idx]
      }
      Q[k] <- s
    }
    C_hat <- matrix(c(Q[1], Q[2], Q[3], Q[4]), nrow = 2, byrow = TRUE)
    if (eps != 0) C_hat <- C_hat + 0 # (placeholder: no reg. used)
    C_hat
  }
}

# pretty equations (you also have describe_state already)
print_equations <- function(state) {
  desc <- describe_state(state)
  cat("M terms:\n", paste0(desc$M, collapse = "\n"), "\n\n", sep = "")
  cat("Q terms:\n", paste0(desc$Q, collapse = "\n"), "\n", sep = "")
}

# example:


set.seed(1)
A <- matrix(rnorm(4), 2, 2)
B <- matrix(rnorm(4), 2, 2)

# 7-term seed (works for any kL,kR >= 2)
stp <- strassen_seed(m = 7, kL = 2, kR = 2, q_len = 4, separable = TRUE)
init_rand <- random_state(m = 7, kL = 2, kR = 2, q_len = 4, separable = TRUE)
m <- 7
n <- 4
batch <- make_batch(n = 4, seed = 42)
out <- mcmc_strassen(m = m, kL = 2, kR = 2, q_len = 4, separable = TRUE,
  init = init_rand, n_iter = 5E5, print_every = 5E4, batch = batch, pow = 2
)
cat(paste0("m=", m, " best loss:", out$best_loss, "\n"))
print_equations(out$best_state)
