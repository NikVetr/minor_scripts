functions {
  // ---------------------------------------------------------
  //  sample_intersection: already defined in your code
  //  (Accepts vector Lambda, real y in [0,1], unit_vector theta,
  //   real sign_flip in {-1,1}, returns a vector on the ellipse-sphere intersection)
  // ---------------------------------------------------------
  vector sample_intersection(vector Lambda, real y, vector theta, real sign_flip) {
    int p = num_elements(Lambda);
    int k = p - 2;  // should match num_elements(theta)

    // 1) Identify indices of min(Lambda) and max(Lambda)
    int i_min = 1;
    int i_max = 1;
    {
      real L_min = Lambda[1];
      real L_max = Lambda[1];
      for (i in 2:p) {
        if (Lambda[i] < L_min) {
          L_min = Lambda[i];
          i_min = i;
        }
        if (Lambda[i] > L_max) {
          L_max = Lambda[i];
          i_max = i;
        }
      }
    }

    // 2) Collect the remaining (p-2) indices
    array[p-2] int remain_index;
    {
      int idx = 1;
      for (i in 1:p) {
        if ((i != i_min) && (i != i_max)) {
          remain_index[idx] = i;
          idx += 1;
        }
      }
    }

    // 3) Extract L_k (remaining) and L_drop (dropped)
    vector[p-2] L_k;
    for (i in 1:(p-2)) {
      L_k[i] = Lambda[ remain_index[i] ];
    }
    vector[2] L_drop;
    L_drop[1] = Lambda[i_min];
    L_drop[2] = Lambda[i_max];

    // 4) Compute scaling factor s
    real sts = dot_self(theta); // sum(theta^2)
    real stes = 0;             // sum(L_k * theta^2)
    for (i in 1:(p-2)) {
      stes += L_k[i] * theta[i]^2;
    }

    real sa  = inv_sqrt(sts);
    real sb1 = inv_sqrt(stes);

    // For each dropped L, compute candidate sb2^2 = (L - 1) / (L * sts - stes)
    vector[2] sb2_sq;
    sb2_sq[1] = (L_drop[1] - 1) / (L_drop[1] * sts - stes);
    sb2_sq[2] = (L_drop[2] - 1) / (L_drop[2] * sts - stes);

    real sb2_1 = (sb2_sq[1] > 0) ? sqrt(sb2_sq[1]) : positive_infinity();
    real sb2_2 = (sb2_sq[2] > 0) ? sqrt(sb2_sq[2]) : positive_infinity();

    // s is the min of sa, sb1, sb2_1, sb2_2
    real s = fmin(sa, fmin(sb1, fmin(sb2_1, sb2_2)));

    // 5) Scale theta by (y * s)
    vector[p-2] x_k = theta * (y * s);

    // 6) Solve the 2x2 linear system for the dropped coordinates
    real sum_xk_sq = dot_self(x_k);
    real sum_Lk_xk_sq = 0;
    for (i in 1:(p-2)) {
      sum_Lk_xk_sq += L_k[i] * x_k[i]^2;
    }

    matrix[2,2] A;
    vector[2] b;
    A[1,1] = 1;
    A[1,2] = 1;
    A[2,1] = L_drop[1];
    A[2,2] = L_drop[2];
    b[1]   = 1 - sum_xk_sq;
    b[2]   = 1 - sum_Lk_xk_sq;

    vector[2] x_d = mdivide_left(A, b);

    // 7) Insert these into the final p-vector
    vector[p] x = rep_vector(0.0, p);

    real x_d1 = sqrt(x_d[1]);
    real x_d2 = sqrt(x_d[2]) * sign_flip;

    x[i_min] = x_d1;
    x[i_max] = x_d2;
    for (i in 1:(p-2)) {
      x[ remain_index[i] ] = x_k[i];
    }

    return x;
  }

  // -------------------------------------------------------------------
  // intersection_constraint:
  //
  // Translated from your R code.  
  // Takes:
  //   - L: vector of length p
  //   - V_curr: a p x pcurr matrix (the columns are existing eigenvectors)
  //   - B_curr: a p x (p - pcurr) matrix for the orthonormal complement 
  //     (or pass a zero-column matrix if pcurr=0 to let it build from identity)
  //   - y in [0,1], pre-sampled from Beta^(1/(prem-2)) if needed
  //   - theta: a unit_vector of dimension (prem-2)
  //   - sign_flip in {-1, 1}
  // Returns:
  //   A p x prem matrix whose:
  //       col 1 = the newly sampled vector (new_vec)
  //       cols 2..prem = the updated basis B_new
  // -------------------------------------------------------------------
  matrix intersection_constraint(
      vector L, 
      matrix V_curr,
      matrix B_curr,
      real y,
      vector theta,
      real sign_flip
  ) {
    int p = num_elements(L);
    int pcurr = cols(V_curr);
    int prem = p - pcurr;

    // 1) Construct a suitable B
    //    If pcurr=0, replicate "B = diag(p)" logic. Otherwise use B_curr.
    matrix[p, prem] B;
    if (pcurr == 0) {
      // Identity is p x p, we only need its first prem columns
      matrix[p, p] I = diag_matrix(rep_vector(1.0, p));
      B = I[, 1:prem];
    } else {
      B = B_curr;   // assumed dimension p x prem
    }

    // 2) Q = B' * diag(L) * B, dimension (prem x prem)
    matrix[p, p] diagL = diag_matrix(L);
    matrix[prem, prem] Q = (B' * diagL) * B;

    // 3) Eigen-decomposition of Q
    //    Stan has eigen_sym() for symmetric matrices
    vector[prem] QD = eigenvalues_sym(Q);
    matrix[prem, prem] QV = eigenvectors_sym(Q);
    // The order of QD is ascending if Q is positive semidefinite. 
    // We'll just use them directly; the "sample_intersection" doesn't require sorting.

    // 4) Sample a point on the intersection for the prem-dim ellipsoid QD
    //    Note: "theta" should be unit_vector[prem-2], "y" in [0,1].
    vector[prem] evec = sample_intersection(QD, y, theta, sign_flip);

    // 5) Transform evec back to original space
    vector[prem] Qvec = QV * evec;   // still prem-dim
    vector[p]    new_vec = B * Qvec; // now p-dim

    // 6) Rank-1 update for B_new
    vector[prem] alpha = B' * new_vec;          // alpha in R^prem
    matrix[p, prem] B_new = B - new_vec * alpha'; // dimension p x prem

    // 7) Re-orthonormalize B_new, keep rank (prem-1)
    //    Stan has qr_thin_Q() which returns a p x prem matrix with orthonormal columns
    matrix[p, prem] Q_temp = qr_thin_Q(B_new);
    // The rank is effectively (prem-1), so the last column is likely near 0 or ill-defined.
    // We'll keep columns [1 : (prem-1)] for the updated basis in R^p:
    matrix[p, prem-1] B_trim = Q_temp[, 1:(prem-1)];

    // 8) Return both "new_vec" and "B_new" in a single matrix:
    //    - col 1 = new_vec
    //    - cols 2..prem = the newly orthonormalized basis
    matrix[p, prem] output;
    output[, 1] = new_vec;
    output[, 2:prem] = B_trim; // fill columns 2..prem

    return output;
  }
}

data {
  // Example data for a minimal run. 
  // In practice, you'd supply real inputs from R or another source.
  int<lower=3> p;           // total dimension
  int<lower=0> pcurr;       // number of existing columns in V_curr
  vector[p] L;              // hyperellipsoid coefficients
  matrix[p, pcurr] V_curr;  // existing basis columns
  matrix[p, (p - pcurr)] B_curr; // orthonormal complement basis (or identity slices if pcurr=0)
  real<lower=0, upper=1> y; // the Beta^(1/(prem-2)) sample in [0,1]
  real sign_flip;           // in {-1, 1}
}

parameters {
  // We want a direction on the (prem-2)-sphere:
  // "prem" = p - pcurr
  unit_vector[(p - pcurr) - 2] theta;
}

model {
  // no priors / no likelihood for a test
}

generated quantities {
  // We'll just call intersection_constraint() once to test it
  matrix[p, (p - pcurr)] out = intersection_constraint(L, V_curr, B_curr, y, theta, sign_flip);
  // The user can extract the new vector and the updated B_new from out
  // in their R code or in another Stan block if needed:
  // new_vec = out[,1], B_new = out[, 2:(p - pcurr)]
}
