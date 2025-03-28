functions {
   /*
    construct_corrmat Arguments:
      - L: a length-p vector of hyperellipsoid coefficients.
      - y_array: an array[p-1] of real values in [0,1], 
                 each is your Beta^(1/(prem-2)) random draw for iteration i.
      - theta_array: an array[p-1] of vector[p-2]. 
          In iteration i, we only use the first (prem_i - 2) elements,
          where prem_i = p - (i - 1).
      - sign_flip_array: an array[p-1] of real in {-1, +1}.
        Each iteration i uses sign_flip_array[i].
      - (Optional) fresh_QR logic can be handled externally if needed.
    
    Returns:
      - A p x p matrix C_samp that is analogous to your R function's output.
  */
matrix construct_corrmat(
    vector L,
    vector y_array,
    vector theta_values,
    vector sign_flip_array
) {
    int p = num_elements(L);
    matrix[p, p] V = rep_matrix(0.0, p, p);
    matrix[p, p] B_curr = diag_matrix(rep_vector(1.0, p));
    
    int pos = 1; // Tracks position in flattened theta_values
    for (i in 1:(p-2)) {  // Iterate up to p-1, leaving the last case for special handling
      int pcurr = i - 1;
      int prem = p - pcurr;
      int used  = p - i - 1;  // Number of sphere coordinates this iteration

      // Extract the relevant part of theta_values
      vector[prem - 2] truncated_theta;
      truncated_theta = segment(theta_values, pos, used);
      pos += used;

      // Build V_sub = p x pcurr
      matrix[p, pcurr] V_sub = (pcurr > 0) ? V[, 1:pcurr] : rep_matrix(0.0, p, 0);

      // Call intersection_constraint(...)
      matrix[p, prem] out = intersection_constraint(
        L,
        V_sub,
        B_curr[, 1:prem],
        y_array[i],
        truncated_theta,
        sign_flip_array[i]
      );
    
      V[, i] = out[, 1];
      B_curr[, 1:(prem - 1)] = out[, 2:prem];
    }
    
    //handle 2 x 2 case... hack it with intersection_constraint for now
    matrix[p, 2] out = intersection_constraint(
      L,
      V[, 1:(p-2)],
      B_curr[, 1:2],
      y_array[p-1],
      segment(theta_values, 1, 1),
      sign_flip_array[p-1]
    );

    
    V[, p-1] = out[, 1];
    V[, p] = out[, 2]; // Fill last column with the remaining single column of B_curr
    

    matrix[p, p] eV = V';
    matrix[p, p] diagL = diag_matrix(L);
    matrix[p, p] C_samp = eV * diagL * eV';

    return C_samp;
  }
  
  // -
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
  // -
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
    // The order of QD is ascending
    matrix[prem, prem] QV_raw = eigenvectors_sym(Q);
    vector[prem] sign_first_row;
    for(iv in 1:prem) {
      sign_first_row[iv] = QV_raw[1,iv] < 0  ? -1 : 1;
    }
    matrix[prem, prem] QV = diag_post_multiply(QV_raw, sign_first_row);
    // we fix the first element of QV to be positive to avoid teleportation

    // 4) Sample a point on the intersection for the prem-dim ellipsoid QD
    //    Note: "theta" should be unit_vector[prem-2], "y" in [0,1].
    vector[prem] evec;
    if (prem >= 3) {
      // Standard case: sample a point on the intersection
      evec = sample_intersection(QD, y, theta, sign_flip);
    } else {
      // Solve 2x2 system: A * x_d = b
      matrix[2, 2] A;
      vector[2] b;
      A[1, 1] = 1;
      A[1, 2] = 1;
      A[2, 1] = QD[1];
      A[2, 2] = QD[2];
      b[1] = 1;
      b[2] = 1;
      
      vector[2] x_d = mdivide_left(A, b);

      // Construct final sampled point (mirrored Beta equivalent)
      evec[1] = sqrt(x_d[1]);
      evec[2] = sqrt(x_d[2]) * sign_flip;
    }
    
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
    
    
  
  // --
  //  sample_intersection: already defined in your code
  //  (Accepts vector Lambda, real y in [0,1], unit_vector theta,
  //   real sign_flip in {-1,1}, returns a vector on the ellipse-sphere intersection)
  // --
  vector sample_intersection(vector Lambda, real y, vector theta, real sign_flip) {
    int p = num_elements(Lambda);
    int k = p - 2;  // should match num_elements(theta)

    // 1) Identify indices of min(Lambda) and max(Lambda)
    int i_min = 1;
    int i_max = p;
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

    real sb2_1 = (sb2_sq[1] > 0) ? sqrt(sb2_sq[1]) : 1E5;
    real sb2_2 = (sb2_sq[2] > 0) ? sqrt(sb2_sq[2]) : 1E5;

    // s is the min of sa, sb1, sb2_1, sb2_2
    real s = fmin(sa, fmin(sb1, fmin(sb2_1, sb2_2)));

    // 5) Scale theta by (y * s)
    vector[p-2] x_k = theta * ((1-y) * s); 
    // we want y to scale x_d, so that sign flips are continuous

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
    
    real det_A = determinant(A);
    
    vector[2] x_d = mdivide_left(A, b);

    // 7) Insert these into the final p-vector
    vector[p] x = rep_vector(0.0, p);

    real x_d1 = sqrt(x_d[1]) * sign_flip;
    real x_d2 = sqrt(x_d[2]);

    x[i_min] = x_d1;
    x[i_max] = x_d2;
    for (i in 1:(p-2)) {
      x[ remain_index[i] ] = x_k[i];
    }

    return x;
  }

}

data {
  int<lower=1> n;           // number of observations
  int<lower=3> p;           // dimension
  matrix[n, p] x;           // data matrix, each row is one observation
  int<lower=0> incl_ll;
}

parameters {
  //-
  // 1) Dirichlet for L
  //    Using new array syntax: 'array[p] real' is replaced by 'simplex[p]'
  //-
  simplex[p] raw_dirichlet;

  //-
  // 2) Signed y in [-1,1]
  //    We'll do a truncated normal(0,0.2) for each entry.
  //    Use new array syntax with lower/upper bounds:
  //-
  array[p - 1] real<lower=-1, upper=1> signed_y_array;

  // 3) The "sphere directions" in hyperspherical 
  // (can't make dynamic tuple of 'unit_vector's)
  vector[((p - 2) * (p - 1)) %/% 2] raw_theta_values;

}

transformed parameters {
 //-
  // 1) Transform raw_dirichlet into a decreasing vector L that sums to p
  vector[p] L;

  // Compute cumulative sum
  L[1] = raw_dirichlet[1];
  for (i in 2:p) {
    L[i] = L[i-1] + raw_dirichlet[i];
  }
  
  // Normalize and scale by p
  real sum_L = sum(L);
  for (i in 1:p) {
    L[i] = (L[i] / sum_L) * p;
  }

  // Construct L by reversing the cumulative sum differences

  //-
  // 2) Build y_array in [0,1] plus sign_flip_array from signed_y_array
  //-
  vector[p - 1] y_array;
  vector[p - 1] sign_flip_array;
  for (i in 1:(p - 1)) {
    y_array[i] = abs(signed_y_array[i])^(1.0/i);       // in [0,1]
    sign_flip_array[i] = signed_y_array[i] < 0  // sign in {-1,+1}
                        ? -1
                        : 1;
  }

  //-
  // 3) Flatten the partial "theta angles" into one big vector "theta_values"
  //    We'll have sum_{i=1..(p-1)}(p - i - 1) = (p-1)(p-2)/2 total entries.
  //-
  vector[((p - 1) * (p - 2)) %/% 2] theta_values;
  {
    int pos = 1;
    for (i in 1:(p - 1)) {
      int used = p - i - 1;
      if (used > 0) {
        // take the first 'used' components of theta_array[i]
        theta_values[pos : (pos + used - 1)] = raw_theta_values[pos : (pos + used - 1)] / 
                                               sqrt(sum(raw_theta_values[pos : (pos + used - 1)]^2));
        pos += used;
      }
      // if used == 0, that iteration doesn't need any sphere coords
    }
  }
  
  // 4) Build the correlation matrix via your custom prior
  matrix[p, p] C_samp = construct_corrmat(L, y_array, theta_values, sign_flip_array);
  
}

model {
  
  // 1) raw_dirichlet => uniform simplex
  raw_dirichlet ~ dirichlet(rep_vector(10, p));
  
  // 2) signed_y_array => truncated normal(0,0.2)
  for (i in 1:(p - 1)) {
    signed_y_array[i] ~ normal(0, 1);
  }
  
  // 3) theta_values => uniform on unit sphere
  raw_theta_values ~ std_normal();
  
  // x[i] ~ Multivariate Normal(0, C_samp)
  if(incl_ll == 1){
    for (i in 1:n) {
      real logp = multi_normal_lpdf(x[i] | rep_vector(0, p), C_samp);
      target += logp;
    }  
  }
}

generated quantities {
  // For convenience, let's store the correlation of the first two variables
  real cor_12 = C_samp[1, 2];
}
