functions {
  /*
    sample_intersection:
      Replicates the core logic of your R function but expects:
        - Lambda: coefficients of the hyperellipsoid (length p).
        - y: a real in [0, 1], which you've already sampled externally 
             as (Beta draw)^(1/k).
        - theta: a unit_vector[k], with k = p-2 (the sphere direction).
        - sign_flip: real in {-1, 1} for flipping the sign 
                     of the second dropped coordinate.
      Returns:
        - x: a unit_vector[p] that lies on the intersection of the unit sphere 
             and the hyperellipsoid with coefficients Lambda.
  */
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
}

data {
  int<lower=3> p;          
  vector[p] Lambda;        
  real<lower=0, upper=1> y;
  real sign_flip;           // in {-1, 1}
}

parameters {
  // A random direction on the (p-2)-sphere:
  unit_vector[p-2] theta;
}

model {
  // No likelihood or priors for this test
}

generated quantities {
  vector[p] x = sample_intersection(Lambda, y, theta, sign_flip);
}
