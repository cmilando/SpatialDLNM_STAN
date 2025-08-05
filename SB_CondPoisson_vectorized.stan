data {
    // First for the spatial component you'll need something similar
  int<lower=1> J; // regions, and each of the params will need a J
  matrix[J, J] Jmat; // Has to be a matrix so you can do math on it
  
  // Basic params
  int<lower=1> N;               // Number of rows of data
  int<lower=1> K;               // Number of model coefficients, NO intercept
  array[N, K, J] real X;        // matrix of model parameters, NO intercept
  array[N, J] int<lower=0> y;   // outcomes
  
  // Now for the strata,  you need some way to signal which rows
  // are in the strata and subset to those
  // -- might be easy just to do a matrix multiplication but with 1s and 0s
  // NOTE IMPORTANT: THESE MUST BE THE SAME FOR EACH STRATA !!!
  matrix[N, N] S;              // the matrix of strata, has to be matrix so you can do math
  int<lower=1> n_strata;
  int<lower=1> max_in_strata;
  array[n_strata, max_in_strata] int S_condensed;
  array[N] int<lower=0> stratum_id;
}

parameters {
  
  // Ok there is a single centering value for each beta
  //  and a single centering value for sigma
  // vector[K] mu;
  // vector<lower=0>[K] sigma;

  // the Leroux value, just one
  // real<lower=0,upper=1> q;
  
  // and then you need Beta* which has dimension K by J
  // this needs to be defined here so it carries through interations
  matrix[K, J] beta_star; 
  
}

model {
  
  // set priors
  // mu ~ normal(0, 1);
  // q ~ normal(0, 1); // limits are given above
  // sigma ~ normal(0, 1); // limits are given above
  //to_vector(z) ~ normal(0, 1);
  
  //
  array[max_in_strata] int all_indices;
  int k_not_zero;
  
  //
  /*
  vector[K] beta_star_sum;
  real n_a;
  real denom;
  vector[K] star_mean;
  vector[K] star_sd;
  */
    
  // updated coefficients
  vector[K] beta;
    
  // element-wise estimates
  vector[N] theta;
  vector[N] xBeta;
  vector[N] denominator;
  
  // ********************************
  // CALCUATE BETA_STAR
  // ********************************
  // ********************************
  // OK SO THE MAIN DIFFERENCE TO GET TO Spatial is to adjust the Betas
  // by their neighbors
  // so
  // from Ballester: B[a,c] = mu + Beta*[a,c]
  // Beta*[a,c] ~ Normal( q / (1 - q + q*n_a) * Sum(B_a,c w/o B), 
  // VARIANCE =                     sigma_c^2  / (1 - q + q*n_a) )
  // so you will have to take the square-root to get sd
  // ** J is region
  // ** K is beta cofficient
  // NOTE: this has to be done in transformed params because the beta_star_sum
  // updates each time, so it would be hard to do this also in parallel
  // BUT MAYBE NOT? Unclear
  /*
  for(j in 1:J) {
    
    // print("********************************************************");
    // print("j = ", j);
    

    for(k in 1:K) {
    
    // So FIRST, get the sum of BETA STAR neighors that aren't the current one
    // this gets the dot product
    // beta_star needs to be WITHIN k because you are going within each coefficient
    // but across space (so across the j dimension), and then ` transpose
    // OK SO JMAT[j,] Is 1xJ
    // AND BETA_STAR is KxJ, so its JxK
    // MEANING THE PRODUCT IS 1 x K
    // print("beta_star = ", beta_star);
    // print("Jmat = ", Jmat[j, ]);
    // to_matrix enables matrix multiplication
    // to vector makes it easy to deal with later and smaller
      beta_star_sum[k] = Jmat[j, ] * beta_star[k,]';  
      // print("beta_star_sum = ", beta_star_sum);
    
      // next get n_a, which is just the sum of this row of Jmat
      // HMM DOES THIS INCLUDE ITSELF? Or is this what the 1 is for?
      // You could always add 1 if so
      n_a = sum(Jmat[j, ]);
      // print("n_a = ", n_a);
        
      // now you should be able to get the mean and var
      denom = 1 - q + q * n_a;
      // print("q = ", q);
    
      
      // Now contruct the beta star mean and sigma
      // remember to square root denom in sigma
      // since in the paper its for the variance
      // element-wise multiplication
      star_mean[k] = (q / denom) * beta_star_sum[k];
      // print("star_mean = ", star_mean);
      
      // element-wise division
      // and set a minimum
      // unsure if the square root is necesssary, they don't do it
      // so it must be a typo in the manuscript
      star_sd[k] = sigma[k] / sqrt(denom);
      // print("star_sd = ", star_sd);
  
      // vectorized draws from normal distributions
      beta_star[k, j] ~ normal(star_mean[k], star_sd[k]);
      // beta_star[, j] = star_mean + z[, j] .* star_sd;
      // print("beta_star = ", beta_star);
    }
  }
  */  
  //////////////////////////////////////////////////////
  // SEND TO REDUCE SUM ////////////////////////////////
  //////////////////////////////////////////////////////

  
  for(j in 1:J) {
 
    // ********************************
    // CALCUATE THETA
    // ********************************
    // if there are some variables that are not spatial you could do this differently
    // for some, so ifelse is_spatial 1/0 eiterh mu or beta_star etc
    // you don't need to do non-centered here since you 
    // are drawing directly from beta_star
    // Note: why doesn't this include the sigma[j] - I think this
    // was not implemented correctly in their winbugs model
    
    // (1) full model
    // beta = mu + beta_star[,j];
    // print("mu = ", mu);
    // print("beta = ", beta);
    
    // (2) independent
    beta = beta_star[,j];
    
    // (3) one beta to rule them all
    //beta = mu;
  
    // ********************************
    // from Armstrong 2014, equation (4)
    // theta = exp(X*beta) / sum( exp(X*beta) for all strata)
    
    // first get get the numerator:
    // ok so X[,,J] is N x K, and beta is K x 1
    // so this turns into N x 1
    // update each sum to control for overflow
    xBeta = exp(fmin(20, fmax(to_matrix(X[, , j]) * beta, -20)));
    // print("xBeta = ", xBeta);
    
    // then I think with matrix math you can get the bottom in one shot
    // S is N x N and xBeta is N x 1
    // so matrix multiplication gives you N x 1
    denominator = S * xBeta;
    // print("denominator = ", denominator);
    
    // now get theta
    // have to use element division
    // gives the THETA for THIS REGION ONLY
    theta = xBeta ./ denominator;
    // print("theta = ", theta);
    
  
    // and get the conditional model
    for (i in 1:n_strata) {
       
       // first get all the strata, which include some spurious 0s
       all_indices = to_array_1d(S_condensed[i, ]);
       
       // now, subset to just the ones that are not 0
      k_not_zero = 0;
       for(k in 1:max_in_strata) {
         if(S_condensed[i, k] > 0) k_not_zero += 1;
       }

       // check that sum theta is ALWAYS = 1 and
       // REMEMBER TO EXCLUDE ANY EMPTY STRATA TO AVOID BIAS
       if(sum(y[all_indices[1:k_not_zero], j]) > 0) {
         
         // // printout
         if(is_nan(sum(theta[all_indices[1:k_not_zero]]))) {
           // print("***********************************");
           // print("beta_star = ", beta_star);
         } else {
           
           // just get the values for this strata
           target += multinomial_lpmf(y[all_indices[1:k_not_zero], j] | 
                                        theta[all_indices[1:k_not_zero]]);
         }
       
       }
    } 
  } // J
}

generated quantities {
  
  // (1) make a new BETA by randomly sampling from mu and beta_star
  matrix[K, J] beta_out;
  
  for(k in 1:K) {
    for(j in 1:J) {
      // beta_out[k,j] = mu[k] + beta_star[k,j]; // probably some additional variance here
      beta_out[k,j] = beta_star[k,j]; // probably some additional variance here
    }
  }
  
  /*
  // (2) apparently you can handle over-dispersion in post-processing as per STATA
  // code, so lets just do that instead of direhclt, which could be another
  // options
  /// Lets see if this works
  /// Converted to STAN from STATA code from Gasp and Armstrong using ChatGPT
  /// and Qc'd
  array[J] real dispersion;
  array[J] real pearson_x2;
  int df_resid = N - K - n_strata;

  array[N] real pred_rescaled;     // rescaled to match stratum totals

  for (j in 1:J) {
    
    // predicted counts before stratum normalization
    vector[N] xBeta_out = exp(to_matrix(X[, , j]) * beta_out[,j]);  

    vector[n_strata] sum_y_stratum = rep_vector(0, n_strata);
    vector[n_strata] sum_pred_stratum = rep_vector(0, n_strata);

    // Sum observed and predicted counts per stratum
    for (n in 1:N) {
      int s = stratum_id[n];
      sum_y_stratum[s] += y[n, j];
      sum_pred_stratum[s] += xBeta_out[n];
    }

    // Rescale predictions so that the sum of the predicted in each group
    // equals the sum of the observed in each group
    for (n in 1:N) {
      int s = stratum_id[n];
      if (sum_pred_stratum[s] > 1e-8) {
        pred_rescaled[n] = xBeta_out[n] * sum_y_stratum[s] / sum_pred_stratum[s];
      } else {
        pred_rescaled[n] = xBeta_out[n];
      }
    }

    // Compute Pearson χ² for this region
    pearson_x2[j] = 0;
    for (n in 1:N) {
      if (pred_rescaled[n] > 1e-8) {
        pearson_x2[j] += square(y[n, j] - pred_rescaled[n]) / pred_rescaled[n];
      }
    }

    // Compute dispersion for this region
    dispersion[j] = pearson_x2[j] / df_resid;
  }
  */
}
