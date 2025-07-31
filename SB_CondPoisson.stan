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
  
  // Ok so instead, there is a single centering value for each beta
  //  and a single centering value for sigma
  vector[K] mu;
  array[K] real<lower=0> sigma;

  // the Leroux value, just one
  real<lower=0,upper=1> q;
  
  // and then you need Beta* which has dimension K by J
  matrix[K, J] beta_star; 
  
}

model {
  
  // set priors
  mu ~ normal(0, 1);
  q ~ normal(0, 1) T[0, 1]; // add limits here
  sigma ~ normal(0, 1) T[0, ]; // add here limits

  // ********************************
  // OK SO THE MAIN DIFFERENCE TO GET TO Spatial is to adjust the Betas
  // by their neighbors
  // so
  // from Ballester: B[a,c] = mu + Beta*[a,c]
  // Beta*[a,c] ~ Normal( q / (1 - q + q*n_a) * Sum(B_a,c w/o B), 
  // VARIANCE =                     sigma_c^2  / (1 - q + q*n_a) )
  // so you will have to take the square-root to get sd
  for(j in 1:J) {
    for(k in 1:K) {
      // ** J is region
      // ** K is beta cofficient
      
      // So FIRST, get the sum of BETA STAR neighors that aren't the current one
      // this gets the dot product
      // beta_star needs to be WITHIN k because you are going within each coefficient
      // but across space (so across the j dimension), and then ` transpose
      real beta_star_sum = Jmat[j, ] * beta_star[k, ]'; 
      
      // next get n_a, which is just the sum of this row of Jmat
      // HMM DOES THIS INCLUDE ITSELF? Or is this what the 1 is for?
      // You could always add 1 if so
      real n_a = sum(Jmat[j, ]);

      // now you should be able to get the mean and var
      // hm how is this sometimes negative?
      real denom = 1 - q + q * n_a;
      
      // Now contruct the beta star mean and sigma
      // remember to square root denom in sigma
      real star_mean = q / denom * beta_star_sum;
      real star_sd = sigma[k] / (sqrt(denom));
      
      // and re-draw
      beta_star[k,j] ~ normal(star_mean, star_sd);
      
    }
  }

  for(j in 1:J) {
    // finally
    // if there are some variables that are not spatial you could do this differently
    // for some, so ifelse is_spatial 1/0 eiterh mu or beta_star etc
    vector[K] beta = mu + beta_star[,j];
    
    // UPDATED TO INCLUDE SIGMA[J] AFTER LOOKING AT THEIR CODE
    // IS THIS CORRECT ?? DIFFERENT FROM WHAT IT SAYS IN THE PAPER
  
    // ********************************
    // from Armstrong 2014, equation (4)
    // theta = exp(X*beta) / sum( exp(X*beta) for all strata)
    
    // first get get the numerator:
    // ok so X[,,J] is N x K, and beta is K x 1
    // so this turns into N x 1
    vector[N] xBeta_raw = to_matrix(X[, ,j]) * beta;
    vector[N] xBeta;
    for(n in 1:N)
      xBeta[n] = exp(fmin(20, fmax(xBeta_raw[n], -20)));
    
    
    // then I think with matrix math you can get the bottom in one shot
    // S is N x N and xBeta is N x 1
    vector[N] denominator = S * xBeta;
    
    // now get theta
    // have to use element division
    vector[N] theta = xBeta ./ denominator;
    
    // and get the conditional model
    for (i in 1:n_strata) {
       
       // first get all the strata, which include some spurious 0s
       array[max_in_strata] int all_indices = to_array_1d(S_condensed[i, ]);
       
       // now, subset to just the ones that are not 0
       int k_not_zero = 0;
       for(k in 1:max_in_strata) {
         if(S_condensed[i, k] > 0) k_not_zero += 1;
       }
       array[k_not_zero] int my_array = all_indices[1:k_not_zero];
      
       // check that sum theta is ALWAYS = 1 and
       // REMEMBER TO EXCLUDE ANY EMPTY STRATA TO AVOID BIAS
       if(sum(y[my_array, j]) > 0) {
      
         // just get the values for this strata
         target += multinomial_lpmf(y[my_array, j] | theta[my_array]);
       
       }
    } 
  } // J
}

generated quantities {
  
  // (1) make a new BETA by randomly sampling from mu and beta_star
  matrix[K, J] beta_out;
  
  for(k in 1:K) {
    for(j in 1:J) {
      beta_out[k,j] = mu[k] + beta_star[k,j]; // probably some additional variance here
    }
  }
  
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
    vector[N] xBeta = exp(to_matrix(X[, , j]) * beta_out[,j]);  

    vector[n_strata] sum_y_stratum = rep_vector(0, n_strata);
    vector[n_strata] sum_pred_stratum = rep_vector(0, n_strata);

    // Sum observed and predicted counts per stratum
    for (n in 1:N) {
      int s = stratum_id[n];
      sum_y_stratum[s] += y[n, j];
      sum_pred_stratum[s] += xBeta[n];
    }

    // Rescale predictions so that the sum of the predicted in each group
    // equals the sum of the observed in each group
    for (n in 1:N) {
      int s = stratum_id[n];
      if (sum_pred_stratum[s] > 1e-8) {
        pred_rescaled[n] = xBeta[n] * sum_y_stratum[s] / sum_pred_stratum[s];
      } else {
        pred_rescaled[n] = xBeta[n];
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
}

