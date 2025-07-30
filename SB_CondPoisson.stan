data {
  
  // Basic params
  int<lower=1> N;          // Number of rows of data
  int<lower=1> K;          // Number of model coefficients, NO intercept
  matrix[N, K] X;          // matrix of model parameters, NO intercept
  array[N] int<lower=0> y;         // fractional outcomes, out i / sum(out i in strata)
  
  // Now for the strata,  you need some way to signal which rows
  // are in the strata and subset to those
  // -- might be easy just to do a matrix multiplication but with 1s and 0s
  matrix[N, N] S;              // the matrix of strata, 
  int<lower=1> n_strata;
  int<lower=1> max_in_strata;
  array[n_strata, max_in_strata] int S_condensed;
  
  // Finally for the spatial component you'll need something similar
  int<lower=1> J; // regions, and each of the params will need a J
}

parameters {
  
  // Ok so instead, there is a single centering value for each beta
  //  and a single centering value for sigma
  vector[K] mu;
  array[K] real<lower=0> sigma;
  array[K] real<lower=0> star_mean;
  
  // the Leroux value, just one
  real<lower=0,upper=1> q;
  
  // and then you need Beta* which has dimension K by J
  matrix[K, J] beta_star; 
  
}

model {
  
  // set priors
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);
  star_mean ~ normal(0, 1);
  

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
      // ** so each has all K, but only some other J
      // so there should be a J index lookup function that gets each one's neighbors
      // so for J = 1, that should return the indices of its neigbhors so 2, 3, 4
      // and 
      //real Beta_Sum = 1; // this will eventually be a result of the lookup
      //real n_a = 1;
      //real star_mean = 0;
      //real star_sd = sigma[k] / (1 - q + q * n_a);
      beta_star[k,j] ~ normal(star_mean[k], sigma[k]);
    }
  }
  
  for(j in 1:J) {
    // finally
    // if there are some variables that are not spatial you could do this differently
    // for some, so ifelse is_spatial 1/0 eiterh mu or beta_star etc
     vector[K] beta = mu + beta_star[,j];
  
    // ********************************
    // from Armstrong 2014, equation (4)
    // theta = exp(X*beta) / sum( exp(X*beta) for all strata)
    
    // first get get the numerator:
    // ok so X is N x K, and beta is K x 1
    // so this turns into N x 1
    vector[N] xBeta = exp(X * beta);
    
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
      
       // REMEMBER TO EXCLUDE ANY EMPTY STRATA TO AVOID BIAS
       if(sum(y[my_array]) > 0) {
      
         // just get the values for this strata
         target += multinomial_lpmf(y[my_array] | theta[my_array]);
       
       }
    } 
  } // J
}

// apparently you can handle dispersion in post-processing as per STATA
// code, so lets just do that instead of direhclt, which could be another
// options
