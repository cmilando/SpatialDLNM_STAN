data {
  
  // Basic params
  int<lower=1> N;          // Number of rows of data
  int<lower=1> K;          // Number of model coefficients, NO intercept
  matrix[N, K] X;          // matrix of model parameters, NO intercept
  array[N] int<lower=0> y;         // fractional outcomes, out i / sum(out i in strata)
  
  // Now for the strata,  you need some way to signal which rows
  // are in the strata and subset to those
  // -- might be easy just to do a matrix multiplication but with 
  //    1s and 0s
  matrix[N, N] S;              // the matrix of strata, 
  int<lower=1> n_strata;
  int<lower=1> max_in_strata;
  array[n_strata, max_in_strata] int S_condensed;
}

parameters {
  vector[K] beta;  // attribute effects 
}

model {
  
  // set priors
  beta ~ normal(0, 5);
  
  // ********************************
  // OK SO THE MAIN DIFFERENCE TO GET TO Spatial is to adjust the Betas
  // by their neighbors
  // 
  
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
       y[my_array] ~ multinomial(theta[my_array]);
     
     }
  } 
}

// apparently you can handle dispersion in post-processing as per STATA
// code, so lets just do that instead of direhclt, which could be another
// options
