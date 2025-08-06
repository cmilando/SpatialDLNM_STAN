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
  matrix[K, J] beta;  // attribute effects 
}

model {

  // now set priors
  // to_vector doesn't work ... so don't do it
  for(k in 1:K) {
      beta[k,] ~ normal(0, 5);
  }

  matrix[N, J] xBeta;
  matrix[N, J] denominator; 
  matrix[N, J] theta;
  
  for(j in 1:J) {
    // from Armstrong 2014, equation (4)
    // theta = exp(X*beta) / sum( exp(X*beta) for all strata)
    
    // first get get the numerator:
    // ok so X is N x K, and beta is K x 1
    // so this turns into N x 1
    // UPDATE added block to keep exp(inf) or exp(-info)
    // UPDATE to the UPDATE: that block messes things up ... so don't do it
    xBeta[, j] = exp(to_matrix(X[,,j]) * beta[,j]);
  
    // then I think with matrix math you can get the bottom in one shot
    // S is N x N and xBeta is 
    denominator[, j] = S * xBeta[, j];
    
    // now get theta
    // have to use element division
    theta[, j] = xBeta[,j] ./ denominator[, j];
  }

  // and get the model
  for(j in 1:J) {
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
       if(sum(y[my_array, j]) > 0) {
          
        //  
          if(is_nan(sum(theta[my_array, j]))) {
            reject("NAN THETA SUM");
          } else {
          
           // just get the values for this strata
           y[my_array, j] ~ multinomial(theta[my_array, j]);
          }
       
       }
    }
  }
}


