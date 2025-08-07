functions {
  real partial_sum_lpmf(array[] int dummy, 
                        int start, int end,
                        data array[] int strata_len,
                        data array[,] int S_condensed,
                        data array[,] int y,
                        data int J,
                        matrix theta) {
    
  // THIS CAN ALL GO TO REDUCE SUM I THINK !
  // https://mc-stan.org/docs/stan-users-guide/parallelization.html
  // first try and reduce_sum here
    
    // inputs you need:
    // -- dummy idx as the first argument
    // -- data array strata_len
    // -- data array S_condensed
    // -- 1:n_strata are given by start:end
    // -- data array y
    // -- array theta

                          
    real lupmf = 0;

    for (i in start:end) {
          
      array[strata_len[i]] int my_array = S_condensed[i, 1:strata_len[i]];
        
      for(j in 1:J) {
         // REMEMBER TO EXCLUDE ANY EMPTY STRATA TO AVOID BIAS
         if(sum(y[my_array, j]) > 0) {
          
           if(is_nan(sum(theta[my_array, j]))) {
             //reject("CHAD TEST: rejecting because sum-theta is nan");
             // this happens because exp(xBeta) is -inf or inf but
             // it is too computationally expensive to check each value
             // or to set limits, so we'll just catch it here
           } else {
             // just get the values for this strata
             // USING LUMPF!!! AS PER THIS
             // https://mc-stan.org/learn-stan/case-studies/reduce_sum_tutorial.html
             lupmf += multinomial_lupmf(y[my_array, j] | theta[my_array, j]);
           }
         
         }
      } // j
    } // n_strata

    return lupmf;
  }
}

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
  
  //
  int<lower=1> grainsize; 
}

transformed data {
  
    array[n_strata] int strata_len;
  
  for(n in 1:n_strata) {

       int k_not_zero = 0;
       for(k in 1:max_in_strata) {
         if(S_condensed[n, k] > 0) k_not_zero += 1;
       }
      strata_len[n] = k_not_zero;
  }
  
    array[n_strata] int dummy;
  for(n in 1:n_strata) {
    dummy[n] = n;
  }
}

parameters {
  matrix[K, J] beta;  // attribute effects 
}

model {

  // now set priors
  // to_vector doesn't work ... so don't do it
  to_vector(beta) ~ normal(0, 5);

  matrix[N, J] xBeta;

  for(j in 1:J) {
    // from Armstrong 2014, equation (4)
    // theta = exp(X*beta) / sum( exp(X*beta) for all strata)
    
    // first get get the numerator:
    // ok so X is N x K, and beta is K x 1
    // so this turns into N x 1
    // UPDATE added block to keep exp(inf) or exp(-info)
    // UPDATE to the UPDATE: that block messes things up ... so don't do it
    xBeta[, j] = exp(to_matrix(X[,,j]) * beta[,j]);
  }
  
  
  // then I think with matrix math you can get the bottom in one shot
  // S is N x N and xBeta is 
  matrix[N, J] denominator = S * xBeta;
    
  // now get theta
  // have to use element division
  matrix[N, J] theta = xBeta ./ denominator;
  
  // and get the model
  target += reduce_sum(partial_sum_lupmf, dummy, grainsize,
                    strata_len, S_condensed, y, J, theta);
}


