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
  vector[K] beta;  // attribute effects 
}

model {

  // now set priors
  beta ~ normal(0, 5);

  // from Armstrong 2014, equation (4)
  // theta = exp(X*beta) / sum( exp(X*beta) for all strata)
  
  // first get get the numerator:
  // ok so X is N x K, and beta is K x 1
  // so this turns into N x 1
  // UPDATE added block to keep exp(inf) or exp(-info)
  // UPDATE ** this really slows things down
  vector[N] xBeta = exp(X * beta);

  // then I think with matrix math you can get the bottom in one shot
  // S is N x N and xBeta is 
  vector[N] denominator = S * xBeta;
  
  // now get theta
  // have to use element division
  vector[N] theta = xBeta ./ denominator;

  // and get the model
  for (i in 1:n_strata) {
     
      array[strata_len[i]] int my_array = S_condensed[i, 1:strata_len[i]];
        
    
     // REMEMBER TO EXCLUDE ANY EMPTY STRATA TO AVOID BIAS
     if(sum(y[my_array]) > 0) {
    
       // just get the values for this strata
       y[my_array] ~ multinomial(theta[my_array]);
     
       }
  
  }
}


