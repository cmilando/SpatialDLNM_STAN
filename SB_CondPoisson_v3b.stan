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
  vector<lower=0>[K] sigma;

  // the Leroux value, just one
  real<lower=0,upper=1> q;
  
  // and then you need Beta* which has dimension K by J
  matrix[K, J] beta_star;  // attribute effects 
}

transformed parameters {
  matrix[K, J] beta;
  matrix[N, J] xBeta;
  matrix[N, J] xBeta_inner;
  matrix[N, J] denominator; 
  matrix[N, J] theta;
  
  for(j in 1:J) {
    
    // get beta
    beta[,j] = mu + beta_star[,j];
    
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
  
}

model {
  
  // --------------------------------------------------------------------------
  // PRIORS
  // --------------------------------------------------------------------------
  mu ~ normal(0, 5);
  q ~ normal(0.5, 1);
  sigma ~ normal(0, 1);
  
  // to_vector doesn't work ... so don't do it
  // for(k in 1:K) {
  //    z[k,] ~ normal(0, 5);
  //}
  
  // --------------------------------------------------------------------------
  // GET BETA STAR
  // --------------------------------------------------------------------------
  real beta_star_sum;
  real n_a;
  real denom;
  real star_mean;
  real star_sd;
  
  for(j in 1:J) {
    for(k in 1:K) {
      // ** J is region
      // ** K is beta cofficient
      
      // So FIRST, get the sum of BETA STAR neighors that aren't the current one
      // this gets the dot product
      // beta_star needs to be WITHIN k because you are going within each coefficient
      // but across space (so across the j dimension), and then ` transpose
      // hmm - this is updated each time, which is probably not correct
      beta_star_sum = Jmat[j, ] * beta_star[k, ]'; 
      
      // next get n_a, which is just the sum of this row of Jmat
      // HMM DOES THIS INCLUDE ITSELF? Or is this what the 1 is for?
      // You could always add 1 if so
      n_a = sum(Jmat[j, ]);
      
      // now you should be able to get the mean and var
      denom = 1 - q + q * n_a;
      
      // Now contruct the beta star mean and sigma
      // remember to square root denom in sigma
      // since in the paper its for the variance
      star_mean = q / denom * beta_star_sum;
      star_sd = sigma[k] / sqrt(denom);
      
      // and re-draw
      beta_star[k,j] ~ normal(star_mean, star_sd);
    }
  }

  // --------------------------------------------------------------------------
  // and get the target
  // --------------------------------------------------------------------------
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
      
         // just get the values for this strata
         y[my_array, j] ~ multinomial(theta[my_array, j]);
       
       }
    }
  }
}

generated quantities {
  
  // (1) make a new BETA by randomly sampling from mu and beta_star
  matrix[K, J] beta_out;
  
  for(k in 1:K) {
    for(j in 1:J) {
      beta_out[k,j] = mu[k] + beta_star[k,j]; // probably some additional variance here
    }
  }
}
