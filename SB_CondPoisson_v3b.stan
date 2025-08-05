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
  
  int<lower=1> model_type;
}

transformed data {
  
  vector[J] n_a;
  
  for(j in 1:J) {
    //// ********************************************************
    // PRECOMPUTE THIS BLOCK IN TRANSFORMED DATA
    // you could even pre-compute this in transformed data
    // next get n_a, which is just the sum of this row of Jmat
    // HMM DOES THIS INCLUDE ITSELF? Or is this what the 1 is for?
    // You could always add 1 if so
    n_a[j] = sum(Jmat[j, ]);
    //// ********************************************************
  }
  
}

parameters {
  // Ok so instead, there is a single centering value for each beta
  //  and a single centering value for sigma
  vector[K] mu;
  vector<lower=0>[K] sigma; // lower limit is 1e-6 so it can NEver be 0

  // the Leroux value, just one
  real<lower=0,upper=1> q;
  
  // and then you need Beta* which has dimension K by J
  matrix[K, J] beta_star;  // attribute effects 
}

model {
  
  // --------------------------------------------------------------------------
  // PRIORS
  // --------------------------------------------------------------------------
  mu ~ normal(0, 5);
  q ~ cauchy(0, 5);
  sigma ~ gamma(2, 2);

  // --------------------------------------------------------------------------
  // GET BETA STAR
  // Note: I tried and tried to get this to work with a z variable for standard
  //       but it just did not work, I think because of the sum
  //       like the star_mean depends on itself.
  // --------------------------------------------------------------------------
  // real beta_star_sum;
  vector[K] beta_star_sum;
  real beta_star_denom;
  vector[K] star_mean;
  vector[K] star_sd;
  
  // ** J is region
  // ** K is beta cofficient
  
  for(j in 1:J) {
    
    // a single value for the denominator applied to 
    // beta star mean and st_dev
    beta_star_denom = 1 - q + q * n_a[j];
    
    // So FIRST, get the sum of BETA STAR neighors that aren't the current one
    // this gets the dot product
    // beta_star needs to be WITHIN k because you are going within each coefficient
    // but across space (so across the j dimension), and then ` transpose
    // hmm - this is updated each time, which is probably not correct
    // Jmat[, j] is  x J
    // beta_star is K x J so transpose is J x K
    // to product is 1 x K;
    // it likes it in vectors so tranpose again? this might be too expensive
    // this works because beta_star has initial values because its a parameter
    beta_star_sum = beta_star * Jmat[, j]; 
  
    // Now contruct the beta star mean and sigma
    // remember to square root denom in sigma
    // since in the paper its for the variance
    // make sure to use element division and multiplication
    star_mean = q / beta_star_denom .* beta_star_sum;
    star_sd = sigma ./ sqrt(beta_star_denom);
  
    // and re-draw across K
    beta_star[,j] ~ normal(star_mean, star_sd);
    
  }

  // --------------------------------------------------------------------------
  // and get the target
  // THIS CAN ALL GO TO REDUCE SUM I THINK !
  // --------------------------------------------------------------------------
  vector[K] beta;
  vector[N] xBeta;
  vector[N] theta_denominator; 
  vector[N] theta;
  
  for(j in 1:J) {
    
    // **************************************************************
    // get beta
    
    // (1) fully spatial model
    //if(model_type == 1) {
      beta = mu + beta_star[,j];
    //}

    // (2) Indepdent
    //if(model_type == 2) {
    //  beta = beta_star[,j];
    //}
    
    // (3) shared
    //if(model_type == 3) {
    //  beta = mu;
    //}
 
    // **************************************************************
    // from Armstrong 2014, equation (4)
    // theta = exp(X*beta) / sum( exp(X*beta) for all strata)
    
    // first get get the numerator:
    // ok so X is N x K, and beta is K x 1
    // so this turns into N x 1
    // UPDATE added block to keep exp(inf) or exp(-info)
    // UPDATE to the UPDATE: that block messes things up ... so don't do it
    xBeta = exp(to_matrix(X[,,j]) * beta);
  
    // then I think with matrix math you can get the bottom in one shot
    // S is N x N and xBeta is 
    theta_denominator = S * xBeta;
    
    // now get theta
    // have to use element division
    theta = xBeta ./ theta_denominator;
    
    for (i in 1:n_strata) {
       
       //// ********************************************************
       // PRECOMPUTE THIS BLOCK IN TRANSFORMED DATA
       // will have to do something like this: 
       // https://mc-stan.org/docs/2_21/stan-users-guide/ragged-data-structs-section.html
       // first get all the strata, which include some spurious 0s
       array[max_in_strata] int all_indices = to_array_1d(S_condensed[i, ]);
       
       // now, subset to just the ones that are not 0
       int k_not_zero = 0;
       for(k in 1:max_in_strata) {
         if(S_condensed[i, k] > 0) k_not_zero += 1;
       }
       array[k_not_zero] int my_array = all_indices[1:k_not_zero];
       //// ********************************************************
      
       // REMEMBER TO EXCLUDE ANY EMPTY STRATA TO AVOID BIAS
       if(sum(y[my_array, j]) > 0) {
        
         if(is_nan(sum(theta[my_array]))) {
           reject("CHAD TEST: rejecting because sum-theta is nan");
         } else {
           // just get the values for this strata
           target += multinomial_lpmf(y[my_array, j] | theta[my_array]);
         }
       
       }
    }
  }
}

