data {
  int<lower=1> N;          // Number of rows of data
  int<lower=1> K;          // Number of model coefficients, includes intercept
  matrix[N, K] X;          // matrix of model parameters, includes intercept
  array[N] int<lower=0> y; // outcomes
}

parameters {
  vector[K] beta;    // beta values
  real<lower=0> phi; // Dispersion parameter (larger = closer to Poisson)
}

model {
  //
  vector[N] mu = exp(X * beta);
  
  // Priors
  beta ~ normal(0, 5);
  phi ~ exponential(1);

  // Likelihood: NegBinomial_2 parameterized with mean and overdispersion
  y ~ neg_binomial_2(mu, phi);
}
