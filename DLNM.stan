data {
  int<lower=1> N;          // Number of rows of data
  int<lower=1> K;          // Number of model coefficients, includes intercept
  matrix[N, K] X;          // matrix of model parameters, includes intercept
  array[N] int<lower=0> y; // outcomes
}

parameters {
  vector<lower=-20,upper=20>[K] beta;  // beta values, bounded to avoid expInf or -Inf
  real<lower=1e-6,upper=10> phi;       // Dispersion parameter (larger = closer to Poisson)
}

transformed parameters {
  // avoids exp(-inf) or exp(inf)
  vector[N] mu_raw = X * beta;
  vector[N] mu;
  for (n in 1:N)
    mu[n] = exp(fmin(20, fmax(mu_raw[n], -20)));
}

model {

  // Priors
  beta ~ normal(0, 1);
  phi ~ gamma(2, 0.1);
  
  // Likelihood: NegBinomial_2 parameterized with mean and overdispersion
  y ~ neg_binomial_2(mu, phi);
}
