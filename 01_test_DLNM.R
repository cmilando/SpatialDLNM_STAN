library(dlnm)
library(splines)
library(rstan)
library(cmdstanr)
library(tidyverse)

cb1.temp <- onebasis(chicagoNMMAPS$temp, fun = 'bs', degree = 1, knots = 21)

# ****
# ok so instead of this, I want to pass this model to STAN
model1 <- glm(death ~ cb1.temp, family=quasipoisson(), chicagoNMMAPS)

summary(model1)
exp(confint(model1))

pred1.temp <- crosspred(cb1.temp, model1, cen=21, by=1)

plot(pred1.temp, "overall")

#
rr <- 1:nrow(chicagoNMMAPS)

# nrows
N = as.integer(nrow(chicagoNMMAPS[rr, ]))

# beta and the intercept
K = as.integer(ncol(cb1.temp) + 1)

# include the intercept
X = matrix(data = cbind(rep(1, times = length(rr)), cb1.temp[rr, ]), 
           nrow = N, ncol = K)
# outcome
y = as.integer(chicagoNMMAPS$death)[rr]

stan_data <- list(
  N = N, 
  K = K, 
  X = X, 
  y = y
)

# Set path to model
stan_model <- cmdstan_model("DLNM.stan")

out1 <- stan_model$sample(
  data = stan_data,
  chains = 1,
  parallel_chains = 1, 
)

draws_array <- out1$draws()

# Convert to data.frame (flattened, easier to use like extract())
draws_df <- posterior::as_draws_df(draws_array)
head(draws_df)

# sick that seems to work
apply(draws_df %>% select(starts_with("beta")), 2, median)
coef(model1)

