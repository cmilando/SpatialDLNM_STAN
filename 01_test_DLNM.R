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


# ****

# int<lower=1> N;
# int<lower=1> K;
# matrix[N, K] X;
# int<lower=0> y[N];

# have to do this to remove NAs
# rr <- (maxlag + 1):nrow(chicagoNMMAPS)

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

init_f <- function() {
  list(beta = rep(5, times = K),
       phi = 1)
}

out1 <- stan_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4, 
  init = init_f
)

 draws_array <- out1$draws()

# Convert to data.frame (flattened, easier to use like extract())
draws_df <- posterior::as_draws_df(draws_array)
head(draws_df)

# sick that seems to work
apply(draws_df %>% select(starts_with("beta")), 2, median)
coef(model1)

# HAHA and then you have to do crosspred by hand lol


## THOUGHTS
# -- TO EXTEND TO MULTIPLE REGIONS
# -- PASS In a 3-D array for beta [N, K, J] where J is regions

# *** the reason INLA takes so long is that there are SO MANY intercepts

# ****
# ***** hmm but you'd have to manually code the gnm
# *** THATS WHERE THE SPREADSHEET COMES IN BABY

# ****
# ***** so you'd need to send in the strata variable somehow, maybe as K
# ***** might be easier to start with time-series




# -- UM JUST LOOK HERE AT PAGE 20
# https://cran.r-project.org/web/packages/CARBayes/vignettes/CARBayes.pdf


# hmmm bc of dlnm you may want to just do this yourself
# and you trust STAN more
# and this Laurexz thing seems easy
# https://academic.oup.com/ije/article/53/3/dyae061/7654027


#
# inla_formula <- visits_n ~ -1 + 
#   cb1 + cb2 + cb3 + cb4 + cb5 + cb6 + 
#   cb7 + cb8 + cb9 +
#   f(strata, model = "iid", hyper = list(prec = list(initial = log(1e-04), fixed = TRUE))) + 
#   f(id_cb1, cb1, model = "bym", graph = list_neig, 
#     hyper = list(prec.unstruct = list(prior = sdunif), 
#                  prec.spatial = list(prior = sdunif)))


# REFERNECES:
# * https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/fit_scotland_bym2.R
# * https://mc-stan.org/learn-stan/case-studies/icar_stan.html




