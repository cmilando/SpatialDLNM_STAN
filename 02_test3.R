library(dlnm)
library(splines)
library(rstan)
library(cmdstanr)
library(tidyverse)
library(foreign) # ENABLES READING THE DATA FILE, WHICH IS A STATA FORMAT

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' Run model
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////
#' 
data <- read.dta("https://raw.github.com/gasparrini/2014_armstrong_BMCmrm_codedata/master/londondataset2002_2006.dta")

# SET THE DEFAULT ACTION FOR MISSING DATA TO na.exclude
# (MISSING EXCLUDED IN ESTIMATION BUT RE-INSERTED IN PREDICTION/RESIDUALS)
options(na.action = "na.exclude")

# SCALE EXPOSURE
data$ozone10 <- data$ozone / 10

# GENERATE MONTH AND YEAR
# changed the factor part here so it wouldn't carry forwards
# subset to reduce the no. of factors
data$month <- months(data$date)
data$year <- format(data$date, format = "%Y")
data$dow <- weekdays(data$date)

# subset for Excel application
data <- subset(data, dow %in% c("Monday", "Tuesday") & 
                 year %in% c("2002", "2003", "2004"))
data$stratum <- interaction(data$year, data$month, data$dow)
length(levels(data$stratum))
data <- data[order(data$date), ]

# FIT UNCONDITIONAL POISSON MODEL
model_upr <- glm(numdeaths ~ ozone10 + temperature + factor(stratum), 
                 data = data, family = poisson)
summary(model_upr)

# FIT A CONDITIONAL POISSON MODEL WITH A YEAR X MONTH X DOW STRATA
library(gnm)
modelcpr1 <- gnm(numdeaths ~ ozone10 + temperature, 
                 data = data, family = poisson, eliminate = factor(stratum))
summary(modelcpr1)
attr(modelcpr1$coefficients, "eliminated")


# yeah this is the one
# https://mc-stan.org/docs/2_20/functions-reference/multinomial-distribution.html
# 

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' SETUP inputs
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////

# nrows
N = as.integer(nrow(data))

# beta and the intercept
K = as.integer(2)

# include the intercept
X1 = as.matrix(data[, c('ozone10', 'temperature')])
head(X)
X2 = X1
X = array(dim = c(dim(X1), 2))
X[,,1] = X1
X[,,2] = X2

# outcome
y1 = as.integer(data$numdeaths)
y2 = y1
y = cbind(y1, y2)
# y <- matrix(y[, 1], nrow = length(y1), ncol = 1)

# create S matrix 
getSmat <- function(strata_vector, include_self = T) {
  
  # strata_vector <- data$stratum 
  strata_matrix <- matrix(as.integer(strata_vector), 
                          nrow = length(strata_vector),
                          ncol = length(strata_vector), 
                          byrow = T)
  
  for(i in 1:length(strata_vector)) {
    strata_matrix[i, ] = 1*(strata_matrix[i, ] == as.integer(strata_vector[i]))
    if(!include_self) {
      strata_matrix[i, i] = 0
    }
  }
  
  return(strata_matrix)
}

#
S <- getSmat(data$stratum, include_self = T)
head(S)


# get strata vars
n_strata <- as.integer(length(unique(data$stratum))) # 72, cool
S_list <- apply(S, 1, function(x) which(x == 1))
max_in_strata <- max(sapply(S_list, length))
S_list <- lapply(S_list, function(l) {
  if(length(l) == max_in_strata) {
    return(l)
  } else {
    diff_n = max_in_strata - length(l)
    return(c(l, rep(0, times = diff_n)))
  }
})
S_condensed <- unique(do.call(rbind, S_list))
dim(S_condensed)

# ok now do the same as include self but with J matrix
# start simple and you can generalize later
# this is an adjacency matrix that DOES NOT include itself
Jvec <- c(1, 2)
Jmat <- matrix(as.integer(c(0, 1, 1, 0)), nrow = 2)
# Jmat <- matrix(Jmat[1, 1], nrow = 1, ncol = 1)
# Jmat

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' Run STAN
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////

stan_data <- list(
  J = as.integer(2),
  Jmat = Jmat,
  N = N, 
  K = K, 
  X = X, 
  y = y,
  S = S,
  n_strata = n_strata,
  max_in_strata = max_in_strata,
  S_condensed = S_condensed
)

# Set path to model
stan_model <- cmdstan_model("SB_CondPoisson.stan")

out1 <- stan_model$sample(
  data = stan_data,
  chains = 2,
  parallel_chains = 2 
)

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' Output
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////


## 
draws_array <- out1$draws()

# Convert to data.frame (flattened, easier to use like extract())
draws_df <- posterior::as_draws_df(draws_array)
data.frame(head(draws_df))

# sick that seems to work
apply(draws_df %>% select(starts_with("beta_out")), 2, median)
coef(modelcpr1)

#
apply(draws_df %>% select(starts_with("q")), 2, summary)

# IT WORKS!!!!! WELL DONE :)