library(dlnm)
library(gnm)
library(splines)
library(rstan)
library(cmdstanr)
library(tidyverse)
library(foreign) # ENABLES READING THE DATA FILE, WHICH IS A STATA FORMAT
library(shinystan)
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

# subset for demonstration
data <- subset(data, dow %in% c("Monday", "Tuesday") & 
                 year %in% c("2002", "2003", "2004"))
data$stratum <- interaction(data$year, data$month, data$dow)
length(levels(data$stratum))
data <- data[order(data$date), ]

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' SETUP inputs
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////

# n regions
J = as.integer(3)

# nrows
N = as.integer(nrow(data))

# beta values, withouth the intercept
K = as.integer(2)

# include the intercept
X1 = as.matrix(data[, c('ozone10', 'temperature')])
X = array(dim = c(dim(X1), 3))
X[,,1] = X1
X[,,2] = X1
X[,,3] = X1

# outcome in two regions
y1 = as.integer(data$numdeaths)
y = cbind(y1, y1, y1)

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

#
stratum_id = as.integer(data$stratum)
stratum_id

# ok now do the same as include self but with J matrix
# start simple and you can generalize later
# this is an adjacency matrix that DOES NOT include itself
Jmat <- matrix(as.integer(c(0, 1, 1, 1, 0, 1, 1, 1, 0)), nrow = J)
Jmat

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' Run STAN
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////

stan_data <- list(
  J = J,
  Jmat = Jmat,
  N = N, 
  K = K, 
  X = X, 
  y = y,
  S = S,
  n_strata = n_strata,
  max_in_strata = max_in_strata,
  S_condensed = S_condensed,
  stratum_id = stratum_id
)

# Set path to model
# https://mc-stan.org/learn-stan/case-studies/reduce_sum_tutorial.html
stan_model <- cmdstan_model("SB_CondPoisson_v3b.stan",
                            cpp_options = list(stan_threads = TRUE))

out1 <- stan_model$sample(
  data = stan_data,
  chains = 2,
  parallel_chains = 2,
  max_treedepth = 8,
  threads_per_chain = 2,
  refresh = 100
)

# 165 for tree-depth 8

# shinystan::launch_shinystan(out1)

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' Output
#' >> ADD POST-PROCESSING HERE !
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////

## 
draws_array <- out1$draws()

# Convert to data.frame (flattened, easier to use like extract())
draws_df <- posterior::as_draws_df(draws_array)
dim(draws_df)
head(draws_df)
data.frame(head(draws_df))

# apply(draws_df %>% select(starts_with("beta")), 2, median)

# Calculate
# I think this 

beta11 <- draws_df$`mu[1]` + draws_df$`beta_star[1,1]`
beta12 <- draws_df$`mu[1]` + draws_df$`beta_star[1,2]`
beta13 <- draws_df$`mu[1]` + draws_df$`beta_star[1,3]`

beta21 <- draws_df$`mu[2]` + draws_df$`beta_star[2,1]`
beta22 <- draws_df$`mu[2]` + draws_df$`beta_star[2,2]`
beta32 <- draws_df$`mu[2]` + draws_df$`beta_star[2,3]`

median(beta11)
median(beta21)

median(beta12)
median(beta22)

median(beta13)
median(beta32)


# comparison
modelcpr1 <- gnm(numdeaths ~ ozone10 + temperature, 
                 data = data, family = quasipoisson, 
                 eliminate = factor(stratum))
coef(modelcpr1)

# q -- hmm this should be 0.9 almost certainly
apply(draws_df %>% select(starts_with("q")), 2, summary)

# IT WORKS!!!!! WELL DONE :)


