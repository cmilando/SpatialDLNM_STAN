library(foreign) # ENABLES READING THE DATA FILE, WHICH IS A STATA FORMAT

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

# create S matrix 
getSmat <- function(strata_vector) {
  
  # strata_vector <- data$stratum 
  strata_matrix <- matrix(as.integer(strata_vector), 
                          nrow = length(strata_vector),
                          ncol = length(strata_vector), 
                          byrow = T)
  
  for(i in 1:length(strata_vector)) {
    strata_matrix[i, ] = 1*(strata_matrix[i, ] == as.integer(strata_vector[i]))
  }
  
  return(strata_matrix)
}



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


# nrows
N = as.integer(nrow(data))

# beta and the intercept
K = as.integer(2)

# include the intercept
X = as.matrix(data[, c('ozone10', 'temperature')])
head(X)

# outcome
y = as.integer(data$numdeaths)

#
S <- getSmat(data$stratum)

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
  
stan_data <- list(
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
stan_model <- cmdstan_model("CondPoisson_v0.stan")

out1 <- stan_model$sample(
  data = stan_data,
  chains = 1,
  parallel_chains = 4, 
)

## 
draws_array <- out1$draws()

# Convert to data.frame (flattened, easier to use like extract())
draws_df <- posterior::as_draws_df(draws_array)
head(draws_df)

# sick that seems to work
apply(draws_df %>% select(starts_with("beta")), 2, median)
coef(modelcpr1)
