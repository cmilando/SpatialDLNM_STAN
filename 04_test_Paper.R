################################################################################
# IMPORTANT: THE MORTALITY DATA provided are simulated FROM MODEL 4 IN THE
# PAPER (SB-DLNM WITH TIME-SERIES DESIGN) due to confidentiality restrictions of 
# the original data set. Original data may be requested to the data owner (The 
# Barcelona Public Health Agency) who should share them in similar terms than 
# those applying to this study. The supplied mortality, with exactly the same 
# structure than the original data set, allows reproducing the code provided so 
# that interested readers may run it as an example of use. However, WinBUGS 
# results from our original analyses are additionally supplied (input/result_paper
# folder) so that readers should be able to reproduce exactly our analyses after 
# those WinBUGS calls (model summaries and plots).
################################################################################

################################################################################
# In this R project we implement Bayesian and Spatial Bayesian Distributed 
# Lag Non-Linear Models (B-DLNM and SB-DLNM) for the case Study of short-term 
# associations between temperature and mortality in the city of Barcelona.
################################################################################

################################################################################
# CODE 1: DATA PREPARATION
# Data preparation: Transforming time-series temperature and mortality data sets
# for B-DLNMs and SB-DLNMs.
################################################################################

# Load libraries

library(sf)
library(spdep)
library(lubridate)
library(dlnm)
library(gnm)
library(splines)
library(rstan)
library(cmdstanr)
library(tidyverse)
library(foreign)
# Load data

load("input/daily_data.RData")

# Set variables defining the dlnm model

dlnm_var <- list(
  var_prc = c(0.50, 0.90),
  var_fun = "ns",
  max_lag = 8,
  lagnk = 2,
  n_reg = 73,
  n_coef = 12)

# Set variables for trend and seasonality

df_seas <- 4
df_trend_10years <- 1 # 1 df every 10 years to control for long-term trends
df_trend <- round(length(2007:2016) / df_trend_10years / 10) # Here we assume the time period for all regions is the same
rm(df_trend_10years)

dlnm_var$df_seas <- df_seas
dlnm_var$df_trend <- df_trend

# Subset data to only summer months of 2007 to 2016

data <- subset(data, month(date) %in% 6:9)
data <- subset(data, year(date) %in% 2007:2016)

# Create crossbasis for each region

# Ensure that the data is ordered by region to maintain alignment between 
# crossbasis and regions

if(is.unsorted(data$region)) {
  stop("data needs to be ordered by region for the next loop")}

list_cb <- lapply(formatC(1:dlnm_var$n_reg, width = 2, flag = "0"), 
      function(i_reg) {
        
        temp <- subset(data, region == i_reg, 
                       select = c("temp", paste0("lag", 1:dlnm_var$max_lag)))
        
        temp_knots <- quantile(temp$temp, dlnm_var$var_prc, na.rm = TRUE)
        temp_boundary <- range(temp, na.rm = TRUE)
        
        cb <- crossbasis(temp,
                         argvar = list(fun = dlnm_var$var_fun,
                                       knots = temp_knots,
                                       Boundary.knots = temp_boundary),
                         arglag = list(fun = "ns",
                                       knots = logknots(dlnm_var$max_lag, 
                                                        dlnm_var$lagnk),
                                       intercept = TRUE))
        
        cb
        
      })


# PREPARE DATA FOR THE CASE-CROSSOVER DESIGN

# Create strata for the case-crossover
# (neighborhood - year - month - day of week)
data$strata <- paste(data$region, 
                         year(data$date), 
                         formatC(month(data$date), width = 2, flag = "0"),
                         wday(data$date, week_start = 1),
                         sep = ":")

# make them both lists
data_l <- split(data, f = data$region)

#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
#' SETUP inputs
#' ============================================================================
#' #' /////////////////////////////////////////////////////////////////////////

# n regions
J = as.integer(length(data_l))

# nrows
N = as.integer(nrow(data_l[[1]]))

# beta values, withouth the intercept
K = as.integer(ncol(list_cb[[1]]))

# include the intercept
X = array(dim = c(dim(list_cb[[1]]), J))
for(j in 1:J) X[,,j] = as.matrix(list_cb[[j]])

# outcome in two regions
y = array(dim = c(nrow(list_cb[[1]]), J))
for(j in 1:J) y[,j] = data_l[[j]]$mort

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
S <- getSmat(factor(data_l[[1]]$strata), include_self = T)
head(S)

# get strata vars
n_strata <- as.integer(length(unique(data_l[[1]]$strata))) # 72, cool
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
stratum_id = as.integer(factor(data_l[[1]]$strata))
stratum_id

## *** THIS BECOMES THE J MATRIX
# ok now do the same as include self but with J matrix
# start simple and you can generalize later
# this is an adjacency matrix that DOES NOT include itself
shapefile_bcn <- read_sf("input/shapefile_bcn.shp")

# Generate a list of spatial structure from the shapefile for use in WinBUGS.
list_neig <- nb2listw(poly2nb(shapefile_bcn))

neighbors <- lapply(list_neig$neighbours,c)
Jmat <- matrix(0, nrow = J, ncol = J)

for(j in 1:J) {
    Jmat[j, neighbors[[j]]] <- 1
}
head(Jmat)
## ******



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
stan_model <- cmdstan_model("SB_CondPoisson.stan")

out1 <- stan_model$sample(
  data = stan_data,
  chains = 1,
  parallel_chains = 1 
)


