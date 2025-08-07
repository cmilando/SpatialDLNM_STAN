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
library(shinystan)
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

saveRDS(dlnm_var, "dlnm_configuration.RDS")

# Set variables for trend and seasonality

# Subset data to only summer months of 2007 to 2016
data <- subset(data, month(date) %in% 6:9)
data <- subset(data, year(date) %in% 2007:2016)

data$year <- factor(lubridate::year(data$date))

# Create crossbasis for each region

# Ensure that the data is ordered by region to maintain alignment between 
# crossbasis and regions

if(is.unsorted(data$region)) {
  stop("data needs to be ordered by region for the next loop")}

list_X <- vector("list", dlnm_var$n_reg)

for(i_reg in 1:dlnm_var$n_reg) {
  
  temp <- subset(data, region == sprintf("%02i",i_reg), 
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
  
  ##
  list_X[[i_reg]] <- cb
}

dim(cb)        
tail(cb)


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

#' ===========================================================================
#' Normal Cond and Poisson



library(gnm)

df_local <- data_l[[39]]
cb_local <- list_X[[39]]

# keep just the non-empty strata
df_local_agg <- df_local %>%
  group_by(strata) %>%
  summarize(
    .groups = 'keep',
    sum_mort = sum(mort)
  ) %>% 
  mutate(keep = 1)

df_local <- left_join(df_local, df_local_agg)

# conditional quasi-poisson
m1.cond <- gnm(mort ~ cb_local, data = df_local,
          family = quasipoisson,
          eliminate = factor(strata),
          subset = keep == 1)
coef(m1.cond)

cp1 <- crosspred(cb_local, m1.cond)
plot(cp1, 'overall')

# quasi-poisson
m1.std <- gnm(mort ~ cb_local + factor(strata), data = df_local,
              family = quasipoisson,
              subset = keep == 1)

# but are the same as the base model
load("past_output/final_simsmatrix_model1_independent_casecrossover.RData")
i_reg <- 39
unique(substr(colnames(winbugs_res), 1, 5))
beta_reg1 <- winbugs_res[,grepl(paste0("^beta\\[", i_reg,","), 
                               colnames(winbugs_res))]

# but are the same as the base model
my_output <- readRDS("draws_df_no_spatial_39.RDS")
cols <- grep(paste0("beta"), names(my_output))
names(my_output)[cols]
beta_reg2 <- my_output[, cols]
head(beta_reg2)

my_output <- readRDS("draws_df_no_spatial.RDS")
i_reg <- 39
cols <- grep(paste0("beta\\[[0-9]+,",i_reg,"\\]"), names(my_output))
names(my_output)[cols]
beta_reg3 <- my_output[, cols]
head(beta_reg3)


# largely similar? 
cbind("paper" = apply(beta_reg1, 2, median), 
      "chad1-39" = apply(beta_reg2, 2, median),
      "chad1" = apply(beta_reg3, 2, median),
      "freq-condP" = coef(m1.cond), 
      "freq-P" = coef(m1.std)[2:13] )


# other things to try
# - maybe the bayes code isn't removing strata?
# - you also know that this works for simpler examples
# - 
