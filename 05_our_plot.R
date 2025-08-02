# Load libraries

library(dlnm)
library(sf)
library(lubridate)

# Load necessary files

load("input/daily_data.RData")
shapefile_bcn <- read_sf("input/shapefile_bcn.shp")
# WARNING: WE LOAD THE OUTPUT FROM MODEL 3 EXECUTED IN THE PAPER WITH THE 
# REAL MORTALITY DATA. IT CAN BE REPLACED FOR ANY OF THE OTHER OUTPUTS OBTAINED
# WITH REAL DATA OR OBTAINED WITH THE PREDICTED DATA FROM CODE 02_run_sbdlm.R
draws_df <- readRDS("draws_df.RDS")

dlnm_var <- list(
  var_prc = c(0.50, 0.90),
  var_fun = "bs",
  degree = 1,
  max_lag = 2,
  lagnk = 1,
  n_reg = 73)

# DATA PREPARATION

# Subset data to only summer months of 2007 to 2016

data <- subset(data, month(date) %in% 6:9)
data <- subset(data, year(date) %in% 2007:2016)

# Define the percentiles of temperature to be calculated 

# Create the temperature values used in the DLNM models 
GLOBAL_KNOTS = c(25, 35)
xCen = 25
x99 = 39
GLOBAL_BOUNDARY = c(16, 40)
x_temp <- 16:40


# Create a list with the cross-basis in each neighbourhood

cb <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  cb <- crossbasis(matrix(rep(x_temp, dlnm_var$max_lag + 1), 
                          ncol = dlnm_var$max_lag + 1),
                   argvar = list(fun = 'bs',
                                 degree = 2,
                                 knots = GLOBAL_KNOTS,
                                 Boundary.knots = GLOBAL_BOUNDARY),
                   arglag = list(fun = "ns",
                                 #degree = 2,
                                 knots = logknots(x = 2, 
                                                  nk = 1),
                                 intercept = F))
  
  
  return(cb)
  
})

#--------------------------------------------------------------
### C) MAP OF THE RISKS AT PERCENTILE 99
### (SAME AS FIGURE 4C IN THE PAPER)
#--------------------------------------------------------------

# Calculate directly the RR of the overall cumulative temperature-mortality 
# associations for all regions
beta_reg_all <- draws_df %>% select(starts_with("beta_out"))

rr <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  # Extract all the iterations of the coefficients of the crossbasis
  # beta_reg <- winbugs_res[,grepl(paste0("^beta\\[", i_reg,","), 
  #                                colnames(winbugs_res))]
  # 
  #   dim(beta_reg) # 1008 (sampled)   x  12 (coefficients)
  beta_reg <- beta_reg_all[,(i_reg * 8 - 7):(i_reg * 8)]
  
  # The RR in each temperature x is the sum of the product of x transformed
  # through the crossbasis function and the coefficients of the crossbasis
  rr <- apply(beta_reg, 1, function(x) {
    sapply(1:length(x_temp), function(i) cb[[i_reg]][i,] %*% x)
  })
  
  dim(rr)
  
  rr
  
})

# Create a function for centring the risk in each region
Center_RR <- function(f.rr, f.cen, f.temp){
  
  cen <- f.temp[which.min(abs(f.temp - f.cen))]
  rr <- apply(f.rr, 2, function(x) x - x[f.temp == cen])
  rr
  
}

# Center the relative risk in each region and extract the point estimate at 
# percentile 99
rr_plot <- sapply(1:dlnm_var$n_reg, function(i_reg) {
  x_plot <- x_temp
  cen_plot <- x_plot[x_plot == xCen]
  rr_plot <- Center_RR(f.rr = rr[[i_reg]], 
                       f.cen = cen_plot,
                       f.temp = x_plot)
  
  # Point estimate as the median of the values at percentile 99
  # it just happens to be 1008 iterations
  median(exp(rr_plot[x_plot == x99,]))
  
})
range(rr_plot)
sort(rr_plot)
# Set the minimum and maximum RR in the plots
rr_max <- round(max(rr_plot) * 10)/10
rr_min <- exp(-log(rr_max))

# Pallete of colours for the maps
pal <- leaflet::colorNumeric(palette = rev(
  c("#A90C38", "#C52A40", "#E24848", "#F16B61", "#F89183", "#FEB6A8", "#FEDAD3",
    "#FFFFFF", "#D3E5F2", "#A8CCE5", "#88B4D5", "#6D9CC3", "#5585B1", "#416F9C", 
    "#2E5A87")), domain = c(log(rr_min), log(rr_max)), reverse = FALSE)

# Plot the map with the risks

plot(shapefile_bcn$geometry, col = pal(log(rr_plot)), 
     main = "C) Map relative risks")

