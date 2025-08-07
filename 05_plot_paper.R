# Load libraries

library(dlnm)
library(sf)
library(lubridate)
library(leaflet)
library(future)

# Load necessary files

load("input/daily_data.RData")
shapefile_bcn <- read_sf("input/shapefile_bcn.shp")
# WARNING: WE LOAD THE OUTPUT FROM MODEL 3 EXECUTED IN THE PAPER WITH THE 
# REAL MORTALITY DATA. IT CAN BE REPLACED FOR ANY OF THE OTHER OUTPUTS OBTAINED
# WITH REAL DATA OR OBTAINED WITH THE PREDICTED DATA FROM CODE 02_run_sbdlm.R
load("past_output/final_simsmatrix_model3_spatial_casecrossover.RData")
load("past_output/final_simsmatrix_model1_independent_casecrossover.RData")

dlnm_var <- list(
  var_prc = c(0.50, 0.90),
  var_fun = "ns",
  max_lag = 8,
  lagnk = 2,
  n_reg = 73,
  n_coef = 12)

# DATA PREPARATION

# Subset data to only summer months of 2007 to 2016

data <- subset(data, month(date) %in% 6:9)
data <- subset(data, year(date) %in% 2007:2016)

# Define the percentiles of temperature to be calculated 
percentiles <- c(seq(0, 1, by = 0.1), 
                 seq(2, 98, by = 1), 
                 seq(99, 100, by = 0.1)) /100

# Create the temperature values used in the DLNM models 
list_temp <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  label_reg <- formatC(i_reg, width = 2, flag = "0")
  temp <- subset(data, region == label_reg, 
                 select = c("temp", paste0("lag", 1:dlnm_var$max_lag)))
  temp_knots <- quantile(temp[["temp"]], dlnm_var$var_prc, na.rm = TRUE)
  temp_boundary <- range(temp, na.rm = TRUE)
  x_temp <- quantile(temp[["temp"]], percentiles, na.rm = TRUE)
  
  return(list(temp_knots = temp_knots, # knots of the exposure-response function  
              temp_boundary = temp_boundary, # range of temperatures
              x_temp = x_temp)) # temperatures in which risk are calculated
  
})

temp_knots <- sapply(list_temp, function(x) x[["temp_knots"]])
temp_boundary <- sapply(list_temp, function(x) x[["temp_boundary"]])
x_temp <- sapply(list_temp, function(x) x[["x_temp"]])
rm(list_temp)

# Create a list with the basis for exposure and the basis for lags 
# in each neighbourhood

basis_all <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  f.temp_knots <- temp_knots[,i_reg]
  f.temp_boundary <- temp_boundary[,i_reg]
  f.x_temp <- x_temp[,i_reg]
  
  # basis temperatures
  Q <- onebasis(f.x_temp, fun = dlnm_var$var_fun, knots = f.temp_knots, 
                Boundary.knots = f.temp_boundary)
  
  # basis lags
  C <- onebasis(0:dlnm_var$max_lag, fun = "ns", 
                knots = logknots(dlnm_var$max_lag, dlnm_var$lagnk), 
                intercept = TRUE)
  
  return(list(basis_exp = Q, basis_lag = C))
  
})

# Create a list with the cross-basis in each neighbourhood

cb <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  cb <- crossbasis(matrix(rep(x_temp[,i_reg], dlnm_var$max_lag + 1), 
                          ncol = dlnm_var$max_lag + 1),
                   argvar = list(fun = dlnm_var$var_fun,
                                 knots = temp_knots[,i_reg],
                                 Boundary.knots = temp_boundary[,i_reg]),
                   arglag = list(fun = "ns",
                                 knots = logknots(dlnm_var$max_lag, dlnm_var$lagnk),
                                 intercept = TRUE))
  
  return(cb)
  
})

Center_RR <- function(f.rr, f.cen, f.temp){
  
  cen <- f.temp[which.min(abs(f.temp - f.cen))]
  rr <- apply(f.rr, 2, function(x) x - x[f.temp == cen])
  rr
  
}

#--------------------------------------------------------------
### C) MAP OF THE RISKS AT PERCENTILE 99
### (SAME AS FIGURE 4C IN THE PAPER)
#--------------------------------------------------------------

# my_output <- readRDS("draws_df_no_spatial.RDS")

# Calculate directly the RR of the overall cumulative temperature-mortality 
# associations for all regions
rr <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  # Extract all the iterations of the coefficients of the crossbasis
  beta_reg <- winbugs_res[,
                          grepl(paste0("^beta\\[", i_reg,","), colnames(winbugs_res))]

  # cols <- grep(paste0("beta_star\\[[0-9]+,",i_reg,"\\]"), names(my_output))
  # cols <- grep(paste0("beta\\[[0-9]+,",i_reg,"\\]"), names(my_output))
  # names(my_output)[cols]
  # beta_reg <- my_output[, cols]
  # head(beta_reg)
  
  # The RR in each temperature x is the sum of the product of x transformed
  # through the crossbasis function and the coefficients of the crossbasis
  rr <- apply(beta_reg, 1, function(x) {
    sapply(1:length(x_temp[,i_reg]), function(i) cb[[i_reg]][i,] %*% x)
  })
  
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
  x_plot <- x_temp[, i_reg]
  cen_plot <- x_plot[percentiles == 0.3]
  rr_plot <- Center_RR(f.rr = rr[[i_reg]], 
                       f.cen = cen_plot,
                       f.temp = x_plot)
  
  # Point estimate as the median of the values at percentile 99
  xx <- median(exp(rr_plot[percentiles == 0.99,]))
  
  # *****
  # uh i think this was added
  min(2, max(xx, 0.5))
  # *****

})

summary(rr_plot)

# Set the minimum and maximum RR in the plots
rr_max <- 2
rr_min <- exp(-log(rr_max))

# Pallete of colours for the maps
pal <- colorNumeric(palette = rev(
  c("#A90C38", "#C52A40", "#E24848", "#F16B61", "#F89183", "#FEB6A8", "#FEDAD3",
    "#FFFFFF", "#D3E5F2", "#A8CCE5", "#88B4D5", "#6D9CC3", "#5585B1", "#416F9C", 
    "#2E5A87")), domain = c(log(rr_min), log(rr_max)), reverse = FALSE)

# Plot the map with the risks
# 
# plot(shapefile_bcn$geometry, col = pal(log(rr_plot)), 
#      main = "C) Map relative risks")


# Compute centroids
centroids <- sf::st_centroid(shapefile_bcn)

# Extract coordinates of centroids
coords <- sf::st_coordinates(centroids)

# Plot the map
plot(shapefile_bcn$geometry, col = pal(log(rr_plot)), 
     main = "C) Map relative risks")

# Add labels
text(coords[, 1], coords[, 2], labels = 1:nrow(shapefile_bcn), cex = 0.7)

which(Jmat[22,] ==1 )

rr_breaks <- seq(0.5, 2, length.out = 100)
log_breaks <- log(rr_breaks)
pal_vals <- pal(length(log_breaks))  # assuming pal is a function like colorRampPalette()
pal_vals
library(fields)
image.plot(
  legend.only = TRUE,
  zlim = log(c(rr_min, rr_max)),             # color scale in log space
  col = pal(seq(log(rr_min), log(rr_max), length.out = 100)),
  axis.args = list(
    at = log(c(0.50, 0.59, 0.71, 0.84, 1.00, 1.19, 1.41, 1.68, 2.00)),  # tick positions (log space)
    labels = c(0.50, 0.59, 0.71, 0.84, 1.00, 1.19, 1.41, 1.68, 2.00)    # tick labels (original RR scale)
  ),
  legend.args = list(
    text = "Relative Risk", side = 3, line = 1, cex = 0.8
  )
)

