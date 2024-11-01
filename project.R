# Required packages
library(tidyverse)  # This includes ggplot2, but we can load it explicitly to be sure
library(ggplot2)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# If any aren't installed, install them first:
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(lubridate)) install.packages("lubridate")
if (!require(rnaturalearth)) install.packages("rnaturalearth")
if (!require(rnaturalearthdata)) install.packages("rnaturalearthdata")
if (!require(sf)) install.packages("sf")

#------------------------------------------------------------------------------------------------------------------------------------------------
# FIT EXPONENTIAL USING MOM AND MLE
# input: numeric vector of interarrival times
# output: list containing rate (numeric), se (numeric), loglik (numeric)
#------------------------------------------------------------------------------------------------------------------------------------------------
fit_exponential_mom <- function(data) {
  rate <- 1/mean(data)
  se <- rate/sqrt(length(data))
  loglik <- sum(dexp(data, rate = rate, log = TRUE))
  return(list(rate = rate, se = se, loglik = loglik))
}

fit_exponential_mle <- function(data) {
  # For exponential, MLE is same as MOM
  rate <- 1/mean(data)
  se <- rate/sqrt(length(data))
  loglik <- sum(dexp(data, rate = rate, log = TRUE))
  return(list(rate = rate, se = se, loglik = loglik))
}

#------------------------------------------------------------------------------------------------------------------------------------------------
# FIT GAMMA USING MOM AND MLE
# input: numeric vector of interarrival times
# output: list containing shape (numeric), rate (numeric), loglik (numeric) for MOM
#         list containing shape (numeric), rate (numeric), se_shape (numeric), se_rate (numeric), loglik (numeric) for MLE
#------------------------------------------------------------------------------------------------------------------------------------------------
fit_gamma_mom <- function(data) {
  mean_x <- mean(data)
  var_x <- var(data)
  
  shape <- mean_x^2/var_x
  rate <- mean_x/var_x
  
  loglik <- sum(dgamma(data, shape = shape, rate = rate, log = TRUE))
  return(list(shape = shape, rate = rate, loglik = loglik))
}

fit_gamma_mle <- function(data) {
  # Use MOM estimates as initial values
  mean_x <- mean(data)
  var_x <- var(data)
  shape_init <- mean_x^2/var_x
  rate_init <- mean_x/var_x
  
  gamma_loglik <- function(params) {
    shape <- params[1]
    rate <- params[2]
    sum(dgamma(data, shape = shape, rate = rate, log = TRUE))
  }
  
  fit <- optim(c(shape_init, rate_init),
               fn = function(p) -gamma_loglik(p),
               method = "BFGS",
               hessian = TRUE)
  
  # Calculate standard errors
  if(!any(is.na(fit$hessian))) {
    fisher_info <- solve(fit$hessian)
    se <- sqrt(diag(fisher_info))
  } else {
    se <- c(NA, NA)
  }
  
  return(list(
    shape = fit$par[1],
    rate = fit$par[2],
    se_shape = se[1],
    se_rate = se[2],
    loglik = -fit$value
  ))
}

#------------------------------------------------------------------------------------------------------------------------------------------------
# FIT PLOT VISUALISATION - MOM OR MLE
# input: interarrival_times (numeric vector), method (string: "MOM" or "MLE")
# output: ggplot object containing histogram, fitted exponential and gamma curves, AIC values
#------------------------------------------------------------------------------------------------------------------------------------------------
create_method_plot <- function(interarrival_times, method) {
  df <- data.frame(interarrival = interarrival_times)
  x_range <- seq(0, max(interarrival_times), length.out = 200)
  
  # Fit distributions based on method
  if(method == "MOM") {
    exp_fit <- fit_exponential_mom(interarrival_times)
    gamma_fit <- fit_gamma_mom(interarrival_times)
    
    exp_dens <- dexp(x_range, rate = exp_fit$rate)
    gamma_dens <- dgamma(x_range, shape = gamma_fit$shape, rate = gamma_fit$rate)
    
    # Calculate AIC
    exp_aic <- -2 * exp_fit$loglik + 2 * 1
    gamma_aic <- -2 * gamma_fit$loglik + 2 * 2
  } else {  # MLE
    exp_fit <- fit_exponential_mle(interarrival_times)
    gamma_fit <- fit_gamma_mle(interarrival_times)
    
    exp_dens <- dexp(x_range, rate = exp_fit$rate)
    gamma_dens <- dgamma(x_range, shape = gamma_fit$shape, rate = gamma_fit$rate)
    
    # Calculate AIC
    exp_aic <- -2 * exp_fit$loglik + 2 * 1
    gamma_aic <- -2 * gamma_fit$loglik + 2 * 2
  }
  
  theoretical_data <- data.frame(
    x = rep(x_range, 2),
    density = c(exp_dens, gamma_dens),
    Distribution = rep(c("Exponential", "Gamma"), each = length(x_range))
  )
  
  # Create the plot
  ggplot() +
    geom_histogram(data = df,
                   aes(x = interarrival, y = ..density..),
                   bins = 30,
                   fill = "lightblue",
                   alpha = 0.7) +
    geom_line(data = theoretical_data,
              aes(x = x, y = density, color = Distribution),
              linewidth = 1) +
    scale_color_manual(values = c("Exponential" = "red", "Gamma" = "blue")) +
    labs(x = "Interarrival Time (hours)",
         y = "Density",
         title = paste("Distribution Fitting using", method),
         subtitle = paste0(method, " Parameter Estimates\n",
                           "Exp Rate: ", round(exp_fit$rate, 4), "\n",
                           "Gamma Shape: ", round(if(method == "MOM") gamma_fit$shape else gamma_fit$shape, 4),
                           ", Rate: ", round(if(method == "MOM") gamma_fit$rate else gamma_fit$rate, 4))) +
    theme_minimal() +
    annotate("text",
             x = max(interarrival_times) * 0.7,
             y = max(density(interarrival_times)$y) * 0.8,
             label = paste("Exponential AIC:", round(exp_aic, 2))) +
    annotate("text",
             x = max(interarrival_times) * 0.7,
             y = max(density(interarrival_times)$y) * 0.7,
             label = paste("Gamma AIC:", round(gamma_aic, 2)))
}



#------------------------------------------------------------------------------------------------------------------------------------------------
# VISUAL MAP FUNCTION
# input: dataframe containing columns: time (POSIXct), latitude (numeric), longitude (numeric), mag (numeric)
# output: ggplot object containing map of Japan with plotted earthquakes
#------------------------------------------------------------------------------------------------------------------------------------------------
create_japan_map <- function(df) {
  # Get Japan map data using natural earth
  japan_sf <- ne_countries(scale = "medium", 
                           country = "japan", 
                           returnclass = "sf")
  
  # Convert earthquake data to sf
  earthquakes_sf <- df %>%
    st_as_sf(coords = c("longitude", "latitude"),
             crs = 4326)
  
  # Create the map
  map_plot <- ggplot() +
    # Add Japan base map
    geom_sf(data = japan_sf, 
            fill = "grey80", 
            color = "grey40") +
    # Add earthquakes
    geom_sf(data = earthquakes_sf,
            aes(color = mag, size = mag),
            alpha = 0.6) +
    # Customize appearance
    scale_color_gradient(
      low = "yellow",
      high = "red",
      name = "Magnitude"
    ) +
    scale_size_continuous(
      range = c(2, 6),
      name = "Magnitude"
    ) +
    # Set map boundaries
    coord_sf(
      xlim = c(129, 146),
      ylim = c(30, 46)
    ) +
    # Labels and theme
    labs(
      title = "Earthquake Locations in Japan",
      subtitle = paste("Period:", 
                       format(min(df$time), "%Y-%m-%d"),
                       "to",
                       format(max(df$time), "%Y-%m-%d"))
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      panel.background = element_rect(fill = "aliceblue")
    )
  
  return(map_plot)
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# FETCHING EARTHQUAKES DATA
# input: start_date (string: "YYYY-MM-DD"), end_date (string: "YYYY-MM-DD"), min_magnitude (numeric)
# output: dataframe containing earthquake data (time, latitude, longitude, magnitude, etc.)
#------------------------------------------------------------------------------------------------------------------------------------------------
fetch_japan_earthquakes <- function(start_date, end_date, min_magnitude = 4.5) {
  # Japan region boundaries
  # More focused on main Japanese islands
  minlat <- 30    # Southern boundary
  maxlat <- 46    # Northern boundary
  minlon <- 129   # Western boundary
  maxlon <- 146   # Eastern boundary
  
  # Construct URL for USGS API
  base_url <- "https://earthquake.usgs.gov/fdsnws/event/1/query"
  
  url <- sprintf(paste0("%s?format=csv",
                        "&starttime=%s",
                        "&endtime=%s",
                        "&minmagnitude=%s",
                        "&minlatitude=%s",
                        "&maxlatitude=%s",
                        "&minlongitude=%s",
                        "&maxlongitude=%s"),
                 base_url, start_date, end_date, min_magnitude,
                 minlat, maxlat, minlon, maxlon)
  
  # Fetch data
  data <- read.csv(url)
  return(data)
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# MAIN ANALYSIS LOGIC
# input: start_date (string: "YYYY-MM-DD"), end_date (string: "YYYY-MM-DD"), min_magnitude (numeric)
# output: list containing:
#         - data (dataframe of raw earthquake data)
#         - interarrival_times (numeric vector)
#         - basic_stats (list of summary statistics)
#         - map (ggplot object)
#         - mom_results (list of MOM fitting results)
#         - mle_results (list of MLE fitting results)
#------------------------------------------------------------------------------------------------------------------------------------------------
analyze_japan_earthquakes <- function(start_date, end_date, min_magnitude) {
  # Fetch data
  cat("Fetching earthquake data for Japan...\n")
  df <- fetch_japan_earthquakes(start_date, end_date, min_magnitude)
  
  if(nrow(df) < 2) {
    stop("Insufficient earthquake data for analysis")
  }
  
  # Convert time and sort
  df$time <- as.POSIXct(df$time, format = "%Y-%m-%dT%H:%M:%S")
  df <- df[order(df$time), ]
  
  # Calculate interarrival times
  times <- diff(df$time)
  interarrival_times <- as.numeric(times, units = "hours")
  
  # Create map plot first
  map_plot <- create_japan_map(df)
  print(map_plot)  # Display the map
  
  # Basic statistics
  stats <- list(
    total_earthquakes = nrow(df),
    mean_magnitude = mean(df$mag),
    max_magnitude = max(df$mag),
    mean_interarrival = mean(interarrival_times),
    median_interarrival = median(interarrival_times),
    sd_interarrival = sd(interarrival_times)
  )
  
  # Create both MOM and MLE plots
  mom_plot <- create_method_plot(interarrival_times, "MOM")
  mle_plot <- create_method_plot(interarrival_times, "MLE")
  
  # Fit all distributions
  mom_exp <- fit_exponential_mom(interarrival_times)
  mom_gamma <- fit_gamma_mom(interarrival_times)
  mle_exp <- fit_exponential_mle(interarrival_times)
  mle_gamma <- fit_gamma_mle(interarrival_times)
  
  # Print basic statistics
  cat("\nBasic Statistics:\n")
  cat("=================\n")
  cat("Total earthquakes:", stats$total_earthquakes, "\n")
  cat("Mean magnitude:", round(stats$mean_magnitude, 2), "\n")
  cat("Maximum magnitude:", stats$max_magnitude, "\n")
  cat("Mean interarrival time:", round(stats$mean_interarrival, 2), "hours\n")
  cat("Median interarrival time:", round(stats$median_interarrival, 2), "hours\n")
  
  # Print results for both methods
  cat("\nMethod of Moments (MOM) Results:\n")
  cat("================================\n")
  cat("Exponential Rate:", round(mom_exp$rate, 4), "\n")
  cat("Gamma Shape:", round(mom_gamma$shape, 4), "\n")
  cat("Gamma Rate:", round(mom_gamma$rate, 4), "\n")
  
  cat("\nMaximum Likelihood Estimation (MLE) Results:\n")
  cat("==========================================\n")
  cat("Exponential Rate:", round(mle_exp$rate, 4), "\n")
  cat("Gamma Shape:", round(mle_gamma$shape, 4), "\n")
  cat("Gamma Rate:", round(mle_gamma$rate, 4), "\n")
  
  # Display distribution plots
  print(mom_plot)
  print(mle_plot)
  
  # Return results
  return(list(
    data = df,
    interarrival_times = interarrival_times,
    basic_stats = stats,
    map = map_plot,  # Now including the map in the returned results
    mom_results = list(
      exponential = mom_exp,
      gamma = mom_gamma,
      plot = mom_plot
    ),
    mle_results = list(
      exponential = mle_exp,
      gamma = mle_gamma,
      plot = mle_plot
    )
  ))
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# "MAIN"
#------------------------------------------------------------------------------------------------------------------------------------------------
japan_results <- analyze_japan_earthquakes("2023-10-18", "2024-10-19", 4.5)


