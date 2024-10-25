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

# Function to fit exponential distribution using MLE
fit_exponential <- function(data) {
  # MLE for exponential is 1/mean
  rate <- 1/mean(data)
  # Standard error calculation
  se <- rate/sqrt(length(data))
  return(list(estimate = rate, sd = se))
}

# Function to fit gamma distribution using MLE
fit_gamma <- function(data) {
  # Initial guess using method of moments
  mean_x <- mean(data)
  var_x <- var(data)
  
  # Initial shape and rate parameters
  shape_init <- mean_x^2/var_x
  rate_init <- mean_x/var_x
  
  # Log-likelihood function for gamma distribution
  gamma_loglik <- function(params) {
    shape <- params[1]
    rate <- params[2]
    sum(dgamma(data, shape = shape, rate = rate, log = TRUE))
  }
  
  # Optimize using optim
  fit <- optim(c(shape_init, rate_init), 
               fn = function(p) -gamma_loglik(p), 
               method = "BFGS")
  
  return(list(shape = fit$par[1], rate = fit$par[2]))
}

# Main function to fit and compare distributions
fit_and_compare_distributions <- function(interarrival_times) {
  # Fit exponential distribution
  exp_fit <- fit_exponential(interarrival_times)
  
  # Fit gamma distribution
  gamma_fit <- fit_gamma(interarrival_times)
  
  # Calculate log-likelihoods
  exp_loglik <- sum(dexp(interarrival_times, 
                         rate = exp_fit$estimate, 
                         log = TRUE))
  
  gamma_loglik <- sum(dgamma(interarrival_times, 
                             shape = gamma_fit$shape, 
                             rate = gamma_fit$rate, 
                             log = TRUE))
  
  # Calculate AIC
  n <- length(interarrival_times)
  exp_aic <- -2 * exp_loglik + 2 * 1  # Exponential has 1 parameter
  gamma_aic <- -2 * gamma_loglik + 2 * 2  # Gamma has 2 parameters
  
  # Return results
  return(list(
    exponential = exp_fit,
    gamma = gamma_fit,
    exp_aic = exp_aic,
    gamma_aic = gamma_aic
  ))
}



# Function to create visualization of the fits
create_plots <- function(interarrival_times, results) {
  # Create data frame for plotting
  df <- data.frame(
    interarrival = interarrival_times
  )
  
  # Generate points for theoretical distributions
  x_range <- seq(0, max(interarrival_times), length.out = 200)
  
  # Create theoretical distribution data
  theoretical_data <- data.frame(
    x = rep(x_range, 2),
    density = c(
      dexp(x_range, rate = results$exponential$estimate),
      dgamma(x_range, shape = results$gamma$shape, rate = results$gamma$rate)
    ),
    Distribution = rep(c("Exponential", "Gamma"), each = length(x_range))
  )
  
  # Create the plot
  ggplot() +
    # Add histogram of actual data
    geom_histogram(data = df, 
                   aes(x = interarrival, y = ..density..), 
                   bins = 30, 
                   fill = "lightblue", 
                   alpha = 0.7) +
    # Add theoretical distributions
    geom_line(data = theoretical_data,
              aes(x = x, y = density, color = Distribution),
              linewidth = 1) +
    # Customize appearance
    scale_color_manual(values = c("Exponential" = "red", "Gamma" = "blue")) +
    labs(x = "Interarrival Time (hours)",
         y = "Density") +
    theme_minimal() +
    # Add AIC values as annotations
    annotate("text", 
             x = max(interarrival_times) * 0.7,
             y = max(density(interarrival_times)$y) * 0.8,
             label = paste("Exponential AIC:", 
                           round(results$exp_aic, 2))) +
    annotate("text",
             x = max(interarrival_times) * 0.7,
             y = max(density(interarrival_times)$y) * 0.7,
             label = paste("Gamma AIC:", 
                           round(results$gamma_aic, 2)))
}



# Function to create Japan map with sf
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

# Function to fetch earthquake data specifically for Japan
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


# Modified analysis function to include new map
analyze_japan_earthquakes <- function(start_date, 
                                      end_date, 
                                      min_magnitude) {
  
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
  
  # Basic statistics
  stats <- list(
    total_earthquakes = nrow(df),
    mean_magnitude = mean(df$mag),
    max_magnitude = max(df$mag),
    mean_interarrival = mean(interarrival_times),
    median_interarrival = median(interarrival_times),
    sd_interarrival = sd(interarrival_times)
  )
  
  # Model comparison
  results <- fit_and_compare_distributions(interarrival_times)
  
  # Create visualizations
  # 1. New map plot with sf
  map_plot <- create_japan_map(df)
  
  # 2. Time series of magnitudes
  time_plot <- ggplot(df, aes(x = time, y = mag)) +
    geom_point(aes(color = mag), alpha = 0.6) +
    scale_color_gradient(low = "yellow", high = "red", name = "Magnitude") +
    labs(title = "Earthquake Magnitudes Over Time",
         x = "Date", y = "Magnitude") +
    theme_minimal()
  
  # 3. Distribution fit plot
  dist_plot <- create_plots(interarrival_times, results) +
    labs(title = "Earthquake Interarrival Times in Japan",
         subtitle = paste("Minimum magnitude:", min_magnitude))
  
  # 4. Monthly frequency plot
  monthly_counts <- df %>%
    mutate(month = floor_date(time, "month")) %>%
    count(month)
  
  monthly_plot <- ggplot(monthly_counts, aes(x = month, y = n)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = "Monthly Earthquake Frequency",
         x = "Month", y = "Number of Earthquakes") +
    theme_minimal()
  
  # Print results
  cat("\nJapan Earthquake Analysis Summary:\n")
  cat("================================\n")
  cat("Time Period:", format(start_date), "to", format(end_date), "\n")
  cat("Minimum magnitude considered:", min_magnitude, "\n")
  cat("\nBasic Statistics:\n")
  cat("Total earthquakes:", stats$total_earthquakes, "\n")
  cat("Mean magnitude:", round(stats$mean_magnitude, 2), "\n")
  cat("Maximum magnitude:", stats$max_magnitude, "\n")
  cat("Mean interarrival time:", round(stats$mean_interarrival, 2), "hours\n")
  cat("Median interarrival time:", round(stats$median_interarrival, 2), "hours\n")
  
  # Model comparison results
  cat("\nModel Comparison Results:\n")
  cat("--------------------------------\n")
  cat("Exponential Distribution:\n")
  cat("Rate:", round(results$exponential$estimate, 4), "\n")
  cat("Standard Error:", round(results$exponential$sd, 4), "\n")
  
  cat("\nGamma Distribution:\n")
  cat("Shape:", round(results$gamma$shape, 4), "\n")
  cat("Rate:", round(results$gamma$rate, 4), "\n")
  
  cat("\nAIC Values:\n")
  cat("Exponential:", round(results$exp_aic, 2), "\n")
  cat("Gamma:", round(results$gamma_aic, 2), "\n")
  
  # Display plots
  print(map_plot)
  print(time_plot)
  print(dist_plot)
  print(monthly_plot)
  
  # Return results
  return(list(
    data = df,
    interarrival_times = interarrival_times,
    basic_stats = stats,
    model_results = results,
    plots = list(
      map = map_plot,
      time_series = time_plot,
      distribution = dist_plot,
      monthly = monthly_plot
    )
  ))
}

# Run the analysis
japan_results <- analyze_japan_earthquakes("2023-10-18", "2024-10-19", 4.5)