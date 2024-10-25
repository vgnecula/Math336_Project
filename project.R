# Required packages
library(tidyverse)
library(lubridate)
library(ggplot2)

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
analyze_japan_earthquakes <- function(start_date = Sys.Date() - years(1), 
                                      end_date = Sys.Date(), 
                                      min_magnitude = 4.5) {
  
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
japan_results <- analyze_japan_earthquakes()
