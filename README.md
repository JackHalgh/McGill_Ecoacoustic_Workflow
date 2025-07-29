# McGill_Ecoacoustic_Workflow

```
#By Jack A. Greenhalgh. June, 2025.
#Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.

#### Part 1. Loading, cleaning, and scaling data ####

# Load packages
library(corrplot)
library(caret)
library(stringr)
library(lubridate)
library(dplyr)
library(hms)
library(ggplot2)

# Load data
M001 <- read.csv("M001_Shore_alpha_acoustic_indices_results.csv")
M002 <- read.csv("M002_Myriophyllum_alpha_acoustic_indices_results.csv")
M003 <- read.csv("M003_Myriophyllum_alpha_acoustic_indices_results.csv")
M004 <- read.csv("M004_Pelagic_alpha_acoustic_indices_results.csv")
M005 <- read.csv("M005_Pelagic_alpha_acoustic_indices_results.csv")

# Merge all data sets into one
merged_data <- rbind(M001, M002, M003, M004, M005)
head(merged_data)

# Extract filename
Filename <- merged_data$filename

# Subset numeric columns from 2 to 61
numeric_data <- merged_data[, 2:61]

# Compute correlation matrix
cor_matrix <- cor(numeric_data, use = "complete.obs")

# Plot correlation matrix
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.cex = 0.7, tl.col = "black", addCoef.col = "black", number.cex = 0.5)
print(cor_matrix)

# Find indices of highly correlated variables (threshold > 0.8)
high_corr_indices <- findCorrelation(cor_matrix, cutoff = 0.8, names = TRUE)

# Remove them from the dataset
filtered_data <- numeric_data[, !colnames(numeric_data) %in% high_corr_indices]

head(filtered_data)

# Apply z-transformation
z_scaled_data <- as.data.frame(scale(filtered_data))

# View the result
head(z_scaled_data)
summary(z_scaled_data)

# Add filename back
z_scaled_data <- cbind(Filename, z_scaled_data)
head(z_scaled_data)

# Extract site information and store in a new column called site
z_scaled_data$Site <- str_extract(z_scaled_data$Filename, "^.*?_.*?(?=_)")
head(z_scaled_data[c("Filename", "Site")])

#Extract datetime
z_scaled_data <- z_scaled_data %>%
  mutate(
    # Extract datetime string (e.g., "20250513_120000") using regex
    datetime_str = str_extract(Filename, "\\d{8}_\\d{6}"),
    
    # Parse into POSIXct format (YYYYMMDD_HHMMSS)
    Datetime = as.POSIXct(datetime_str, format = "%Y%m%d_%H%M%S")
  )

# Round timestamps to the nearest 10
z_scaled_data <- z_scaled_data %>%
  mutate(
    Datetime = round_date(Datetime, unit = "10 minutes")
  )

head(z_scaled_data)

#### Part 2. Tapestry plots with 1 index ####

# Define base output directory
base_output_dir <- "C:/Users/jgreenhalgh/OneDrive - McGill University/Gault Data/May - July 2025"

# Min-max scaling function
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Split main dataframe by Site
site_dfs <- split(z_scaled_data, z_scaled_data$Site)

# Loop through each site
for (site_name in names(site_dfs)) {
  
  site_data <- site_dfs[[site_name]]
  
  # Skip if 'Datetime' column is missing
  if (!"Datetime" %in% names(site_data)) {
    warning(paste("Skipping site:", site_name, "- missing 'Datetime' column"))
    next
  }
  
  # Ensure Datetime is POSIXct
  if (!inherits(site_data$Datetime, "POSIXct")) {
    site_data$Datetime <- as.POSIXct(site_data$Datetime)
  }
  
  # Identify numeric columns to scale
  cols_to_exclude <- c("Datetime", "time_posix")
  numeric_vars <- names(site_data)[sapply(site_data, is.numeric)]
  cols_to_plot <- setdiff(numeric_vars, cols_to_exclude)
  
  # Apply min-max scaling
  site_data_scaled <- site_data
  site_data_scaled[cols_to_plot] <- lapply(site_data_scaled[cols_to_plot], min_max_scale)
  
  # Add date and time_of_day columns
  site_data_scaled <- site_data_scaled %>%
    mutate(
      date = as.Date(Datetime),
      time_of_day = as_hms(Datetime)
    )
  
  # Create subfolder for this site
  site_output_dir <- file.path(base_output_dir, site_name)
  dir.create(site_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Loop through each variable to plot
  for (var in cols_to_plot) {
    # Generate grayscale HEX codes
    site_data_scaled$HEX_codes <- rgb(
      red = site_data_scaled[[var]],
      green = site_data_scaled[[var]],
      blue = site_data_scaled[[var]],
      maxColorValue = 1
    )
    
    # Create heatmap
    p <- ggplot(site_data_scaled, aes(x = time_of_day, y = date)) +
      geom_tile(aes(fill = HEX_codes), color = "black") +
      scale_fill_identity() +
      labs(
        title = paste("Heatmap of", var, "at", site_name),
        x = "Time of Day",
        y = "Date"
      ) +
      scale_x_time(
        breaks = scales::breaks_width("2 hours"),
        labels = scales::time_format("%H:%M")
      ) +
      scale_y_date(
        date_breaks = "1 week",
        date_labels = "%b %d"
      ) +
      theme_bw()
    
    # Save to file
    ggsave(
      filename = file.path(site_output_dir, paste0("Heatmap_", var, ".png")),
      plot = p,
      width = 10,
      height = 10,
      dpi = 300
    )
  }
}

#### Part 3. Tapestry plots with 3 indices ####

# Load necessary packages
library(lubridate)
library(dplyr)
library(hms)
library(ggplot2)

# Extract timestamp string from Filename using regex
z_scaled_data$datetime_string <- sub(".*_(\\d{8}_\\d{6})\\.WAV", "\\1", z_scaled_data$Filename)

# Convert to POSIXct datetime
z_scaled_data$Datetime <- as.POSIXct(z_scaled_data$datetime_string, format = "%Y%m%d_%H%M%S")

# Round to nearest 10 minutes
z_scaled_data$Datetime <- round_date(z_scaled_data$Datetime, unit = "10 minutes")

# Extract date and time parts
z_scaled_data <- z_scaled_data %>%
  mutate(
    date = as.Date(Datetime),
    time = format(Datetime, "%H:%M:%S"),
    time_of_day = as_hms(Datetime),
    datetime_string = format(Datetime, "%Y%m%d_%H%M%S")  # regenerate based on rounded time
  )

# Identify numeric columns to scale (excluding Datetime)
cols_to_scale <- names(z_scaled_data)[sapply(z_scaled_data, is.numeric)]
cols_to_scale <- setdiff(cols_to_scale, "Datetime")  # exclude Datetime

# Min-max scaling function
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Apply scaling
z_scaled_data_scaled <- z_scaled_data
z_scaled_data_scaled[cols_to_scale] <- lapply(z_scaled_data_scaled[cols_to_scale], min_max_scale)

# Generate HEX codes using selected indices
z_scaled_data_scaled$HEX_codes <- rgb(
  red = z_scaled_data_scaled$EAS,
  green = z_scaled_data_scaled$HFC,
  blue = z_scaled_data_scaled$Ht,
  maxColorValue = 1
)

# Regenerate time_posix from rounded datetime_string (for safety)
z_scaled_data_scaled$time_posix <- as.POSIXct(z_scaled_data_scaled$datetime_string, format = "%Y%m%d_%H%M%S")

# Create the plot
ggplot(z_scaled_data_scaled, aes(x = time_of_day, y = date)) + 
  geom_tile(aes(fill = HEX_codes), color = "black") + 
  labs(x = "Time of day", y = "Date") + 
  scale_x_time(
    breaks = scales::breaks_width("2 hours"),
    labels = scales::time_format("%H:%M")
  ) + 
  scale_y_date(
    date_breaks = "1 week",
    date_labels = "%b %Y"
  ) + 
  theme_bw() + 
  scale_fill_identity()

```


### Automatic BirdNET extraction ###

```
setwd("D:/BirdNET output")

# Define column types
column_types <- cols(
  `Start (s)` = col_double(),
  `End (s)` = col_double(),
  `Low Freq (Hz)` = col_double(),
  `High Freq (Hz)` = col_double(),
  `Confidence` = col_double(),
  `Species` = col_character(),
  .default = col_character()
)

# Get folder names in working directory
folders <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

# Loop through each folder
for (folder in folders) {
  folder_path <- file.path(getwd(), folder)
  
  # List all .csv files in the folder
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Skip if no CSV files found
  if (length(csv_files) == 0) next
  
  # Read and combine all .csv files
  combined_data <- lapply(csv_files, function(file) {
    read_csv(file, col_types = column_types, show_col_types = FALSE) %>%
      mutate(source_file = basename(file))
  }) %>% bind_rows()
  
  # Create dynamic filename based on folder name
  output_filename <- paste0(folder, "_Combined_Results.csv")
  output_path <- file.path(folder_path, output_filename)
  
  # Write the combined data to CSV
  write_csv(combined_data, output_path)
  
  # Optional: load the result and drop unnecessary columns (adjust indices as needed)
  temp_data <- read.csv(output_path)
  if (ncol(temp_data) >= 18) {
    temp_data <- temp_data[, -c(10:18)]
  }
  assign(folder, temp_data)
  
  cat("Processed folder:", folder, "\n")
}
```
