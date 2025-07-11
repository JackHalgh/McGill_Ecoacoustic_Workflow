# McGill_Ecoacoustic_Workflow

```
#By Jack A. Greenhalgh. June, 2025.
#Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.

#### Part 1. Loading, cleaning, and scaling data ####

M001 <- read.csv("M001_Shore_alpha_acoustic_indices_results.csv")
head(M001)

library(corrplot)
library(caret)

# Extract filename

Filename <- M001$filename

# Subset numeric columns from 2 to 61
numeric_data <- M001[, 2:61]

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

#### Part 2. Tapestry plots with 1 index ####

# Load necessary packages
library(lubridate)
library(dplyr)
library(hms)
library(ggplot2)

# Set working directory for saving plots (use forward slashes)
output_dir <- "C:/Users/jgreenhalgh/OneDrive - McGill University/Gault Data/May - July 2025/M001"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)


# Identify numeric columns to scale (excluding time-related ones)
cols_to_exclude <- c("Datetime", "time_posix")
numeric_vars <- names(z_scaled_data)[sapply(z_scaled_data, is.numeric)]
cols_to_plot <- setdiff(numeric_vars, cols_to_exclude)

# Apply min-max scaling
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
z_scaled_data_scaled <- z_scaled_data
z_scaled_data_scaled[cols_to_plot] <- lapply(z_scaled_data_scaled[cols_to_plot], min_max_scale)

# Ensure time columns exist
if (!"date" %in% colnames(z_scaled_data_scaled)) {
  z_scaled_data_scaled <- z_scaled_data_scaled %>%
    mutate(
      date = as.Date(Datetime),
      time_of_day = as_hms(Datetime)
    )
}

# Loop through each variable and generate + save plots
for (var in cols_to_plot) {
  # Generate grayscale HEX color based on the current variable
  z_scaled_data_scaled$HEX_codes <- rgb(
    red = z_scaled_data_scaled[[var]],
    green = z_scaled_data_scaled[[var]],
    blue = z_scaled_data_scaled[[var]],
    maxColorValue = 1
  )
  
  # Create the plot
  p <- ggplot(z_scaled_data_scaled, aes(x = time_of_day, y = date)) +
    geom_tile(aes(fill = HEX_codes), color = "black") +
    scale_fill_identity() +
    labs(
      title = paste("Heatmap of", var),
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
  
  # Save the plot as PNG
  ggsave(
    filename = paste0("Heatmap_", var, ".png"),
    plot = p,
    width = 10,
    height = 10,
    dpi = 300
  )
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
  red = z_scaled_data_scaled$BioEnergy,
  green = z_scaled_data_scaled$BioEnergy,
  blue = z_scaled_data_scaled$BioEnergy,
  maxColorValue = 1
)

# Regenerate time_posix from rounded datetime_string (for safety)
z_scaled_data_scaled$time_posix <- as.POSIXct(z_scaled_data_scaled$datetime_string, format = "%Y%m%d_%H%M%S")

# 🖼️ Create the plot
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
