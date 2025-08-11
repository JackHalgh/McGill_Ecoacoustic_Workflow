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

### BirdNET 

```
library(dplyr)
library(readr)

##### Part 1. Loading and combining BirdNET output #####

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

##### Part 2. Top 50 species #####

# Get list of all CSV files in the working directory
csv_files <- list.files(pattern = "\\.csv$")

# Read all CSV files into a list of data frames
data_list <- lapply(csv_files, read.csv)

# Optionally, name each element of the list after the file (without .csv extension)
names(data_list) <- tools::file_path_sans_ext(csv_files)

library(dplyr)
library(forcats)
library(stringr)

library(dplyr)
library(forcats)
library(stringr)

# Function to extract habitat from dataset name
extract_habitat <- function(name) {
  case_when(
    str_detect(name, "River") ~ "River",
    str_detect(name, "Maple_Beech") ~ "Maple_Beech",
    str_detect(name, "Oak") ~ "Oak",
    str_detect(name, "Wetland") ~ "Wetland",
    str_detect(name, "Beaver_Pond") ~ "Beaver_Pond",
    str_detect(name, "Lake_Shore") ~ "Lake_Shore",
    TRUE ~ "Unknown"
  )
}

# Loop through each data frame in the list
for (name in names(data_list)) {
  df <- data_list[[name]]
  
  # Get top 50 species and reorder
  top_50_species <- df %>%
    count(Common.name, sort = TRUE) %>%
    slice_head(n = 50) %>%
    pull(Common.name)
  
  top_data <- df %>%
    filter(Common.name %in% top_50_species) %>%
    mutate(Common.name = fct_infreq(Common.name))
  
  # Compute counts per species for labels
  species_counts <- top_data %>%
    count(Common.name) %>%
    mutate(label = paste0("n = ", n))
  
  # Extract habitat name from the dataset name
  habitat <- extract_habitat(name)
  species_counts <- species_counts %>%
    mutate(Habitat = habitat)
  
  # Write to CSV using the name of the dataset
  output_filename <- paste0(name, "_Top_50.csv")
  write.csv(species_counts, output_filename, row.names = FALSE)
}

##### Part 3. NMDS ####

library(dplyr)
library(readr)

# List all exported _Top_50.csv files in the directory
top50_files <- list.files(pattern = "_Top_50\\.csv$")

# Read and combine all files into a single data frame
combined_top50 <- top50_files %>%
  lapply(read_csv) %>%
  bind_rows()

# Write combined data to a single CSV
write_csv(combined_top50, "All_Sites_Top_50_Merged.csv")

library(dplyr)
library(tidyr)
library(vegan)

# Step 1: Summarize total counts per Habitat and Species (in case of duplicates)
summary_data <- combined_top50 %>%
  group_by(Habitat, Common.name) %>%
  summarise(total_n = sum(n), .groups = "drop")

# Step 2: Pivot to wide format (species matrix)
species_matrix <- summary_data %>%
  pivot_wider(names_from = Common.name, values_from = total_n, values_fill = list(total_n = 0))

# Step 3: Prepare matrix for NMDS
community_matrix <- as.data.frame(species_matrix[,-1])
rownames(community_matrix) <- species_matrix$Habitat

# Step 4: Run NMDS using Bray-Curtis distance
nmds <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

# Step 5: Plot NMDS results
plot(nmds, type = "t", main = "NMDS of Bird Communities by Habitat")


#### Part 4 Plotting NMDS ####

library(ggplot2)
library(dplyr)
library(ggrepel)  

# Extract NMDS site scores (habitats)
site_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  mutate(Label = rownames(.), Type = "Habitat")

# Extract NMDS species scores
species_scores <- as.data.frame(scores(nmds, display = "species")) %>%
  mutate(Label = rownames(.), Type = "Species")

# Combine both into one data frame
nmds_plot_data <- bind_rows(site_scores, species_scores)

# Separate species and habitat for plotting
species_only <- filter(nmds_plot_data, Type == "Species")
habitat_only <- filter(nmds_plot_data, Type == "Habitat")

ggplot() +
  # Plot habitat points
  geom_point(data = habitat_only, aes(x = NMDS1, y = NMDS2), color = "steelblue", size = 4, shape = 16) +
  
  # Label habitats (optional: just points are often clear enough)
  geom_text_repel(data = habitat_only, aes(x = NMDS1, y = NMDS2, label = Label),
                  size = 4, fontface = "bold", color = "steelblue") +
  
  # Label species (text only, black, repel to avoid overlap, no lines)
  geom_text_repel(data = species_only, aes(x = NMDS1 + 0.05, y = NMDS2 + 0.05, label = Label),
                  size = 3, color = "black", segment.color = NA) +
  
  theme_bw() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

```

### Plotting BirdNET call count and species richness ###

```
# Load necessary library
library(dplyr)

# Set working directory
setwd("C:/Users/jgreenhalgh/OneDrive - McGill University/Gault Data/All_Top_50")

# List all CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# Read each CSV into a list of data frames
data_list <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)

# Name each list element by the file name (without extension)
names(data_list) <- tools::file_path_sans_ext(basename(csv_files))

# Combine all data frames into one (handles mismatched columns)
combined_data <- bind_rows(data_list)

combined_data <- combined_data[-71138, ]

combined_data <- combined_data %>%
  select(-sensitivity, -min_conf)

# Check result
head(combined_data)

library(lubridate)

combined_data <- combined_data %>%
  mutate(
    # Extract the YYYYMMDD_HHMMSS part using regex
    datetime_str = sub("(\\d{8}_\\d{6}).*", "\\1", source_file),
    # Parse it into POSIXct datetime
    Datetime = ymd_hms(gsub("_", " ", datetime_str))
  ) %>%
  select(-datetime_str)  # Remove helper column

# Check result
head(combined_data)


library(dplyr)
library(ggplot2)
library(lubridate)

combined_data %>%
  mutate(Month = floor_date(Datetime, "month")) %>%
  group_by(Habitat, Month) %>%
  summarise(Count = n()) %>%
  ggplot(aes(x = Month, y = Count, color = Habitat)) +
  geom_line() +
  labs(x = "Month", y = "Number of Detections") +
  theme_bw()


#### Total call count ####

# Aggregate call count per Habitat, Date, and Hour (excluding NA Hour)
call_count_data <- combined_data %>%
  mutate(
    Date = as.Date(Datetime),
    Hour = hour(Datetime)  # Extract hour from Datetime
  ) %>%
  filter(!is.na(Hour)) %>%    # Remove rows with NA Hour
  group_by(Habitat, Date, Hour) %>%
  summarise(CallCount = n(), .groups = "drop") %>%  # Count calls instead of species richness
  mutate(HourLabel = sprintf("%02d:00", Hour))

# Plot heatmap of call count with formatted Hour labels
p_callcount <- ggplot(call_count_data, aes(y = Date, x = factor(HourLabel, levels = sprintf("%02d:00", 0:23)), fill = CallCount)) +
  geom_tile(color = "white") +
  facet_wrap(~ Habitat, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(
    y = "Date",
    x = "Hour of day",
    fill = "Call count"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

print(p_callcount)

ggsave("call_count_heatmap_with_time.jpeg", plot = p_callcount, device = "jpeg",
       width = 14, height = 10, units = "in", dpi = 300)


#### Species richness ####

# Aggregate species richness per Habitat, Date, and Hour (excluding NA Hour)
richness_data <- combined_data %>%
  mutate(
    Date = as.Date(Datetime),
    Hour = hour(Datetime)  # Extract hour from Datetime
  ) %>%
  filter(!is.na(Hour)) %>%    # Remove rows with NA Hour
  group_by(Habitat, Date, Hour) %>%
  summarise(SpeciesRichness = n_distinct(Scientific.name), .groups = "drop") %>%
  mutate(HourLabel = sprintf("%02d:00", Hour))

# Plot heatmap with formatted Hour labels
p_richness <- ggplot(richness_data, aes(y = Date, x = factor(HourLabel, levels = sprintf("%02d:00", 0:23)), fill = SpeciesRichness)) +
  geom_tile(color = "white") +
  facet_wrap(~ Habitat, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(
    y = "Date",
    x = "Hour of day",
    fill = "Species richness"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

print(p_richness)

ggsave("species_richness_heatmap_with_time.jpeg", plot = p_richness, device = "jpeg",
       width = 14, height = 10, units = "in", dpi = 300)

```




