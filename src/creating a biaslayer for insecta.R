# --------------------------------------------------------------------------------------
# Set working directory to the script's location
# --------------------------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Working directory:", getwd(), "\n")

######################################
### Creating bias layer for insects ###
######################################

# Reference: https://scottrinnan.wordpress.com/2015/08/31/how-to-construct-a-bias-file-with-r-for-use-in-maxent-modeling/

# Step 1: Load libraries
library(data.table)  # For fast data handling
library(terra)      # For raster data manipulation
library(sf)          # For spatial data manipulation
library(dplyr)       # For data wrangling


# Step 2: Load location data and filter missing values
# Read from the zipped file directly
# Define the file path
file_path <- "0037565-241126133413365/occurrence.txt"
# Read the file using fread with selected columns
locations <- fread(file_path, select = c("decimalLongitude", "decimalLatitude"))
locations_clean <- locations[!is.na(locations$decimalLongitude) & !is.na(locations$decimalLatitude), ]
nrow(locations_clean)

# Convert the cleaned data frame to an sf object
locations_sf <- st_as_sf(locations_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Step 3: Load the raster data (e.g., CHELSA annual precipitation) and convert it to a mask
climdat <- (rast('CHELSA_annPrecip.tif')*0)+1 
climdat <- aggregate(climdat, fact = 10, fun = "mean", expand = TRUE)

# Step 4: Ensure CRS alignment
locations_sf <- st_transform(locations_sf, crs = crs(climdat, proj = TRUE))

# Step 5: Rasterize the observation points
# Count points per grid cell and create a bias layer
bias_raster <- rasterize(locations_sf, climdat, field = 1, fun = 'sum', background = 0)
bias_raster <- bias_raster*climdat
bias_raster_log <- log(bias_raster+1)
plot(bias_raster)
plot(bias_raster_log)

#sense check
nrow(locations_clean)
global(bias_raster, fun = "sum", na.rm = TRUE)

# Save the bias raster (optional)
writeRaster(bias_raster, "10km_bias_layer.tif", overwrite = TRUE)
writeRaster(bias_raster_log, "10km_bias_layer_log.tif", overwrite = TRUE)

#smooth
library(terra)

# Step 1: Load the log-transformed raster
bias_raster_log <- rast("1km_bias_layer_log.tif")

# Step 2: Convert raster to points (include zeros, exclude NA values)
points_log <- as.data.frame(bias_raster_log, xy = TRUE, na.rm = TRUE)  # Exclude NA values only
names(points_log) <- c("x", "y", "z")  # Rename columns for consistency

# Step 3: Create a matrix of points with coordinates and values
xyz <- as.matrix(points_log)  # Columns: x, y, z

# Step 4: Perform IDW interpolation
# Radius defines the search area for points
# Power controls the weighting (default is 2)
smoothed_raster_log <- interpIDW(
  x = bias_raster_log,  # Target raster (provides resolution and extent)
  y = xyz,              # Input points as a matrix
  radius = 5,           # Neighborhood search radius in raster units
  power = 2,            # Weighting power (default is 2)
  smooth = 1,           # Smoothing parameter
  maxPoints = 10,       # Maximum number of points to use for interpolation
  minPoints = 1,        # Minimum number of points for interpolation
  fill = NA             # Value to fill empty cells (e.g., NA for oceans)
)

# Step 5: Mask the smoothed raster to preserve original NA areas
final_smoothed_raster_log <- mask(smoothed_raster_log, bias_raster_log)

# Step 6: Visualize the smoothed raster
plot(final_smoothed_raster_log, main = "Smoothed Log Bias Layer (Using interpIDW)")

# Step 7: Save the smoothed raster
writeRaster(final_smoothed_raster_log, "1km_bias_layer_log_smoothed_idw.tif", overwrite = TRUE)
