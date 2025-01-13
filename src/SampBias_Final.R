# --------------------------------------------------------------------------------------
# Set working directory to the script's location
# --------------------------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Working directory:", getwd(), "\n")

# --------------------------------------------------------------------------------------
# Load necessary libraries
# --------------------------------------------------------------------------------------
library(sampbias)
library(data.table)
library(dismo)
library(raster)
library(blockCV)
library(terra)
library(ecospat)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(GeoThinneR)
library(prettymapr)
library(leaflet)
library(tmap)
library(ggspatial)
library(ggplot2)

# --------------------------------------------------------------------------------------
# Load WWF ecoregions shapefile
# --------------------------------------------------------------------------------------
cat("Loading WWF ecoregions shapefile...\n")
wwf_eco <- vect("./official/wwf_terr_ecos.shp")

# --------------------------------------------------------------------------------------
# Load and filter hornet occurrence data
# --------------------------------------------------------------------------------------
cat("Reading and filtering hornet occurrence data...\n")
data <- fread("0036526-241126133413365.csv")
hornet_data <- data[, c(13, 22, 23), with = FALSE] %>%
  na.omit()  # Remove rows with NA values

# Inspect data
head(hornet_data)
unique(hornet_data$scientificName)

# Convert to spatial format
hornet_sf <- st_as_sf(hornet_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Plot hornet data
fig <- ggplot() +
  annotation_map_tile(type = "osm") +
  geom_sf(data = hornet_sf, aes(color = "red"), size = 0.5) +
  theme_minimal() +
  labs(title = "Hornet Locations", color = "Hornet Data") +
  theme(legend.position = "none")
#fig

# --------------------------------------------------------------------------------------
# Filter Vespa velutina data
# --------------------------------------------------------------------------------------
cat("Filtering Vespa velutina data...\n")
velutina_data <- hornet_data[grep("^Vespa velutina", hornet_data$scientificName), ]
unique(velutina_data$scientificName)

# Convert to spatial format
velutina_data_sf <- st_as_sf(velutina_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Plot Vespa velutina data
fig_velutina <- ggplot() +
  annotation_map_tile(type = "osm") +
  geom_sf(data = velutina_data_sf, aes(color = "red"), size = 0.5) +
  theme_minimal() +
  labs(title = "Vespa velutina Locations", color = "Hornet Data") +
  theme(legend.position = "none")

# Save plot
ggsave(filename = "vespa_velutina_map.pdf", plot = fig_velutina, width = 8, height = 6, units = "in")

# --------------------------------------------------------------------------------------
# Clean Vespa velutina data and filter by range
# --------------------------------------------------------------------------------------
cat("Cleaning and filtering Vespa velutina data...\n")

# Ensure valid geometries for WWF ecoregions
wwf_eco_sf <- st_make_valid(st_as_sf(wwf_eco))

# Add longitude and latitude columns
velutina_data_sf <- velutina_data_sf %>%
  mutate(
    decimalLongitude = st_coordinates(.)[, 1],
    decimalLatitude = st_coordinates(.)[, 2]
  )

# Filter out points outside defined range (e.g., Americas)
velutina_points_filtered <- velutina_data_sf %>%
  filter(decimalLongitude > -30)

# Simplify WWF ecoregions for faster processing
wwf_eco_sf_simplified <- st_simplify(wwf_eco_sf, dTolerance = 0.01)

# Ensure CRS compatibility
if (!st_crs(velutina_points_filtered) == st_crs(wwf_eco_sf_simplified)) {
  velutina_points_filtered <- st_transform(velutina_points_filtered, st_crs(wwf_eco_sf_simplified))
}

# Filter out points outside ecoregions
valid_intersections <- st_intersects(velutina_points_filtered, wwf_eco_sf_simplified)
keep_indices <- lengths(valid_intersections) > 0
velutina_points_cleaned <- velutina_points_filtered[keep_indices, ]

# --------------------------------------------------------------------------------------
# Plot cleaned Vespa velutina data
# --------------------------------------------------------------------------------------
cat("Plotting cleaned Vespa velutina data...\n")

fig_cleaned <- ggplot() +
  geom_sf(data = wwf_eco_sf_simplified, fill = "lightgrey", color = "darkgrey", size = 0.2) +
  geom_sf(data = velutina_points_cleaned, color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Cleaned Vespa velutina Occurrence Data") +
  theme(legend.position = "none")

# Save plot to file
ggsave("cleaned_hornet_occurrence_map.png", plot = fig_cleaned, width = 10, height = 7, dpi = 300)

# --------------------------------------------------------------------------------------
# Intersect cleaned Vespa velutina points with WWF ecoregions
# --------------------------------------------------------------------------------------
cat("Intersecting Vespa velutina points with WWF ecoregions...\n")

# Perform spatial intersection
intersection <- st_intersects(velutina_points_cleaned, wwf_eco_sf_simplified)

# Add ECO_NAME to each point based on the intersection
velutina_points_cleaned <- velutina_points_cleaned %>%
  mutate(ECO_NAME = sapply(intersection, function(idx) {
    if (length(idx) > 0) {
      wwf_eco_sf_simplified$ECO_NAME[idx[1]]  # Take the first match
    } else {
      NA
    }
  }))

# Filter WWF ecoregions to retain only those with at least one occurrence
ecoregions_with_occurrences <- velutina_points_cleaned %>%
  st_drop_geometry() %>%
  filter(!is.na(ECO_NAME)) %>%
  group_by(ECO_NAME) %>%
  summarise(count = n()) %>%
  filter(count > 0) %>%
  pull(ECO_NAME)

# Subset WWF ecoregions based on valid ECO_NAMEs
filtered_wwf_eco <- wwf_eco_sf_simplified %>%
  filter(ECO_NAME %in% ecoregions_with_occurrences)

# --------------------------------------------------------------------------------------
# Plot the filtered ecoregions with Vespa velutina occurrences (no legend for ECO_NAME)
# --------------------------------------------------------------------------------------
cat("Plotting filtered ecoregions with Vespa velutina occurrences (no legend)...\n")

fig_filtered <- ggplot() +
  geom_sf(data = filtered_wwf_eco, fill = "lightblue", color = "darkblue", size = 0.2) +
  geom_sf(data = velutina_points_cleaned, color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Filtered Ecoregions with Vespa velutina Occurrences") +
  theme(legend.position = "none")  # Remove the legend entirely

# Save the plot
ggsave("filtered_ecoregions_with_occurrences_no_legend.png", plot = fig_filtered, width = 10, height = 7, dpi = 300)

cat("Filtered ecoregions and plot without legend created successfully.\n")

# --------------------------------------------------------------------------------------
# Split filtered WWF ecoregions into European and Asian shapefiles
# --------------------------------------------------------------------------------------
cat("Splitting filtered WWF ecoregions into European and Asian regions...\n")

# Extract centroids for each polygon to determine their locations
centroids <- st_centroid(filtered_wwf_eco)

# Extract longitude and latitude from centroids
centroids_coords <- st_coordinates(centroids)

# Add longitude and latitude to the data
filtered_wwf_eco <- filtered_wwf_eco %>%
  mutate(
    longitude = centroids_coords[, 1],
    latitude = centroids_coords[, 2]
  )

# Define spatial criteria for Europe and Asia
europe_eco <- filtered_wwf_eco %>%
  filter(longitude <= 60)

asia_eco <- filtered_wwf_eco %>%
  filter(longitude > 60)

# Save the European and Asian shapefiles
st_write(europe_eco, "filtered_wwf_eco_europe.shp", append = FALSE)
st_write(asia_eco, "filtered_wwf_eco_asia.shp", append = FALSE)

cat("European and Asian shapefiles created successfully.\n")

# --------------------------------------------------------------------------------------
# Simplify and buffer European and Asian shapefiles correctly
# --------------------------------------------------------------------------------------
cat("Simplifying and buffering European and Asian shapefiles with proper projection...\n")

# Load the European and Asian shapefiles
europe_eco <- st_read("filtered_wwf_eco_europe.shp")
asia_eco <- st_read("filtered_wwf_eco_asia.shp")

# Define a suitable projected CRS (e.g., Robinson projection for global use)
projected_crs <- "+proj=robin +datum=WGS84"

# Reproject the shapefiles to the projected CRS
europe_projected <- st_transform(europe_eco, crs = projected_crs)
asia_projected <- st_transform(asia_eco, crs = projected_crs)

# Dissolve polygons into single unified polygons
europe_unified <- st_union(europe_projected) %>%
  st_make_valid()  # Ensure geometry is valid

asia_unified <- st_union(asia_projected) %>%
  st_make_valid()

# Apply a buffer (1-degree equivalent in projected CRS, approximately 111 km)
buffer_distance <- 111000  # 1 degree in meters
europe_buffered <- st_buffer(europe_unified, dist = buffer_distance)
asia_buffered <- st_buffer(asia_unified, dist = buffer_distance)

# Reproject back to the original geographic CRS (EPSG:4326)
europe_final <- st_transform(europe_buffered, crs = 4326)
asia_final <- st_transform(asia_buffered, crs = 4326)

# Simplify geometries for performance (optional)
europe_final_simplified <- st_simplify(europe_final, dTolerance = 0.01)
asia_final_simplified <- st_simplify(asia_final, dTolerance = 0.01)

# Save the unified and buffered shapefiles
st_write(europe_final_simplified, "europe_buffered_unified_fixed.shp", append = FALSE)
st_write(asia_final_simplified, "asia_buffered_unified_fixed.shp", append = FALSE)

cat("Unified and buffered shapefiles created successfully with correct projection.\n")

# --------------------------------------------------------------------------------------
# Filter ALL HORNET DATA to retain only those overlapping with asia_final_simplified
# --------------------------------------------------------------------------------------
cat("Filtering ALL HORNET DATA points for overlaps with Asia region...\n")

hornet_data_sf <- st_as_sf(
  hornet_data,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326  # Assuming WGS84 CRS
)

# Ensure geometries are valid for intersection
hornet_data_sf <- st_make_valid(hornet_data_sf)
asia_final_simplified <- st_make_valid(asia_final_simplified)

# Perform the spatial intersection to retain only overlapping points
asia_overlap <- st_intersects(hornet_data_sf, asia_final_simplified, sparse = FALSE)

# Filter the points that overlap with Asia
hornet_asia_points <- hornet_data_sf[asia_overlap, ]

# Extract coordinates from the filtered points
coords <- st_coordinates(hornet_asia_points)

# Convert the result to a data.frame with the specified columns
hornet_asia_df <- hornet_asia_points %>%
  st_drop_geometry() %>%
  mutate(
    species = "velutina_asia",  # Add the species column with constant value
    decimalLongitude = coords[, 1],  # Extract longitude
    decimalLatitude = coords[, 2]   # Extract latitude
  ) %>%
  select(species, decimalLongitude, decimalLatitude)  # Retain required columns only

# Inspect the resulting dataset
head(hornet_asia_df)

# --------------------------------------------------------------------------------------
# Run sampling bias analysis FOR ASIA
# --------------------------------------------------------------------------------------
cat("Running sampling bias analysis for Asia...\n")

samp_mask <- calculate_bias(x = hornet_asia_df, res = 0.25, terrestrial = TRUE)
summary(samp_mask)
#plot(samp_mask)

proj <- project_bias(samp_mask)
#map_bias(proj, type = "log_sampling_rate")

all_bias <- proj[[4]]
all_bias_log10 <- log10(all_bias + 1)  # Convert to log scale
plot(all_bias_log10)

# Save raster to file
output_filename <- file.path(getwd(), "all_bias_log10_025degree.tif")
writeRaster(all_bias_log10, filename = output_filename, overwrite = TRUE)

# --------------------------------------------------------------------------------------
# Filter ALL HORNET DATA to retain only those overlapping with europe_final_simplified
# --------------------------------------------------------------------------------------
cat("Filtering ALL HORNET DATA points for overlaps with Europe region...\n")

hornet_data_sf <- st_as_sf(
  hornet_data,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326  # Assuming WGS84 CRS
)

# Filter dataset to retain only rows where 'scientificName' contains 'crabro'
crabro_data_sf_filtered <- hornet_data_sf[grepl("crabro", hornet_data_sf$scientificName, ignore.case = TRUE), ]

# Inspect unique scientific names
unique(crabro_data_sf_filtered$scientificName)

# Ensure geometries are valid for intersection
hornet_data_sf <- st_make_valid(crabro_data_sf_filtered)
europe_final_simplified <- st_make_valid(europe_final_simplified)

# Perform the spatial intersection to retain only overlapping points
europe_overlap <- st_intersects(crabro_data_sf_filtered, europe_final_simplified, sparse = FALSE)

# Filter the points that overlap with Europe
crabro_europe_points <- crabro_data_sf_filtered[europe_overlap, ]

# Extract coordinates from the filtered points
coords <- st_coordinates(crabro_europe_points)

# Convert the result to a data.frame with the specified columns
crabro_europe_df <- crabro_europe_points %>%
  st_drop_geometry() %>%
  mutate(
    species = "crabro_europe",  # Add the species column with constant value
    decimalLongitude = coords[, 1],  # Extract longitude
    decimalLatitude = coords[, 2]   # Extract latitude
  ) %>%
  select(species, decimalLongitude, decimalLatitude)  # Retain required columns only

# Inspect the resulting dataset
head(crabro_europe_df)

# --------------------------------------------------------------------------------------
# Create map for European Crabro distribution
# --------------------------------------------------------------------------------------
cat("Creating map for European Crabro distribution...\n")

# Load Europe map
europe_map <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")

# Create the map with points from crabro_europe_df
pdf("Crabro_Europe_Distribution.pdf", width = 10, height = 10)  # Open PDF device
ggplot() +
  geom_sf(data = europe_map, fill = "lightgray", color = "black", size = 0.2) +  # Europe map
  geom_point(data = crabro_europe_df, 
             aes(x = decimalLongitude, y = decimalLatitude), 
             color = "red", size = 2) +  # Points
  coord_sf(xlim = c(-10, 40), ylim = c(35, 70), expand = FALSE) +  # Adjust map extent
  labs(title = "Distribution of Crabro in Europe",
       x = "Longitude", y = "Latitude") +
  theme_minimal()
dev.off()  # Close PDF device

# --------------------------------------------------------------------------------------
# Run sampling bias analysis FOR EUROPE
# --------------------------------------------------------------------------------------
cat("Running sampling bias analysis for Europe...\n")
samp_mask_eu <- calculate_bias(x = crabro_europe_df, res = 0.25,terrestrial = TRUE)
summary(samp_mask_eu)
#plot(samp_mask_eu)

proj_eu <- project_bias(samp_mask_eu)
#map_bias(proj_eu, type = "log_sampling_rate")

all_bias_eu <- proj_eu[[4]]
all_bias_log10_eu <- log10(all_bias_eu + 1)  # Convert to log scale
plot(all_bias_log10_eu)

# Save raster to file
output_filename <- file.path(getwd(), "all_bias_log10_eu_025degree.tif")
writeRaster(all_bias_log10_eu, filename = output_filename, overwrite = TRUE)

# --------------------------------------------------------------------------------------
# Combine Asia and Europe bias rasters
# --------------------------------------------------------------------------------------
cat("Combining Asia and Europe bias rasters...\n")

# Load the two .tif files
europe <- rast("all_bias_log10_eu_025degree.tif")
asia <- rast("all_bias_log10_025degree.tif")

# Define a global raster grid with the same resolution
global_template <- rast("CHELSA_annPrecip.tif")
global_template[!is.na(global_template)] <- 1
global_template_025 <- aggregate(global_template, fact = c(30, 30), fun = "any")
writeRaster(global_template_025, "global_template_025_mask.tif", overwrite = TRUE)

# Resample both rasters to align them to the global grid
europe_aligned <- resample(europe, global_template, method = "near")
asia_aligned <- resample(asia, global_template, method = "near")

# Combine the rasters using merge()
combined_raster <- merge(europe_aligned, asia_aligned)

# Check the combined raster
plot(combined_raster)
# Save the combined raster
writeRaster(combined_raster, "combined_global_bias.tif", overwrite = TRUE)
# Write the hornet_data to a text file
write.table(hornet_data, 
            file = "hornet_data.txt", 
            sep = "\t",            # Use tab as the separator
            row.names = FALSE,     # Do not include row numbers
            quote = FALSE)         # Do not quote character strings



