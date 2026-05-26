library(terra)
library(sf)

# =========================================================
# User-defined paths
# =========================================================

# The paths below should be adapted by the user before running the script.

base_dir <- "path/to/input_data"
out_dir  <- "path/to/output_results"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

chm_2014_low_path  <- file.path(base_dir, "chm_2014_low_resolution.tif")
chm_2022_low_path  <- file.path(base_dir, "chm_2022_low_resolution.tif")
chm_2024_high_path <- file.path(base_dir, "chm_2024_high_resolution.tif")

crowns_2022_low_path  <- file.path(base_dir, "crowns_low_resolution.shp")
crowns_2024_high_path <- file.path(base_dir, "crowns_high_resolution.shp")

# =========================================================
# Load rasters
# =========================================================

chm_2014_low  <- terra::rast(chm_2014_low_path)
chm_2022_low  <- terra::rast(chm_2022_low_path)
chm_2024_high <- terra::rast(chm_2024_high_path)


# =========================================================
# Load crown polygons
# =========================================================

crowns_2022_low  <- sf::st_read(crowns_2022_low_path)
crowns_2024_high <- sf::st_read(crowns_2024_high_path)


# =========================================================
# Assign CRS
# =========================================================

target_crs <- "EPSG:3006"

terra::crs(chm_2014_low)  <- target_crs
terra::crs(chm_2022_low)  <- target_crs
terra::crs(chm_2024_high) <- target_crs

sf::st_crs(crowns_2022_low)  <- target_crs
sf::st_crs(crowns_2024_high) <- target_crs


# =========================================================
# Convert crown polygons to terra vectors
# =========================================================

crowns_2022_low_v  <- terra::vect(crowns_2022_low)
crowns_2024_high_v <- terra::vect(crowns_2024_high)

# =========================================================
# Define summary functions
# =========================================================

# 95th percentile
p95 <- function(x, ...) {
  quantile(x, probs = 0.95, na.rm = TRUE, names = FALSE)
}

# Mean value
mean_fun <- function(x, ...) {
  mean(x, na.rm = TRUE)
}

# =========================================================
# Extract heights for low-resolution crowns (2014-2022)
# =========================================================

h2014_low <- terra::extract(
  chm_2014_low,
  crowns_2022_low_v,
  fun = p95,
  na.rm = TRUE
)

h2022_low <- terra::extract(
  chm_2022_low,
  crowns_2022_low_v,
  fun = p95,
  na.rm = TRUE
)

h2014_low_mean <- terra::extract(
  chm_2014_low,
  crowns_2022_low_v,
  fun = mean_fun,
  na.rm = TRUE
)

h2022_low_mean <- terra::extract(
  chm_2022_low,
  crowns_2022_low_v,
  fun = mean_fun,
  na.rm = TRUE
)

crowns_2022_low$h2014 <- h2014_low[, 2]
crowns_2022_low$h2022 <- h2022_low[, 2]
crowns_2022_low$h2014_mean <- h2014_low_mean[, 2]
crowns_2022_low$h2022_mean <- h2022_low_mean[, 2]

crowns_2022_low$g_tot <- crowns_2022_low$h2022 - crowns_2022_low$h2014

# =========================================================
# Extract heights for high-resolution crowns (2014-2024)
# =========================================================

h2014_high <- terra::extract(
  chm_2014_low,
  crowns_2024_high_v,
  fun = p95,
  na.rm = TRUE
)

h2024_high <- terra::extract(
  chm_2024_high,
  crowns_2024_high_v,
  fun = p95,
  na.rm = TRUE
)

h2014_high_mean <- terra::extract(
  chm_2014_low,
  crowns_2024_high_v,
  fun = mean_fun,
  na.rm = TRUE
)

h2024_high_mean <- terra::extract(
  chm_2024_high,
  crowns_2024_high_v,
  fun = mean_fun,
  na.rm = TRUE
)

crowns_2024_high$h2014 <- h2014_high[, 2] / 10
crowns_2024_high$h2024 <- h2024_high[, 2]
crowns_2024_high$h2014_mean <- h2014_high_mean[, 2] / 10
crowns_2024_high$h2024_mean <- h2024_high_mean[, 2]

crowns_2024_high$g_tot <- crowns_2024_high$h2024 - crowns_2024_high$h2014

# =========================================================
# Convert low-resolution CHM values from decimeters to meters
# =========================================================

crowns_2022_low$h2014 <- crowns_2022_low$h2014 / 10
crowns_2022_low$h2022 <- crowns_2022_low$h2022 / 10
crowns_2022_low$h2014_mean <- crowns_2022_low$h2014_mean / 10
crowns_2022_low$h2022_mean <- crowns_2022_low$h2022_mean / 10

crowns_2022_low$g_tot <- crowns_2022_low$g_tot / 10

# =========================================================
# Remove unrealistic growth values
# =========================================================

# Remove negative growth values
crowns_2022_low$g_tot[crowns_2022_low$g_tot < 0] <- NA
crowns_2024_high$g_tot[crowns_2024_high$g_tot < 0] <- NA

# Remove extremely high growth values 
crowns_2022_low$g_tot[crowns_2022_low$g_tot > 6] <- NA
crowns_2024_high$g_tot[crowns_2024_high$g_tot > 6] <- NA

# =========================================================
# Calculate annual growth
# =========================================================

crowns_2022_low$g_ann <- crowns_2022_low$g_tot / 8
crowns_2024_high$g_ann <- crowns_2024_high$g_tot / 10

# =========================================================
# Check summary statistics
# =========================================================

summary(crowns_2022_low$g_tot)
summary(crowns_2022_low$g_ann)

summary(crowns_2024_high$g_tot)
summary(crowns_2024_high$g_ann)

# =========================================================
# Export results
# =========================================================

sf::st_write(
  crowns_2022_low,
  file.path(out_dir, "crowns_low_growth95.shp"),
  delete_layer = TRUE
)

sf::st_write(
  crowns_2024_high,
  file.path(out_dir, "crowns_high_growth95.shp"),
  delete_layer = TRUE
)