library(terra)
library(sf)
library(raster)
library(ForestTools)

# =========================================================
# User-defined paths
# =========================================================

# The paths below should be adapted by the user before running the script.

base_dir <- "path/to/input_data"
out_dir  <- "path/to/output_results"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

chm_path      <- file.path(base_dir, "chm_low_segmentation_input.tif")
treetops_path <- file.path(base_dir, "treetops_low_resolution.shp")

min_h <- 3  # Same height threshold as used during preprocessing

# =========================================================
# Load CHM
# =========================================================

chm_low_pre <- terra::rast(chm_path)

# Convert to RasterLayer for mcws and assign CRS
chm_low_r <- raster::raster(chm_low_pre)
raster::crs(chm_low_r) <- terra::crs(chm_low_pre)

# =========================================================
# Load treetops
# =========================================================

treetops_low <- sf::st_read(treetops_path, quiet = TRUE)

# Ensure EPSG:3006 is assigned
sf::st_crs(treetops_low) <- 3006

if (!terra::same.crs(chm_low_pre, terra::vect(treetops_low))) {
  stop("CHM and treetops do not have the same CRS")
}

# Convert treetops to SpatialPoints for mcws
ttops_low_sp <- as(treetops_low, "Spatial")

# =========================================================
# Run mcws for low-resolution CHM
# =========================================================

cat("Running mcws for low-resolution CHM...\n")
crowns_low <- mcws(
  treetops  = ttops_low_sp,
  CHM       = chm_low_r,
  minHeight = min_h,
  format    = "polygons"
)

# =========================================================
# Add ID and export crown polygons
# =========================================================

crowns_low$tree_id <- seq_len(nrow(crowns_low))
crowns_low_sf <- sf::st_as_sf(crowns_low)
sf::st_crs(crowns_low_sf) <- 3006
crowns_low_sf$area_m2 <- as.numeric(sf::st_area(crowns_low_sf))

sf::st_write(
  crowns_low_sf,
  file.path(out_dir, "crowns_low_resolution.shp"),
  delete_layer = TRUE
)

cat("Done! Number of crown segments, low resolution:", nrow(crowns_low_sf), "\n")