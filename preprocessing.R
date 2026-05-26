library(terra)
library(sf)
library(raster)
library(ForestTools)
library(dplyr)
library(igraph)

# =========================================================
# User-defined paths
# =========================================================

# The paths below should be adapted by the user before running the script.
# All processing parameters are kept in the script to make the workflow transparent.

base_dir <- "path/to/input_data"
out_dir  <- "path/to/output_results"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

chm_low_path  <- file.path(base_dir, "chm_low_resolution.tif")
chm_high_path <- file.path(base_dir, "chm_high_resolution.tif")
aoi_path      <- file.path(base_dir, "study_area.shp")

# =========================================================
# Data import, CRS and initial checks
# =========================================================

# Load CHM rasters
chm_low  <- terra::rast(chm_low_path)
chm_high <- terra::rast(chm_high_path)


# Target CRS (SWEREF 99 TM)
target_crs <- "EPSG:3006"

terra::crs(chm_high) <- target_crs
terra::crs(chm_low)  <- target_crs

# Expected raster resolution
expected_res_high <- c(0.5, 0.5)  # meter
expected_res_low  <- c(1.0, 1.0)  # meter

# Verify raster resolution and CRS
stopifnot(all.equal(terra::res(chm_high), expected_res_high))
stopifnot(all.equal(terra::res(chm_low),  expected_res_low))

# =========================================================
# Data type and unit conversion
# =========================================================

# Remove categories / color table from low-resolution CHM
levels(chm_low) <- NULL

# Unit conversion: high-resolution CHM already in meters,
# low-resolution CHM converted from decimeters to meters
chm_high_m <- chm_high
chm_low_m  <- chm_low / 10

# Assign descriptive layer names
names(chm_high_m) <- "height_m"
names(chm_low_m)  <- "height_m"

# =========================================================
# AOI buffering and CHM clipping
# =========================================================

# Load study area (AOI)
aoi <- sf::st_read(aoi_path)

# Ensure AOI uses the same CRS as the CHMs
aoi <- sf::st_transform(aoi, target_crs)

# Create buffer around study area (meters)
buffer_dist <- 5
aoi_buf <- sf::st_buffer(aoi, dist = buffer_dist)

# Convert AOI to terra vector format
aoi_buf_v <- terra::vect(aoi_buf)

# Crop and mask high-resolution CHM to buffered AOI
chm_high_clip <- terra::crop(chm_high_m, aoi_buf_v)
chm_high_clip <- terra::mask(chm_high_clip, aoi_buf_v)

# Crop and mask low-resolution CHM to buffered AOI
chm_low_clip <- terra::crop(chm_low_m, aoi_buf_v)
chm_low_clip <- terra::mask(chm_low_clip, aoi_buf_v)

# Height threshold for vegetation filtering
height_thresh <- 3

# Create masks for pixels >= 3 m in original CHMs
orig_mask_high <- chm_high_clip >= height_thresh
orig_mask_low  <- chm_low_clip  >= height_thresh

# =========================================================
# Gaussian smoothing of CHMs before VWF
# =========================================================

# 3x3 Gaussian kernel
gauss3 <- matrix(c(
  1,2,1,
  2,4,2,
  1,2,1
), nrow = 3, byrow = TRUE)

# Normalize kernel weights to sum to 1
gauss3 <- gauss3 / sum(gauss3)

# Smooth high-resolution CHM before VWF
chm_high_smooth <- terra::focal(
  chm_high_clip,
  w = gauss3,
  fun = "sum",
  na.policy = "omit",
  pad = TRUE
)

# Smooth low-resolution CHM before VWF
chm_low_smooth <- terra::focal(
  chm_low_clip,
  w = gauss3,
  fun = "sum",
  na.policy = "omit",
  pad = TRUE
)

# Remove clipped CHMs from memory
rm(chm_high_clip, chm_low_clip); gc()

# =========================================================
# Vegetation masking after smoothing
# =========================================================

# CHMs for VWF treetop detection:
# smoothed CHMs thresholded to NA below 3 m
# Saved temporarily for internal processing

chm_high_smooth_thr <- terra::ifel(
  chm_high_smooth < height_thresh, NA, chm_high_smooth, filename = file.path(out_dir, "chm_high_VWF.tif"),
  overwrite = TRUE) 


chm_low_smooth_thr <- terra::ifel(
  chm_low_smooth < height_thresh, NA, chm_low_smooth, filename = file.path(out_dir, "chm_low_VWF.tif"),
  overwrite = TRUE
)

# CHMs for crown segmentation:
# pixels must be >= 3 m in both smoothed and original CHMs
# Saved as final segmentation inputs

seg_mask_high <- (chm_high_smooth >= height_thresh) & orig_mask_high
seg_mask_low  <- (chm_low_smooth  >= height_thresh) & orig_mask_low

chm_high_seg <- terra::ifel(seg_mask_high, chm_high_smooth, NA,
                            filename = file.path(out_dir, "chm_high_segmentation_input.tif"),
                            overwrite = TRUE)

chm_low_seg <- terra::ifel(seg_mask_low, chm_low_smooth, NA,
                           filename = file.path(out_dir, "chm_low_segmentation_input.tif"),
                           overwrite = TRUE)


# =========================================================
# Variable Window Filter (VWF) function
# =========================================================

# Function for combined GLOBAL + FINE VWF detection
# including treetop thinning

run_vwf_once <- function(chm_terra,
                         a_g, b_g,
                         a_f, b_f,
                         minHeight = 3,
                         thin_dist = 0.5) {
  
  # Convert terra raster to RasterLayer for ForestTools
  chm_r <- raster::raster(chm_terra)
  raster::crs(chm_r) <- terra::crs(chm_terra)
  
  # GLOBAL window detection
  tops_global <- ForestTools::vwf(
    CHM       = chm_r,
    winFun    = function(h_m) a_g * h_m + b_g,
    minHeight = minHeight
  )
  
  # FINE window detection
  tops_fine <- ForestTools::vwf(
    CHM       = chm_r,
    winFun    = function(h_m) a_f * h_m + b_f,
    minHeight = minHeight
  )
  
  # Return empty result if no treetops are detected
  if ((is.null(tops_global) || nrow(tops_global) == 0) &&
      (is.null(tops_fine)   || nrow(tops_fine)   == 0)) {
    return(list(
      treetops = NULL,
      stats    = data.frame(n_tops = 0)
    ))
  }
  
  # Merge detections and label source
  crs_target <- sf::st_crs(3006)
  
  tops_g <- if (!is.null(tops_global) && nrow(tops_global) > 0) {
    x <- sf::st_as_sf(tops_global)
    sf::st_crs(x) <- crs_target
    x$source <- "global"
    x
  } else {
    NULL
  }
  
  tops_f <- if (!is.null(tops_fine) && nrow(tops_fine) > 0) {
    x <- sf::st_as_sf(tops_fine)
    sf::st_crs(x) <- crs_target
    x$source <- "fine"
    x
  } else {
    NULL
  }
  
  if (!is.null(tops_g) && !is.null(tops_f)) {
    toppar_sf <- rbind(tops_g, tops_f)
  } else if (!is.null(tops_g)) {
    toppar_sf <- tops_g
  } else if (!is.null(tops_f)) {
    toppar_sf <- tops_f
  } else {
    return(list(
      treetops = NULL,
      stats    = data.frame(n_tops = 0)
    ))
  }
  
  # Extract CHM height at each detected treetop
  toppar_sf$h_chm <- terra::extract(chm_terra, terra::vect(toppar_sf))[, 2]
  
  # Assign priority: fine > global
  toppar_sf$prio <- ifelse(toppar_sf$source == "fine", 2L, 1L)
  
  # Thinning using buffer distance and connected components
  buf  <- sf::st_buffer(toppar_sf, thin_dist)
  grp  <- sf::st_intersects(buf)
  memb <- igraph::components(igraph::graph_from_adj_list(grp))$membership
  toppar_sf$grp <- memb
  
  toppar_thin <- toppar_sf |>
    dplyr::group_by(grp) |>
    dplyr::arrange(dplyr::desc(prio), dplyr::desc(h_chm)) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  
  # Remove temporary grouping variable
  toppar_thin$grp <- NULL
  
  n_tops <- nrow(toppar_thin)
  
  list(
    treetops = toppar_thin,
    stats    = data.frame(n_tops = n_tops)
  )
}
# =========================================================
# Final VWF parameters and processing
# =========================================================

# Final VWF parameters
a_g_final  <- 0.045
b_g_final  <- 0.90
a_f_final  <- 0.035
b_f_final  <- 0.80
thin_final <- 0.5
min_h      <- 3

# Run VWF on high-resolution CHM
res_high <- run_vwf_once(
  chm_terra  = chm_high_smooth_thr,
  a_g        = a_g_final, b_g = b_g_final,
  a_f        = a_f_final, b_f = b_f_final,
  minHeight  = min_h,
  thin_dist  = thin_final
)
treetops_high <- res_high$treetops
print(res_high$stats)

# Run VWF on low-resolution CHM
res_low <- run_vwf_once(
  chm_terra  = chm_low_smooth_thr,
  a_g        = a_g_final, b_g = b_g_final,
  a_f        = a_f_final, b_f = b_f_final,
  minHeight  = min_h,
  thin_dist  = thin_final
)
treetops_low <- res_low$treetops
print(res_low$stats)

# =========================================================
# Export treetops and summary statistics
# =========================================================

# Add unique tree IDs
if (!is.null(treetops_high)) {
  treetops_high$tree_id <- seq_len(nrow(treetops_high))
}

if (!is.null(treetops_low)) {
  treetops_low$tree_id  <- seq_len(nrow(treetops_low))
}

# Export treetops as shapefiles
if (!is.null(treetops_high)) {
  sf::st_write(
    treetops_high,
    file.path(out_dir, "treetops_high_resolution.shp"),
    delete_layer = TRUE
  )
}

if (!is.null(treetops_low)) {
  sf::st_write(
    treetops_low,
    file.path(out_dir, "treetops_low_resolution.shp"),
    delete_layer = TRUE
  )
}

# Combine summary statistics
stats_df <- rbind(
  data.frame(dataset = "high", res_high$stats),
  data.frame(dataset = "low",  res_low$stats)
)

# Export summary statistics as CSV
write.csv(
  stats_df,
  file.path(out_dir, "treetop_stats.csv"),
  row.names = FALSE
)