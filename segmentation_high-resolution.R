library(terra)
library(sf)
library(raster)
library(ForestTools)
library(dplyr)

# =========================================================
# User-defined paths
# =========================================================

# The paths below should be adapted by the user before running the script.

base_dir <- "path/to/input_data"
out_dir  <- "path/to/output_results"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

chm_path      <- file.path(base_dir, "chm_high_segmentation_input.tif")
treetops_path <- file.path(base_dir, "treetops_high_resolution.shp")

# =========================================================
# Processing parameters
# =========================================================

min_h    <- 3
tile_size <- 1000
overlap_buf <- 10

# =========================================================
# Load CHM and treetops
# =========================================================

chm_high_pre  <- terra::rast(chm_path)
treetops_high <- sf::st_read(treetops_path)

# Ensure EPSG:3006 is assigned
sf::st_crs(treetops_high) <- 3006

# =========================================================
# Create tile grid covering the CHM extent
# =========================================================

e <- terra::ext(chm_high_pre)   # xmin, xmax, ymin, ymax

x_breaks <- seq(e[1], e[2], by = tile_size)
y_breaks <- seq(e[3], e[4], by = tile_size)

# Ensure the final tile fully covers the CHM extent
if (tail(x_breaks, 1) < e[2]) x_breaks <- c(x_breaks, e[2])
if (tail(y_breaks, 1) < e[4]) y_breaks <- c(y_breaks, e[4])

crowns_list <- list()
tile_id <- 1L

# =========================================================
# Process tiles
# =========================================================

for (ix in seq_len(length(x_breaks) - 1)) {
  for (iy in seq_len(length(y_breaks) - 1)) {
    
    xmin <- x_breaks[ix]
    xmax <- x_breaks[ix + 1]
    ymin <- y_breaks[iy]
    ymax <- y_breaks[iy + 1]
    
    # Skip tile if any extent value is missing
    if (any(is.na(c(xmin, xmax, ymin, ymax)))) next
    
    tile_core_ext <- terra::ext(xmin, xmax, ymin, ymax)
    
    # Create extended tile extent with overlap buffer
    tile_ext <- terra::ext(
      xmin - overlap_buf,
      xmax + overlap_buf,
      ymin - overlap_buf,
      ymax + overlap_buf
    )
    
    # Crop CHM to tile extent
    chm_tile <- try(terra::crop(chm_high_pre, tile_ext), silent = TRUE)
    if (inherits(chm_tile, "try-error") || is.null(chm_tile) || terra::ncell(chm_tile) == 0) {
      next
    }
    
    # Create polygon for extended tile extent
    
    tile_poly_v  <- terra::as.polygons(tile_ext)
    terra::crs(tile_poly_v) <- terra::crs(chm_high_pre)
    
    tile_poly_sf <- sf::st_as_sf(tile_poly_v)
    tile_poly_sf <- sf::st_transform(tile_poly_sf, sf::st_crs(treetops_high))
    
    # Create polygon for core tile extent
    
    tile_core_poly_v  <- terra::as.polygons(tile_core_ext)
    terra::crs(tile_core_poly_v) <- terra::crs(chm_high_pre)
    
    tile_core_sf <- sf::st_as_sf(tile_core_poly_v)
    tile_core_sf <- sf::st_transform(tile_core_sf, sf::st_crs(treetops_high))
    
    
    # Select treetops located within the extended tile
    
    idx_mat <- sf::st_intersects(treetops_high, tile_poly_sf, sparse = FALSE)
    if (nrow(idx_mat) == 0) next
    idx <- idx_mat[, 1]
    
    tt_tile <- treetops_high[idx, ]
    if (nrow(tt_tile) == 0) {
      next
    }
    
    # Convert objects for mcws
    chm_tile_r <- raster::raster(chm_tile)
    raster::crs(chm_tile_r) <- terra::crs(chm_high_pre)
    
    tt_tile_sp <- as(tt_tile, "Spatial")
    
    cat("Running tile", tile_id,
        "x:", xmin, "-", xmax,
        "y:", ymin, "-", ymax,
        "treetops (full/overlap):", nrow(tt_tile), "\n")
    
    # =========================================================
    # Run mcws for current tile
    # =========================================================
    crowns_tile <- try(
      mcws(
        treetops  = tt_tile_sp,
        CHM       = chm_tile_r,
        minHeight = min_h,
        format    = "polygons"
      ),
      silent = TRUE
    )
    
    if (!inherits(crowns_tile, "try-error") &&
        !is.null(crowns_tile) &&
        nrow(crowns_tile) > 0) {
      
      crowns_tile_sf <- sf::st_as_sf(crowns_tile)
      sf::st_crs(crowns_tile_sf) <- 3006
      
      # Calculate centroids for crown polygons
      centroids <- sf::st_centroid(crowns_tile_sf)
      
      # Retain crowns whose centroids fall within the core tile
      inside_mat <- sf::st_intersects(centroids, tile_core_sf, sparse = FALSE)
      keep_idx <- if (nrow(inside_mat) > 0) inside_mat[, 1] else rep(FALSE, nrow(crowns_tile_sf))
      crowns_tile_sf <- crowns_tile_sf[keep_idx, ]
      
      if (nrow(crowns_tile_sf) > 0) {
        
        crowns_tile_sf <- sf::st_join(
          crowns_tile_sf,
          tt_tile[, c("tree_id")],
          join = sf::st_contains,
          left = TRUE
        )
        
        # Remove duplicate crown geometries
        crowns_tile_sf <- crowns_tile_sf[!duplicated(sf::st_as_text(sf::st_geometry(crowns_tile_sf))), ]
        
        crowns_tile_sf$tile_id <- tile_id
        crowns_list[[length(crowns_list) + 1L]] <- crowns_tile_sf
      }
    }
    
    # Remove temporary objects between tiles
    rm(chm_tile, chm_tile_r, tt_tile, tt_tile_sp,
       crowns_tile, crowns_tile_sf, centroids,
       tile_poly_v, tile_poly_sf, tile_core_poly_v, tile_core_sf)
    gc()
    
    tile_id <- tile_id + 1L
    
  }
}

# =========================================================
# Merge tiles and remove duplicate crowns
# =========================================================

if (length(crowns_list) == 0) {
  stop("No crown polygons were created. Check input data and parameters")
}

crowns_high_all <- do.call(rbind, crowns_list)

# Calculate crown area
crowns_high_all$area_m2 <- as.numeric(sf::st_area(crowns_high_all))

# Retain one crown polygon per tree_id
crowns_high_all <- crowns_high_all |>
  dplyr::group_by(tree_id) |>
  dplyr::slice_max(order_by = area_m2, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

cat("After deduplication - number of crowns:", nrow(crowns_high_all), "\n")
cat("Duplicate tree_id values:", sum(duplicated(crowns_high_all$tree_id)), "\n")

# =========================================================
# Create summary table per tile
# =========================================================

tile_counts <- crowns_high_all |>
  sf::st_drop_geometry() |>
  dplyr::count(tile_id, name = "crowns")

tile_coords <- expand.grid(
  ix = seq_len(length(x_breaks) - 1),
  iy = seq_len(length(y_breaks) - 1)
)

tile_coords <- tile_coords |>
  dplyr::mutate(
    tile_id = dplyr::row_number(),
    xmin = x_breaks[ix],
    xmax = x_breaks[ix + 1],
    ymin = y_breaks[iy],
    ymax = y_breaks[iy + 1]
  ) |>
  dplyr::select(tile_id, xmin, xmax, ymin, ymax)

tile_table <- tile_counts |>
  dplyr::left_join(tile_coords, by = "tile_id") |>
  dplyr::select(tile_id, xmin, xmax, ymin, ymax, crowns) |>
  dplyr::arrange(tile_id)

write.csv(
  tile_table,
  file.path(out_dir, "crowns_per_tile.csv"),
  row.names = FALSE
)

# =========================================================
# Export crown polygons
# =========================================================

sf::st_write(
  crowns_high_all,
  file.path(out_dir, "crowns_high_resolution.shp"),
  delete_layer = TRUE
)

cat("Done! Number of crown polygons (high resolution, tiled processing):", nrow(crowns_high_all), "\n")