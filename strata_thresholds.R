library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# =========================================================
# User-defined paths
# =========================================================

# This script can be used for both low- and high-resolution datasets.
# Input files, field names and year references should be adapted
# depending on the dataset being processed.

base_dir <- "path/to/input_data"
out_dir  <- "path/to/output_results"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

crowns_path    <- file.path(base_dir, "crowns_high_growth95.shp")
elevation_path <- file.path(base_dir, "DEM.tif")
moist_path     <- file.path(base_dir, "soil_moisture.tif")

thresholds_path <- file.path(out_dir, "thresholds_high.csv")
cuts_path       <- file.path(out_dir, "strata_cuts_high.csv")
strata_path     <- file.path(out_dir, "strata_groups_high.gpkg")

# =========================================================
# Field names and processing parameters
# =========================================================

height_field   <- "h2024" # change year depending on dataset
growth_field   <- "g_ann"
dia_field      <- "cr_dia"
flatness_field <- "flatness"
area_field     <- "area_m2"

area_flag_limit <- 130
tpi_radius_m    <- 40

# =========================================================
# Load data
# =========================================================
crowns    <- st_read(crowns_path, quiet = TRUE)
elevation <- rast(elevation_path)
moist     <- rast(moist_path)

target_crs <- "EPSG:3006"

crs(elevation) <- target_crs
crs(moist)     <- target_crs
st_crs(crowns) <- target_crs

# =========================================================
# Polygon ID
# =========================================================
crowns$poly_id <- 1:nrow(crowns)
crowns_v <- vect(crowns)

# =========================================================
# Mean elevation
# =========================================================
zone_rast <- terra::rasterize(crowns_v, elevation, field = "poly_id")
zonal_vals <- terra::zonal(elevation, zone_rast, fun = "mean", na.rm = TRUE)

zonal_df <- as.data.frame(zonal_vals)
names(zonal_df) <- c("poly_id", "E_mean")

crowns <- crowns %>%
  left_join(zonal_df, by = "poly_id")

# =========================================================
# Soil moisture
# =========================================================
zone_rast_moist <- terra::rasterize(crowns_v, moist, field = "poly_id")

zone_vals  <- terra::values(zone_rast_moist, mat = FALSE)
moist_vals <- terra::values(moist, mat = FALSE)

moist_df <- data.frame(
  poly_id = zone_vals,
  mark_code = moist_vals
) %>%
  filter(!is.na(poly_id), !is.na(mark_code)) %>%
  group_by(poly_id) %>%
  summarise(
    mark_code = as.numeric(names(which.max(table(mark_code)))),
    .groups = "drop"
  )

crowns <- crowns %>%
  left_join(moist_df, by = "poly_id")

crowns <- crowns %>%
  mutate(
    mark_klass = case_when(
      mark_code == 1 ~ "dry-mesic",
      mark_code == 2 ~ "mesic-moist",
      mark_code == 3 ~ "moist-wet",
      mark_code == 4 ~ NA_character_,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(mark_klass))

# =========================================================
# Elevation classes
# =========================================================
cuts <- quantile(crowns$E_mean, probs = c(1/3, 2/3), na.rm = TRUE)

crowns <- crowns %>%
  mutate(
    e_band = case_when(
      E_mean < cuts[1]  ~ "EL",
      E_mean <= cuts[2] ~ "EM",
      E_mean > cuts[2]  ~ "EH",
      TRUE ~ NA_character_
    )
  )

# =========================================================
# Topographic Position Index (TPI)
# =========================================================
w_tpi <- focalMat(elevation, d = tpi_radius_m, type = "circle")

elev_mean_local <- focal(
  elevation,
  w = w_tpi,
  fun = mean,
  na.rm = TRUE,
  fillvalue = NA
)

tpi_rast <- elevation - elev_mean_local

zone_rast_tpi <- terra::rasterize(crowns_v, tpi_rast, field = "poly_id")
zonal_tpi <- terra::zonal(tpi_rast, zone_rast_tpi, fun = "mean", na.rm = TRUE)

zonal_tpi_df <- as.data.frame(zonal_tpi)
names(zonal_tpi_df) <- c("poly_id", "TPI_40")

crowns <- crowns %>%
  left_join(zonal_tpi_df, by = "poly_id")

tpi_cuts <- quantile(crowns$TPI_40, probs = c(1/3, 2/3), na.rm = TRUE)

crowns <- crowns %>%
  mutate(
    tpi_klass = case_when(
      TPI_40 < tpi_cuts[1]  ~ "TL",
      TPI_40 <= tpi_cuts[2] ~ "TM",
      TPI_40 > tpi_cuts[2]  ~ "TH",
      TRUE ~ NA_character_
    )
  )

# =========================================================
# Strata groups
# =========================================================
crowns <- crowns %>%
  mutate(
    group = if_else(
      !is.na(e_band) & !is.na(mark_klass) & !is.na(tpi_klass),
      paste(e_band, mark_klass, tpi_klass, sep = "_"),
      NA_character_
    )
  ) %>%
  filter(!is.na(group))

# =========================================================
# Crown diameter and flatness (change year depending on dataset)
# =========================================================
crowns <- crowns %>%
  mutate(
    cr_dia = 2 * sqrt(.data[[area_field]] / pi),
    flatness = if_else(
      is.finite(h2024) & h2024 > 0 & is.finite(h2024_mean),
      (h2024 - h2024_mean) / h2024,
      NA_real_
    ),
    flatness_valid = if_else(
      is.finite(flatness) & flatness > 0,
      1L, 0L
    ),
    underseg_flag = if_else(.data[[area_field]] > area_flag_limit, 1L, 0L)
  )

crowns_stats <- crowns

# =========================================================
# Export strata layer
# =========================================================

st_write(
  crowns, 
  strata_path, 
  delete_dsn = TRUE
)

# =========================================================
# Calculate thresholds per group
# =========================================================

thresholds <- crowns_stats %>%
  st_drop_geometry() %>%
  group_by(group) %>%
  summarise(
    h_thr   = quantile(.data[[height_field]], 0.75, na.rm = TRUE),
    d_thr   = quantile(.data[[dia_field]], 0.75, na.rm = TRUE),
    g_thr   = quantile(.data[[growth_field]], 0.15, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

write_csv(thresholds, thresholds_path)

# =========================================================
# Export strata class boundaries
# =========================================================

cuts_df <- tibble(
  e_low_mid    = as.numeric(cuts[1]),
  e_mid_high   = as.numeric(cuts[2]),
  tpi_low_mid  = as.numeric(tpi_cuts[1]),
  tpi_mid_high = as.numeric(tpi_cuts[2])
) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

write_csv(cuts_df, cuts_path)

cat("\nThresholds and strata cuts exported:\n")
cat(thresholds_path, "\n")
cat(cuts_path, "\n")