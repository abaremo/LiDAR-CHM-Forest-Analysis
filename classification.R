library(sf)
library(dplyr)
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

crowns_path     <- file.path(base_dir, "strata_groups_high.gpkg")
thresholds_path <- file.path(base_dir, "thresholds_high.csv")
candidate_path  <- file.path(out_dir, "old_candidates_high.shp")

# =========================================================
# Fixed biological thresholds
# =========================================================

height_floor   <- 12.0
dia_floor      <- 4.5
growth_cap     <- 0.06
flatness_floor <- 0.08

# =========================================================
# Load data
# =========================================================

crowns <- st_read(crowns_path, quiet = TRUE)
thresholds <- read_csv(thresholds_path, show_col_types = FALSE)

# =========================================================
# Classify old-tree candidates
# =========================================================

crowns <- crowns %>%
  left_join(thresholds, by = "group") %>%
  mutate(
    pass_hard_filters = if_else(
      flatness_valid == 1 &
        h2024 >= height_floor &   # change year depending on dataset
        cr_dia >= dia_floor &
        g_ann <= growth_cap &
        flatness >= flatness_floor,
      1L, 0L
    ),
    
    crit_height = if_else(
      pass_hard_filters == 1 & h2024 >= h_thr,   # change year depending on dataset
      1L, 0L
    ),
    
    crit_dia = if_else(
      pass_hard_filters == 1 & cr_dia >= d_thr,
      1L, 0L
    ),
    
    crit_growth = if_else(
      pass_hard_filters == 1 & g_ann <= g_thr,
      1L, 0L
    ),
    
    score_old = crit_height + crit_dia + crit_growth,
    
    class_old = case_when(
      pass_hard_filters == 0 ~ "Non_candidate",
      score_old == 3 ~ "A_full",
      score_old == 2 ~ "B_partial",
      TRUE ~ "Non_candidate"
    ),
    
    cand_old = if_else(
      class_old %in% c("A_full", "B_partial"),
      1L, 0L
    )
  )

# =========================================================
# Export A and B candidates
# =========================================================

old_candidates <- crowns %>%
  filter(class_old %in% c("A_full", "B_partial")) %>%
  select(
    poly_id,
    group,
    h2024, # change year depending on dataset
    g_ann,
    cr_dia,
    flatness,
    flatness_valid,
    underseg_flag,
    crit_height,
    crit_dia,
    crit_growth,
    score_old,
    class_old
  ) %>%
  rename(
    flat_val = flatness_valid,
    underseg = underseg_flag,
    c_h = crit_height,
    c_d = crit_dia,
    c_g = crit_growth
  )

print(table(old_candidates$class_old, useNA = "ifany"))
print(table(old_candidates$score_old, useNA = "ifany"))

st_write(old_candidates, candidate_path, delete_layer = TRUE)

cat("\nScript completed.\n")
cat("Candidate layer exported to:\n", candidate_path, "\n")


table(crowns$class_old, crowns$underseg_flag)