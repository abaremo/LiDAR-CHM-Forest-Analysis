# Tree crown segmenation and classification using LiDAR-derived CHMs
This repository contains R scripts developed as part of a master´s thesis at the University of Gävle. The workflow uses airborne LiDAR-derived canopy height models (CHMs) for preprosessing, tree crown segmentation, growth estimation, and classification of Scots pine trees with structural characteristics assiciated with older trees.

## Workflow overview

1. Preprocessing of CHMs
2. Variable Window Filter (VWF) treetop detection
3. Crown segmentation using MCWS
4. Growth estimation from multi-temporal CHMs
5. Environmental stratification
6. Threshold-based classification of candidate old Scots pine trees

### Repository structure

- `preprocessing.R` – preprocessing of CHMs, Gaussian smoothing and treetop detection using VWF
- `segmentation_low-resolution.R` – MCWS crown segmentation for low-resolution CHM
- `segmentation_high-resolution.R` – tiled MCWS crown segmentation for high-resolution CHM
- `growth_calculation.R` – extraction of height and growth metrics from multi-temporal CHMs
- `strata_thresholds.R` – environmental stratification and calculation of group-specific thresholds
- `classification.R` – classification of candidate old Scots pine trees
## Input data

The workflow requires airborne LiDAR-derived canopy height models (CHMs), a digital elevation model (DEM), and soil moisture data. Input paths must be adapted by the user before running the scripts.

## Authors

Clara Spik and Alexandra Bäremo  
University of Gävle, Sweden
