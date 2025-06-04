README – Image Analysis and FRAP Analysis Codebase

Overview
This codebase provides MATLAB scripts and functions for the analysis of fluorescence recovery after photobleaching (FRAP) experiments. It supports image preprocessing, nuclear segmentation and tracking, intensity extraction, kinetic parameter estimation, and reconstruction of missing nuclei.

Workflow Summary:
1. Load `.czi` microscopy image files and convert them to `.mat` format.
2. Segment nuclei and extract nuclear and cytoplasmic intensity profiles.
3. Fit recovery curves to estimate nuclear import/export rates.
4. Generate plots and save results for downstream analysis.

Requirements
- MATLAB R2020a or later  
- Bio-Formats plugin (for reading `.czi` files) -  This code was run using the 2017 version of Bio-Formats plugin.
- Image Processing Toolbox  
- Optimization Toolbox  

File Structure

Main Scripts (Run in the order listed)

- openFilesWithSamePrefix_savedata.m  
  Reads `.czi` files (possibly split into parts) and saves the image data as `.mat` files.

- script_gathermatfiles.m  
  Scans target directories to gather `.mat` files and exports metadata into an Excel file `Bleaching_stats.xlsx`.

- script_mats_analysis.m  
  Batch-processing script that analyzes multiple `.mat` files. Calls `run_analyze_bleach` for each dataset.

- run_analyze_bleach.m  
  Wrapper script that initiates analysis of bleaching experiments.

- analyze_bleach_czi.m  
  Loads `.mat` files, extracts image stacks and timestamps, calls `find_nucleizoom_new_comp_set2`, and computes intensity profiles for bleached and neighboring nuclei.

- find_nucleizoom_new_comp_set2.m  
  Core function for segmenting, labeling, and tracking nuclei and cytoplasm across frames.

- reconstruct_missing_nucleus.m  
  Reconstructs missing nuclear masks in selected frames using motion interpolation.

- dl_vs_t.m  
  Extracts nuclear and cytoplasmic intensity time series using masks over image data.

- fit_bleach.m  
  Fits the fluorescence recovery curve to estimate:
  - `k_in`, `k_out`: nuclear import/export rates  
  - `c0`: initial nuclear intensity  
  - `r2`: coefficient of determination  
  Optionally generates plots comparing model fit with raw data.

Utility Functions

- fraccircshift.m  
  Applies fractional circular shifts to an image matrix (used for subpixel motion compensation).

- traceobject.m  
  Extracts object boundaries from region properties.

- persistent_closest_neighbors.m  
  Identifies nuclei consistently located near the bleached nucleus across all frames.

- plot_neighbors.m  
  Plots neighboring nuclei around the bleached nucleus for visualization.

- openczi.m  
  Opens `.czi` files using Bio-Formats.

- ftn_FRAP.m  
  FRAP model fitting function.

- track_nuclei_labels_exact.m  
  Assigns consistent labels to nuclei across frames by registering centroids. Includes a lookback mechanism (default = 3).

Usage Instructions

1. Set the base path for image data in `openFilesWithSamePrefix_savedata.m` and run it to generate `.mat` files.
2. Run `script_gathermatfiles.m` to process the files and compile metadata.
3. Use `script_mats_analysis.m` to perform segmentation, intensity extraction, fitting, and visualization.
4. To generate and save plots, set `yesplot = true` in relevant scripts.
5. Processed results are stored in a structure "Soln_1cpt" and saved in a mat file "Mat/1componentfits". Exported info is also stored in `czi mat results` directory, inside folders like `free gfp`, `2 gfp`, or `normal gfp`, based on the experimental type.
6. If you wish to fit the 1x1x embryos to a 2-component model, run "script_fit_twocompFRAP_cact"
7. To reproduce the plots from the FRAP figure in the paper, run "script_postanaly_CactLT"
