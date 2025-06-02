
README - Cactus Embryo Cross-sectional Intensity Analysis (Updated)
=================================================================

Overview
--------
This codebase includes MATLAB scripts and supporting functions for analyzing cross-sectional fluorescence intensities from long-term imaging of Cactus in Drosophila embryos.

Key functions:
1. Reads microscopy data preprocessed via FIJI/ImageJ using Bio-Formats.
2. Segments nuclei and cytoplasmic regions across timepoints.
3. Aligns intensity time series from multiple embryos.
4. Averages profiles and visualizes data with standard deviation bands.
5. Outputs figures and save the data.


Key Scripts
-----------
- script_cactlt_analysis.m           : Main pipeline script for full dataset processing.
- run_analyze_cactlt.m               : Analyzes nuclear and cytoplasmic channels for all embryos.
- average_intensity_figure_embryos.m : Aligns and averages signals from multiple embryos.
- analyze_cactlt_nuc.m / analyze_cactlt_cyto.m : Analyze individual channels.
- script_gathercactltfiles.m         : Aggregates processed MAT files into one structure.
- Functions/                         : Contains utility functions used across the pipeline:
    * ftn_gradfit6.m, ftn_gradfit7.m, ftn_gradfit8_plot_slope.m
    * borderFinder.m, circfit.m, domainMeas.m, fit_gaussian2.m
    * calcbg, epsDU,

Requirements
------------
- MATLAB R2020a or newer
- Toolboxes:
  * Image Processing Toolbox
  * Signal Processing Toolbox
- Bio-Formats for ImageJ/FIJI
  * Version used: Bio-Formats 6.12.0

Usage
-----
1. Open MATLAB and navigate to the root directory of this repository.
2. Run "script_gathercactltfiles" to compile a list of filenames
3. Run `script_cactlt_analysis.m` to begin processing.
4. To generate summary plots from aligned data, run `average_intensity_figure_embryos.m`.

Output
------
Figures showing averaged nuclear and cytoplasmic intensity with error shading are saved in both `.fig` and `.jpg` formats.

