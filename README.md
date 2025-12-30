# ACCRSC_rogersEtAl2025 ReadMe
This respository contains the code and paths to structured data to reproduce all analyses and figures shown in "Differential modulation of aversive signaling by expectation across the cingulate cortex," by Sophie A. Rogers, Corinna S. Oswell, Lindsay L. Ejoh, and Gregory Corder (University of Pennsylvania). All scripts were written by S.A.R.
## Table of contents:
1. Data structures
   * shockStabilityACC.mat
   * shockStabilityRSC.mat
   * tfcACC.mat
   * tfcRSC.mat
2. Code
   * shockStabilityAnalysisPipeline.m
   * tfcAnalysisPipeline.m
   * special functions
## Data structures
All structures are organized with pre-processed (longitudinally registered, detrended for photobleaching, aligned to trial start times, etc.) data as follows:
 - shockExp
      * calcium
          * matrix of stimulus-aligned trial and session concatenated fluorescence traces for every ROI detected across animals. sampling rate = 20Hz
      * freezing (for TFC only)
          * matrix of stimulus-aligned trial and session concatenated freezing scores from ezTrack (Pennington et al., 2019). sampling rate = 15Hz
      * nCells
          * array of number of cells per animal in order as added to the calcium matrix. the first entry is zero to permit later operations, then animals 1-n.
      * tensors
          * a cell containing a tensor of samples x ROIs x trials size for each animal containing calcium fluorescence over each trial.
### All data are in this [Google Drive folder](https://drive.google.com/drive/folders/1WfiPiSTn8dxITaw50XrYa7ysv36WHvn7?usp=sharing)
#### shockStability.mat
This structure contains the data for the ACC-implanted animals in Figs. 1 and 2.
#### shockStabilityRSC.mat
This structure contains the data for the RSC-implanted animals in Figs. 1 and 2.
#### rogers2025data2.mat
This structure contains the data for the ACC-implanted animals in Figs. 3-5.
#### rogers2025CalciumData.mat
This structure contains the data for the RSC-implanted animals in Figs. 3-5.
## Code
### mainAnalysisPipeline.m
1. Load data. Choose which to load
2. Extract data for each animal
3a. Single cell dynamics over trials for shock stability experiments
    * Amplitudes, peak times, rise times, decays
    * Identifies pre-shock, shock, and post-shock neurons
    * Plots:
        * population average activity over trials
        * pooled fraction of cell types
        * mean ± SEM fraction cell types over animals
        * cdf plots of dynamic statistics
        * average activity of each cell type in each session
        * heatmaps of peaktimes in trial 1 vs. trial 8
        * overlaps across trials
3b. Same as 4a for the TFC experiments
4. Calculate and plot neural synchrony
5. Calculate and plot pooled PCA, temporal drift
    * Plots:
        * distance between trials of each session
        * distance between vs. within trials
        * distance as a function of temporal lag
        * trajectory speed / acceleration
6. Plot speed/acceleration over trials
7. Calculate and plot trajectories for individual animals
8. For TFC, analyze & visualize behavior







