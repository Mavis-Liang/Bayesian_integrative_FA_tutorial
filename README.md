# Tutorial

Paper: TBD
Tutorial (in e-book): [https://mavis-liang.github.io/Bayesian_integrative_FA_tutorial/](https://mavis-liang.github.io/Bayesian_integrative_FA_tutorial_book/)

This repository contains the raw code for simulations and applications. There are also some results uploaded. Note that the nutrition data is not publicly available.

## Organization of this Repository

- Main_code
  - `installation.r`: how each package is installed.
  - `sim_scenarioXX.r`: the code to run each model to each senarios, Tetris exlucded. Exucution of post-processing (with `post_xxxx.r` in the `function` folder) and calculation of RV/FN also included. Scenario 1 is based on MOM-SS (corresponding to the Scenario 2 in the manuscript); Scenario 2 is based on BMSFA; Scenario 3 is based on PFA (corresponding to the Scenario 1 in the manuscript); Scenario 4 mimics nutrition data based on Tetris (corresponding to the Scenario 4 in the manuscript); Scenario 5 mimics gene expression data based on Tetris (corresponding to the Scenario 5 in the manuscript); Scenario 6 is based on SUFA (corresponding to the Scenario 3 in the manuscript).
  - `sim_Tetris.r`: runs Tetris. Because Tetris takes far more time to run, we run it alone every time.
  - `plot_box.r` and `tab_NumFactor.R`: the code that generates the box plots and table for numbers of factors to display the simulation results shown in the manuscript.
  - `run_nutrition.r` and `run_curatedData.r`: fit the each models to the real nutrition and gene-expression data.
  - `nutrition_analysis.R` and `CuratedOvarian_analysis.R`: pre-processing and post-processing (determining numbers of factors, calculating MSE, visualizations) for the real-data applications.

- Functions
  - `gen_scenarioXX.r`: data generations for each simulation scenarios. The functions are used in `sim_scenarioXX.r` in the `main_code` folder, so that it can be executed easily when I run 50 repeatitions (generated data are randomly different across repeatitions).
  - `post_xxxx.r`: post-process (such as varimax, OP, or calculating the mean) the output (usually MCMC chains) from the respect model to obtain point estimates (loadings and covariances) that can be compared across methods. Used in `sim_scenarioXX.r`.
  - `measurements.r`: calculating RV, FN of different estimated quantities with the ground true. The function is executed in `sim_scenarioXX.r`.
  - `RV_pd.R`: Calculating a second definition of RV - it defines specifically for positive semi-definite matrices (matrices that can be represented as XX^T). We don't use this definition in our manuscript because our loadings some times do not satisfy the positive semi-definite criteria and the estimated loadings and true can have different numbers of columns.

- RDS
  Results are stored as RDS files. To read the files in R, use `readRDS("file_path.filename.rds"). Within each RDS, we can see a list containing the true generated data (if this is the results of simulations); estimated loadings and covariances of each method; the metrics (RV and FN).

- Figs
  heatmaps, networks plots, etc.
