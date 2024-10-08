This repository contains the data and R code used to perform analyses for "Soil resource acquisition strategy modulates global plant nutrient and water economics" (Cheaib et al., 2024).
The manuscript is currently submitted to New Phytologist. Folder contents are as follows:

data = An initial file data/data_clean_beta.csv compiles all species nutrient acquisition strategies data, climate, isotopes, and soil data. 
Relevant calculations for the manuscript are calculated in scripts/Beta_calculus.R, which gives a .csv file (data/data_C3.csv) that is then used for all analyses. 
A metadata file for data/data_clean_beta.csv is also included in this folder

scripts = contains R scripts used for carbon cost calculations (scripts/Beta_calculus.R, and scripts/calc_optimal_vcmax.R), data analysis (scripts/stat_analysis_Beta_Final), and supplemental plots (scripts/Whittaker_diagram.R) 
functions = contains functions called for carbon cost calculations
