# Data
Data for weeks 1-12 of the HPS are available at https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/. These may be downloaded individually, or (on a Mac/Unix system) by running the following shell commands in the data/ directory.

```
 for i in {1..12}
 do if [ $i -lt 10 ]
 then wget https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk$i/HPS_Week0$i'_PUF_CSV.zip'
 else wget https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk$i/HPS_Week$i'_PUF_CSV.zip'
 fi
 done
 unzip HPS_Week*PUF_CSV.zip
```

The file data/generate_empirical_pop.R then takes these .csv files as input and produces a single data frame as output.

For the empirical simulation, the script data/generate_empirical_pop.R outputs the population data file.

# Code and Instructions for use
1. Run the file code/generate_samples.R to produce informative sub-samples for use in the empirical simulation.
2. The code/fit_*.R files run the empirical simulation for the Gaussian and binary data, respectively.
  * fit_gbulm.R returns posterior predictions
  * the remaining models return the MCMC draws and must be processed using binary_posterior_preds.R and gtulm_posterior_preds.R
3. The files process_empirical_sim_results_gauss.R and process_empirical_sim_results_binary.R then take the posterior predictions and generate the associated tables and figures presented in Section 3 of the paper.
4. The file binary_data_analysis.R runs the data analysis on the full HPS phase 1 sample and process_data_analysis_results.R can then be run to get posterior predictions and to produce the tables and figures found in Section 4.
5. The simulation and data analyisis files source helper_functions.R, which contains helper functions for taking an informative sample and calculating direct estimates.
