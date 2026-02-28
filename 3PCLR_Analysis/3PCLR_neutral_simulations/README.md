# Neutral simulation of 3PCLR scores

This script uses a demographic model and runs 3PCLR with fixed drift rates
- main pipeline in 
`run_all_3pclr_parallel.sh`

## helper scripts

- demographic model (params from relate demographic inference): 
`demography_3pclr_var.py`

- script to invoke 3PCLR for each replicate: 
`run_3pclr_var.sh`

- script for maf filter:
`filterout_maf0.01.sh`


## Plotting and evaluation

- script to evaluate the 5\%, 1\% and 0.5\% cutoff and QQ plots:
`3PCLR_SimComparison.R`
