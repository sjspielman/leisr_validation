### LEISR validation

This repository accompanies the manuscript `Relative evolutionary rate inference in HyPhy with LEISR` by SJ Spielman\* and SLK Pond.

+ The Python script `simulate_compare_methods.py` performs all simulations and inferences.
	+ Dependencies: `pyvolve`, `phyphy`
+ The R script `process_and_plot.R` processes rate inferences and makes figures that appear in MS
	+ Dependencies: `tidyverse`, `cowplot`, `purrr`, `broom` 
+ `simulations_rates` directory contains all simulated data and rate inferences. Note that trees used for simulation were generated with the `ape` R package and are saved in the files `rtree<25,50,100.tre>`.

+ `test1e5.fna` and `test1e5.fna.site-rates.json` are a simulated (nucleotide) alignment of 10,000 sequences, with 100 sites, which is confirmed to run in `LEISR` (as indicated in MS). 

\*Please contact `stephanie.spielman<at>temple.edu` with questions.
