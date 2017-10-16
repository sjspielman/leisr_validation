### LEISR validation

This repository accompanies the manuscript `Relative evolutionary rate inference in HyPhy with LEISR` by SJ Spielman\* and SLK Pond.

+ `simulate_compare_methods.py` performs all simulations and inferences.
	+ Dependencies: `pyvolve`, `phyphy`
	+ Note that the package [`phyphy`](https://github.com/sjspielman/phyphy) is currently *in development* and exists to faciliate HyPhy inference. You will have to edit HyPhy inference code to match your system setup, in the context of this package. The package is slated for official release in late 2017, early 2018, at which point code here will be updated to reflect the final usage of `phyphy`.
+ `process_and_plot.R` processes rate inferences and makes figures that appear in MS
	+ Dependencies: `tidyverse`, `cowplot`, `purrr`, `broom` 
+ `simulations_rates` directory contains all simulated data and rate inferences.

\*Please contact `stephanie.spielman<at>temple.edu` with questions.