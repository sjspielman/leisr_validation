### LEISR validation

This repository accompanies the manuscript `Relative evolutionary rate inference in HyPhy with LEISR` by SJ Spielman\* and SLK Pond, currently in press at [PeerJ](https://peerj.com/).

Please see the "Publication Release" under the releases tab for the version of this repository that accompanies the published manuscript.

#### Dependencies

+ `requirements.txt` may be used to install all Python dependencies with the code `pip install -r requirements.txt` (with a prepended `sudo` as needed), called from this directory.

+ The script `install_R_requirements.R` can be run with `Rscript install_R_requirements.R` to install all R dependencies.

+ [`HyPhy`](http://hyphy.org) version >= 2.3.8 and [`Rate4Site`](https://www.tau.ac.il/~itaymay/cp/rate4site.html) must be installed and available in your system path.


#### Repository contents 

+ `simulation_analysis` contains code and data to reproduce simulations and subsequent inference, visualization
	+ The Python script `simulate_compare_methods.py` performs all simulations and inferences.
		+ Dependencies: `pyvolve`, `phyphy`
	+ The R script `process_and_plot.R` processes rate inferences and makes figures that appear in MS
		+ Dependencies: `tidyverse`, `cowplot`, `broom` 
	+ `simulations_rates` directory contains all simulated data and rate inferences. Note that trees used for simulation were generated with the `ape` R package and are saved in the files `rtree<25,50,100.tre>`.

+ `empirical_analysis` contains code and data to reproduce inference and visualization on empirical dataset. See `README.md` within the directory for details.

+ `test1e5.fna` and `test1e5.fna.site-rates.json` are a simulated (nucleotide) alignment of 10,000 sequences, with 100 sites, which is confirmed to run in `LEISR` (as indicated in MS). 



\*Please contact `stephanie.spielman<at>temple.edu` with questions.
