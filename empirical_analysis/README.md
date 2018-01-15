+ Original files `HRH1_aligned.fasta` and `HRH1.tre` are taken from the paper [Measuring evolutionary rates of proteins in a structural context](https://f1000research.com/articles/6-1845/v1). 

+ The script `clean.py` remove columns which are gaps in the first sequence from the alignment to produce the file `HRH1_aligned_cleaned.fasta`, which is used for input to LEISR and Rate4Site

+ Use `python run_leisr.py` to infer LEISR rates (CSV output in `leisr.csv`) and `bash run_r4s.sh` to infer Rate4Site rates (produces `r4sOrig.res` which was manually converted to CSV `r4s.csv`.

+ Use `plot.R` to obtain scatterplot and *R^2* 