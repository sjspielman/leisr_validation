from phyphy import *

l = LEISR(type = "protein", tree = "HRH1.tre", alignment = "HRH1_aligned_cleaned.fasta", model = "WAG")
l.run_analysis()
x = Extractor(l)
x.extract_csv("leisr.csv")