"""
    Stephanie J. Spielman
    Code to simulate alignments and infer rates w/ LEISR and R4S.
    
    Note that this script requires packages pyvolve and phyphy (both available in PyPI).
"""

from phyphy import *
import os
import subprocess
import re
from pyvolve import *



"""
    R code to make simulation trees:
    library(ape)
    for (n in c(25, 50, 100)){
        t <- rtree(n)
        t <- compute.brlen(t, rgamma, 0.25, 0.25)
        write.tree(t, paste0("rtree", n, ".tre"))
    } 
"""

def rs4_to_csv(infile, outfile):
    out = "site,rate\n"
    with open(infile, "r") as f:
        lines = f.readlines()
    for line in lines:
        find = re.match(r"^[1-9]", line)
        if find is not None:
            line2 = re.split("\s+", line)
            out += line2[0] + "," + line2[2] + "\n"
    with open(outfile, "w") as f:
        f.write(out.strip())


output_directory = "simulation_rates/"

for rep in range(10):
    print "Replicate", rep
    m = Model("WAG", alpha = 0.4, num_categories = 20)
    

    for i in [25, 50, 100]:

        print "Running for", i, "taxa"
        
    
        seqfile = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + ".fasta"
        ratefile = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_rates.txt"
        infofile = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_info.txt"
        treefile = output_directory + "rtree" + str(i) + ".tre"
    
        hyphy_homo_csv   = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_homo_hyphy.csv"
        hyphy_gamma_csv   = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_gamma_hyphy.csv"
        hyphy_gdd_csv   = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_gdd_hyphy.csv"
        ml_homo_csv = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_ML_homo_r4s.csv"
        ml_gamma_csv = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_ML_gamma_r4s.csv"
        beb_gamma_csv = output_directory + "sim" + str(i) + "taxa_rep" + str(rep) + "_BEB_gamma_r4s.csv"
        
        
        print "     Simulating"
        t = read_tree(file = treefile)
        p = Partition(size=100, models = m)
        e = Evolver(tree = t, partitions = p)
        e(seqfile = seqfile, ratefile = ratefile, infofile = infofile)
        
        ### Hyphy inference #####
        print "     Running hyphy"
        print "          No rv"
        r = LEISR(alignment = seqfile, tree = treefile, type = "protein", output = "sim.json", model = "WAG")
        r.run_analysis()
        p = Extractor(r)
        p.extract_csv(hyphy_homo_csv)    
        os.system("rm sim.json")

        print "          Gamma"        
        r = LEISR(alignment = seqfile, tree = treefile, type = "protein", output = "sim.json", model = "WAG", rate_variation = "Gamma")
        r.run_analysis()
        p = Extractor(r)
        p.extract_csv(hyphy_gamma_csv)    
        os.system("rm sim.json")
      
        
        #### r4s inference, note I  also use 4 categories because inference fails with default 16, and bonus this achieves more comparable performance with hyphy anyways. #####
        print "     Running r4s"
        # ML with homogenous bl optimization
        subprocess.call("rate4site -s " + seqfile + " -t " + treefile + " -y out.txt -Mw -im -bh > /dev/null", shell=True)
        try:
            rs4_to_csv("out.txt", ml_homo_csv)
        except:
            pass
        os.system("rm -f out.txt")
        os.system("rm -f r4s.res")
        os.system("rm -f TheTree.txt")   
        
        # ML with gamma bl optimization
        subprocess.call("rate4site -s " + seqfile + " -t " + treefile + " -y out.txt -Mw -im -k 4 > /dev/null", shell=True)
        try:
            rs4_to_csv("out.txt", ml_gamma_csv)
        except:
            pass
        os.system("rm -f out.txt")
        os.system("rm -f r4s.res")
        os.system("rm -f TheTree.txt")
                

        
