"""
    Randomly select 100 alignments from Jack et al. 2016
"""

import os
from Bio import AlignIO
import random

get = 100
minseq = 25

original_trees = "raxml/"
original_alignments = "phy/"

final = "enzyme_data/"

alns = [x for x in os.listdir(original_alignments) if x.endswith("phy")]
random.shuffle(alns)


x = 0
for aln in alns:
    name = aln.split(".phy")[0]
    
    with open(original_alignments + aln, "r") as f:
        seqs = AlignIO.read(f, "phylip-relaxed")
        
    if len(seqs) >= minseq:
        x += 1
    
        treename = "RAxML_bestTree." + name + ".txt"
        with open(original_trees + treename, "r") as f:
            treestring = f.read()
        
        treestring = treestring.replace("lcl|","")
    
    
        output = name + ".fna"
        with open(final + output, "w") as f:
            for record in seqs:
                f.write(">" + str(record.id).replace("lcl|","") + "\n" + str(record.seq) + "\n")
            f.write("\n" + treestring)    
    if x == get:
        break
    



# 132L_A.fna:301
# 13PK_A.fna:301
# 1A4L_A.fna:301
# 1A50_B.fna:301
# 1A8H_A.fna:154
# 1A8S_A.fna:301
# 1AF7_A.fna:301
# 1AJ8_A.fna:301
# 1B73_A.fna:301
# 1B8F_A.fna:301
# 1B93_A.fna:301
# 1BF2_A.fna:301
# 1BJO_A.fna:301
# 1BQC_A.fna:75
# 1BZC_A.fna:83
# 1C0K_A.fna:301
# 1C2T_A.fna:301
# 1C82_A.fna:79
# 1CBG_A.fna:301
# 1CEL_A.fna:301
# 1CF2_O.fna:301
# 1CTN_A.fna:148
# 1CTT_A.fna:301
# 1CZ1_A.fna:301
# 1DBT_A.fna:301
# 1DEK_A.fna:35
# 1DJ0_A.fna:301
# 1DPG_A.fna:301
# 1E2A_A.fna:301
# 1E3V_A.fna:216
# 1E6E_A.fna:301
# 1E7Q_A.fna:301
# 1EG7_A.fna:301
# 1EUU_A.fna:203
# 1EXP_A.fna:158
# 1FF3_A.fna:301
# 1FPS_A.fna:301
# 1FR8_A.fna:121
# 1G6T_A.fna:301
# 1G72_A.fna:301
# 1G79_A.fna:301
# 1G8F_A.fna:301
# 1GA8_A.fna:131
# 1GCB_A.fna:301
# 1GDH_A.fna:301
# 1GP1_A.fna:301
# 1GPJ_A.fna:301
# 1GTP_A.fna:301
# 1H3I_A.fna:39
# 1HR6_B.fna:301
# 1HV9_A.fna:301
# 1IVH_A.fna:301
# 1J79_A.fna:301
# 1JOF_A.fna:88
# 1K30_A.fna:66
# 1L1D_A.fna:301
# 1L7N_A.fna:201
# 1LJL_A.fna:301
# 1MBB_A.fna:301
# 1MHY_D.fna:33
# 1MQW_A.fna:301
# 1NBA_A.fna:301
# 1NIR_A.fna:157
# 1NU3_A.fna:164
# 1OAC_A.fna:301
# 1OH9_A.fna:301
# 1OPM_A.fna:128
# 1OS7_A.fna:301
# 1PMA_B.fna:301
# 1POW_A.fna:301
# 1QAZ_A.fna:43
# 1QFM_A.fna:301
# 1QI9_A.fna:97
# 1QJE_A.fna:122
# 1QRG_A.fna:301
# 1REQ_A.fna:301
# 1S3I_A.fna:266
# 1S95_A.fna:77
# 1SML_A.fna:301
# 1SMN_A.fna:183
# 1STC_E.fna:301
# 1TAH_B.fna:301
# 1TLP_E.fna:34
# 1TRK_B.fna:301
# 1TYF_A.fna:301
# 1UAQ_A.fna:301
# 1UQT_A.fna:301
# 1W0H_A.fna:86
# 1WGI_A.fna:301
# 1YCF_A.fna:301
# 1YON_A.fna:301
# 1YVE_I.fna:145
# 1ZE1_A.fna:301
# 206L_A.fna:143
# 2DOR_A.fna:301
# 2JCW_A.fna:301
# 3CSM_A.fna:301
# 3ECA_A.fna:301
# 3PVA_A.fna:301
# 7ATJ_A.fna:301
    