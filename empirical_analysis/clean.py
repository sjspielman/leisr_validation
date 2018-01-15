"""
    Clean alignment to remove gaps from first sequence, so that leisr and r4s produce rates for same sites.
"""

from Bio import AlignIO

aln = AlignIO.read("HRH1_aligned.fasta", "fasta")
refseq = str(aln[0].seq)
gaps = []
for x in range(len(refseq)):
    if refseq[x] == "-":
        gaps.append(x)
    else:
        print refseq[x]

with open("HRH1_aligned_cleaned.fasta", "w") as f:

    for record in aln:
        seq = str(record.seq)
        cleanseq = ""
        for x in range(len(seq)):
            if x not in gaps:
                cleanseq += seq[x]
            else:
                print x
                continue
        assert(len(cleanseq) == (len(seq) - len(gaps))), "Gap sites not removed"
        f.write(">" + str(record.id) + "\n" + cleanseq + "\n")