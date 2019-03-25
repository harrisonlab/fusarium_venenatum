#!/usr/bin/python
#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import os
import argparse
import re
import math
import numpy as np
from sets import Set
from collections import defaultdict
from collections import Counter
from operator import itemgetter
from Bio import SeqIO

ap = argparse.ArgumentParser()

ap.add_argument('--sig', required=True,type=str,help='Gregs Deseq2 output file showing LFC and P values in tsv format')
conf = ap.parse_args()

with open(conf.sig) as f:
    DEG_lines = f.readlines()

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

problem_genes = []
for line in DEG_lines[1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    if len(split_line) < 16:
        continue
    # print gene_id
    LFCs = itemgetter(2,4,6,8,10,12,14)(split_line)
    # print LFCs
    Pvals = itemgetter(3,5,7,9,11,13,15)(split_line)
    # print Pvals
    if any([p == "" for p in Pvals]):
        problem_genes.append(gene_id)
        continue
    if (
        any([float(p) < 0.01 for p in Pvals]) and
        any([abs(float(x)) > 2 for x in LFCs])
        ):
        print line
        # x=1
# print(",".join(problem_genes))
