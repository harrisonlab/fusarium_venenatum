#!/usr/bin/python

'''

'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--fpkm',required=True,type=str,help='tab delimited file containing geneID, length, then fpkm by treatment')
ap.add_argument('--min',required=True,type=str,help='minimum fpkm threshold')
ap.add_argument('--max',required=False,type=str,help='maximum fpkm threshold')


conf = ap.parse_args()

with open(conf.fpkm) as f:
    fpkm_lines = f.readlines()

min = conf.min

max = conf.max

print fpkm_lines[0].rstrip()
for line in fpkm_lines[1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    if all(float(x) > float(min) for x in split_line[3:]):
        if max != '':
            if all(float(x) < float(max) for x in split_line[3:]):
                print line
        else:
            print line
    # else:
    #     print 'monkeys'
