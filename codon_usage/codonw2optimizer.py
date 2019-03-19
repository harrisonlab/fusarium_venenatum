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

ap.add_argument('--codonw',required=True,type=str,help='codonw .cutot file')


conf = ap.parse_args()

with open(conf.codonw) as f:
    weight_lines = f.readlines()


#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------

aa_list = [
    'Phe', 'Sre', 'Tyr', 'Cys', 'Leu', 'TER', 'Trp',
    'Pro', 'His', 'Arg', 'Gln', 'Ile', 'Thr', 'Asn',
    'Ser', 'Lys', 'Arg', 'Met', 'Val', 'Ala', 'Asp', 'Gly',
    'Glu']
element_list = []
for line in weight_lines[0:-2]:
    line = line.rstrip()
    split_line = line.split(" ")
    for element in split_line:
        if element:
            if element.isdigit():
                # print float(element)
                if float(element) < 10:
                    # print element
                    element_list.append(element)
            elif not any(x in element for x in aa_list):
                # print element
                if re.match("^[a-zA-Z]+.*", element):
                    element = element[0:3]
                    # print element
                element_list.append(element)

def pairwise(it):
    it = iter(it)
    while True:
        yield next(it), next(it)

codon_dict = {
'UUU':'Phe','UUC':'Phe',
'UUA':'Leu','UUG':'Leu','CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu',
'AUU':'Ile','AUC':'Ile','AUA':'Ile',
'AUG':'Met',
'GUU':'Val', 'GUC':'Val','GUA':'Val','GUG':'Val',
'UCU':'Ser','UCC':'Ser','UCA':'Ser','UCG':'Ser',
'CCU':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
'ACU':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
'GCU':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
'UAU':'Tyr','UAC':'Tyr',
'UAA':'Ter','UAG':'Ter',
'CAU':'His','CAC':'His',
'CAA':'Gln','CAG':'Gln',
'AAU':'Asn','AAC':'Asn',
'AAA':'Lys','AAG':'Lys',
'GAU':'Asp','GAC':'Asp',
'GAA':'Glu','GAG':'Glu',
'UGU':'Cys','UGC':'Cys',
'UGA':'Ter',
'UGG':'Trp',
'CGU':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
'AGU':'Ser','AGC':'Ser',
'AGA':'Arg','AGG':'Arg',
'GGU':'Gly','GGC':'Gly', 'GGA':'Gly','GGG':'Gly'
}

out_list = []
for a, b in pairwise(element_list):
    # print a
    aa = codon_dict[a]
    # print("\t".join([aa,a,b]))
    out_list.append(": ".join([a,b]))

print("; ".join(out_list))



# for element in element_list:
#     print element
