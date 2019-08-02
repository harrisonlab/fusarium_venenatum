#!/usr/bin/python

'''
Generation of a scaffolded illumina assembly from satsuma alignment
of an illumina assembly to a reference genome.
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

ap.add_argument('--satsuma',required=True,type=str,help='satsuma alignment file')
# ap.add_argument('--illumina',required=True,type=str,help='A fasta file of the assembled contigs')
# ap.add_argument('--reference',required=True,type=str,help='fasta file of the reference genome')
# ap.add_argument('--gff',required=True,type=str,help='gff file of illumina gene models to be transferred onto the new genome.')

conf = ap.parse_args()

with open(conf.satsuma) as f:
    satsuma_lines = f.readlines()



#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------

class Scaffold_obj(object):
    """
    A set of contigs aligned to a single reference chromosome.
    Attributes:
        ref_name: A string representing the chromsome name.
        contig_obj_dict: diction of contig object associated with this chormsomes
        scaffold_gapped_seq: Scaffolded sequence with 100 Ns at contig boundaries.
    """

    def __init__(self):
        """Return a Annot_obj whose name is *gene_name*"""
        self.ref_name = ''
        self.contig_obj_dict = defaultdict(list)
        self.scaffolded_gapped_seq = ''

    # def set_conditions(self, gff_elements):
    #     """Reset conditions_tested to a list of conditions"""
    #     self.contig = gff_elements[0]
    #     self.start = gff_elements[3]
    #     self.stop = gff_elements[4]
    #     self.strand = gff_elements[6]
    #     gene_features = gff_elements[8].split(';')
    #     gene_id = gene_features[0]
    #     gene_id = gene_id.replace('ID=', '')
    #     self.gene_name = gene_id

class Contig_obj(object):
    """
    Information on query contigs and the locations that they align to.
    Attributes:
        contig_name: A string representing the chromsome name.
        hit_contigs: set of contig names
        hits_dict: a dictionary of hit_obj associated with this query contig
        scaffold_gapped_seq: Scaffolded sequence with 100 Ns at contig boundaries.
    """

    def __init__(self):
        """Return a Annot_obj whose name is *gene_name*"""
        self.contig_name = ''
        self.hit_contigs = set()
        self.hits_obj_list = []
        self.collapsed_hits_dict = defaultdict(list)

    def check_consistency(self):
        # if len(self.hit_contigs) == 1:
        consistency_check_list = []
        for hit in self.hits_obj_list:
            consistency_check_list.append("".join([hit.ref_contig, hit.orientation]))
        print set(consistency_check_list)

    def collapse_alignments(self):
        # if len(self.hit_contigs) == 1:
        consistency_check_list = []
        collapsed_ref = ''
        collapsed_start = ''
        collapsed_stop = ''
        collapsed_orientation = ''
        for hit in self.hits_obj_list:
            if (hit.ref_contig == collapsed_ref) and (hit.orientation == collapsed_orientation):
                # print "monkeys"
                if hit.orientation == '+':
                    collapsed_hit_obj.ref_stop = hit.ref_stop
                    collapsed_hit_obj.query_stop = hit.query_stop
                if hit.orientation == '-':
                    collapsed_hit_obj.ref_start = hit.ref_start
                    collapsed_hit_obj.query_stop = hit.query_start
                    # print "\t".join([hit.ref_start, hit.ref_stop])
            else:
                # print "badgers"
                if collapsed_ref != '':
                    collapsed_hit_obj.query_length = int(collapsed_hit_obj.query_stop) - int(collapsed_hit_obj.query_start)
                    collapsed_hit_obj.ref_length = int(collapsed_hit_obj.ref_stop) - int(collapsed_hit_obj.ref_start)
                    self.collapsed_hits_dict[hit.ref_contig].append(collapsed_hit_obj)
                collapsed_ref = hit.ref_contig
                collapsed_orientation = hit.orientation
                collapsed_hit_obj = Hits_obj(
                                    "\t".join([
                                            hit.query_contig,
                                            hit.query_start,
                                            hit.query_stop,
                                            hit.ref_contig,
                                            hit.ref_start,
                                            hit.ref_stop,
                                            "NA",
                                            hit.orientation
                                            ]))
        collapsed_hit_obj.query_length = int(collapsed_hit_obj.query_stop) - int(collapsed_hit_obj.query_start)
        collapsed_hit_obj.ref_length = int(collapsed_hit_obj.ref_stop) - int(collapsed_hit_obj.ref_start)
        self.collapsed_hits_dict[hit.ref_contig].append(collapsed_hit_obj)



class Hits_obj(object):
    """
    Information on a hit:
        query_contig:
        ref_contig:
        query_location: set of contig names
        ref_location:
        orientation:
    """

    def __init__(self, line):
        """"""
        split_line = line.split()
        self.query_contig = split_line[0]
        self.query_start = split_line[1]
        self.query_stop = split_line[2]
        self.ref_contig = split_line[3]
        self.ref_start = split_line[4]
        self.ref_stop = split_line[5]
        self.identity = split_line[6]
        self.orientation = split_line[7]
        self.query_length = int()
        self.ref_length = int()

#-----------------------------------------------------
# Step 3
# Import variables & load input files
#-----------------------------------------------------


contig_dict = defaultdict(list)

for line in satsuma_lines:
    line = line.rstrip()
    # print line
    split_line = line.split()
    query = split_line[0]
    ref = split_line[3]
    if contig_dict[query]:
        contig_obj = contig_dict[query][0]
    else:
        contig_obj = Contig_obj()
        contig_obj.contig_name = query
    contig_obj.hit_contigs.add(ref)
    hit_obj = Hits_obj(line)
    contig_obj.hits_obj_list.append(hit_obj)
    contig_dict[query] = [contig_obj]

# print contig_dict

scaffold_obj_dict = defaultdict(list)
for key in contig_dict.keys():
    # print key
    # print contig_dict[key]
    contig_obj = contig_dict[key][0]
    # print contig_obj
    hits = list(contig_obj.hit_contigs)
    # print ", ".join(hits_set)
    # hits =
    # hits = ", ".join()
    # print("\t".join([key, ", ".join(hits)]))
    # locations = []
    # for hit in contig_obj.hits_obj_list:
        # locations.append(hit.ref_start)
    # print "\t".join(locations)
    # contig_obj.check_consistency()
    contig_obj.collapse_alignments()
    for key in contig_obj.collapsed_hits_dict.keys():
        collapsed_hits = contig_obj.collapsed_hits_dict[key]
        for hit in collapsed_hits:
            print "\t".join([hit.query_contig, hit.query_start, hit.query_stop, str(hit.query_length), hit.ref_contig, hit.ref_start, hit.ref_stop, hit.orientation, str(hit.ref_length)])
            if scaffold_obj_dict[hit.ref_contig]:
                scaffold_obj = scaffold_obj_dict[hit.ref_contig][0]
            else:
                scaffold_obj = Scaffold_obj()
                scaffold_obj.ref_contig = hit.ref_contig
            scaffold_obj.contig_obj_dict[hit.ref_start] = [hit]
            scaffold_obj_dict[hit.ref_contig] = [scaffold_obj]

for scaffold in scaffold_obj_dict.keys():
    scaffold_obj = scaffold_obj_dict[scaffold][0]
    start_pos_list = scaffold_obj.contig_obj_dict.keys()
    start_pos_sorted = sorted([int(x) for x in start_pos_list])
    for start in start_pos_sorted:
        # print start
        hit = scaffold_obj.contig_obj_dict[str(start)][0]
        print "\t".join([hit.query_contig, hit.query_start, hit.query_stop, str(hit.query_length), hit.ref_contig, hit.ref_start, hit.ref_stop, hit.orientation, str(hit.ref_length)])
