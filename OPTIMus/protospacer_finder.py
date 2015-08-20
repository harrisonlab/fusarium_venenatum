#!/usr/bin/python

'''
This tool can be used to add a note attribute to a gff feature in a
gffutils database. Gff features matching a given ID will have notes
added to them. IDs can be searched in other attributes as well such
as Names using the "--attribute" option
'''

import sys,argparse
import gffutils
from itertools import chain


#######################################
#           Load sys. args.           #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('--in_db',required=True,type=str,help='The gffutils database to be searched')
ap.add_argument('--out_db',required=True,type=str,help='The name of the output gff file')
ap.add_argument('--id_file',required=True,type=str,help='A text file containing gene IDs or transcript IDs for annotation')
ap.add_argument('--str',required=True,type=str,help='Text to add to the notes section for each matched feature.')
ap.add_argument('--attribute',required=False,default='ID',type=str,help='The attribute to search for the IDs within, default is ID')

conf = ap.parse_args() #sys.argv


#######################################
#        Create output database       #
#                                     #
#                                     #
#######################################

f = conf.in_db
in_db = gffutils.FeatureDB(f)


for gene in in_db:
	gene_name
	contig_name
	gene_sense
	for exon in gene:
		seq = fasta_file[exon_coordinates]
		search_for_protospacer in seq
		protospacer_start
		protospacer_stop
		start_in_contig
		stop_in_contig
		protospacer_strand
		protospacaer_name = gene_name + itteration
		protospacer_sense
		print(gene_name	contig_name	protospacer_name protspacer_seq	strand_start	strand_stop)
