#!/usr/bin/python

'''
This program is used to build information on all the genes predicted in
an annotated genome. These commands take information on location of genes
& suppliment this information with information on interproscan domains
and swissprot annotations.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
import numpy as np
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--genome',required=True,type=str,help='A fasta file of the assembled contigs')
ap.add_argument('--genes_gff',required=True,type=str,help='A gff file of the genes from the FoC')
ap.add_argument('--InterPro',required=True,type=str,help='The Interproscan functional annotation .tsv file')
ap.add_argument('--Swissprot',required=True,type=str,help='A parsed table of BLAST results against the Swissprot database. Note - must have been parsed with swissprot_parser.py')
ap.add_argument('--Antismash',required=True,type=str,help='Output of Antismash parsed into a tsv file of gene names, contig, secmet function and cluster ID')
ap.add_argument('--Smurf',required=True,type=str,help='Output of Smurf parsed into a tsv file of gene names, contig, secmet function and cluster ID')
ap.add_argument('--vitamins',required=True,type=str,help='Table of blast hits produced by Greg showing hits from Fg homologs identified as vitamin pathway genes')
ap.add_argument('--TFs',required=True,type=str,help='Tab seperated of putative transcription factors and their domains as identified by interpro2TFs.py')
ap.add_argument('--orthogroups_PH1', required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--orthogroups_GR1', required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--DEGs',required=True,type=str,help='List of the total set of differentially expressed genes')
ap.add_argument('--fpkm',required=True,type=str,help='Calculated fpkm per condition for each gene in tsv format')

# ap.add_argument('--Fv_OrthoMCL_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
# ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')
# ap.add_argument('--OrthoMCL_path',required=True,type=str,nargs='+',help='The identifiers of pathogenic strains used in the orthology analysis')
# ap.add_argument('--OrthoMCL_nonpath',required=True,type=str,nargs='+',help='The identifiers of non-pathogenic strains used in the orthology analysis')



conf = ap.parse_args()

with open(conf.genome) as f:
    contig_lines = f.readlines()
with open(conf.genes_gff) as f:
    unordered_gff_lines = f.readlines()
with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()
with open(conf.Swissprot) as f:
    swissprot_lines = f.readlines()
with open(conf.Antismash) as f:
    antismash_lines = f.readlines()
with open(conf.Smurf) as f:
    smurf_lines = f.readlines()
with open(conf.vitamins) as f:
    vitamin_lines = f.readlines()
with open(conf.TFs) as f:
    TF_lines = f.readlines()
with open(conf.orthogroups_PH1) as f:
    PH1_orthogroup_lines = f.readlines()
with open(conf.orthogroups_GR1) as f:
    GR1_orthogroup_lines = f.readlines()
with open(conf.DEGs) as f:
    DEG_lines = f.readlines()
with open(conf.fpkm) as f:
    fpkm_lines = f.readlines()

column_list=[]


#-----------------------------------------------------
# Define classes
#-----------------------------------------------------

class Expression_obj(object):
    """Information on a genes fpkm under various conditions and whether a gene
        is differentially expressed under different conditions:
    """
    def __init__(self, gene_id):
        """Return a Expression_obj whose name is *gene_id*"""
        self.gene_id = gene_id
        # self.DEG_conditions = set()
        self.DEG_conditions = str()
        # self.fpkm_dict = defaultdict(list)
        self.condition_list = []
        self.fpkm_list = []
        self.treatments = []
        self.mean_fpkm = []
        self.FC_values = []
        self.FC_values_adj = []
        self.FC_dict = defaultdict(str)
        self.DEG_conditions_adj = []

    def add_DEG(self, DEG_condition):
        """"""
        # self.DEG_conditions.add(DEG_condition)
        self.DEG_conditions = DEG_condition
        # print DEG_condition

    def filter_DEG(self):
        """"""
        if self.DEG_conditions != "":
            conditions = self.DEG_conditions
            split_condits = conditions.split(";")
            # print split_condits
            for condit in split_condits:
                FC = self.FC_dict[condit]
                # print self.gene_id
                # print condit
                # print FC
                if float(FC) >= 2:
                    self.DEG_conditions_adj.append(("DEG:" + condit).upper())
                if float(FC) <= -2:
                    self.DEG_conditions_adj.append(("DEG:" + condit).lower())


    def add_fpkm(self, conditions_list, fpkm_list):
        """"""
        for condition, fpkm in zip(conditions_list, fpkm_list):
            self.condition_list.append(condition)
            fpkm = int(np.round_(float(fpkm),  decimals=0))
            self.fpkm_list.append(str(fpkm))

    def add_LFC(self):
        "Add LFC information based upon fpkm data already in the dictionary"
        fpkm_dict = defaultdict(list)
        for condition_rep, fpkm in zip(self.condition_list, self.fpkm_list):
            condition = condition_rep[:-2]
            fpkm_dict[condition].append(int(fpkm))
        fpkm_list = []
        condition_list = []
        for condition in fpkm_dict.keys():
            # condition_list.append(condition)
            fpkm_mean = np.mean(fpkm_dict[condition])
            fpkm_dict[condition].append(fpkm_mean)
            self.treatments.append(condition)
            self.mean_fpkm.append(str(np.round_(fpkm_mean,  decimals=1)))
        LFC_conditions = [
        ('Control', '02793'),
        ('Control', 'F55'),
        ('Control', '10170'),
        ('Control', 'MWT'),
        ('Control', 'MOL'),
        ('Control', 'MKO'),
        ('Control', 'TJ')
        ]
        for x, y in LFC_conditions:
            x_fpkm_mean = fpkm_dict[x][-1]
            y_fpkm_mean = fpkm_dict[y][-1]
            if y_fpkm_mean == x_fpkm_mean:
                FC = 1
                FC_adj = 0
            elif y_fpkm_mean >= x_fpkm_mean:
                FC = np.divide(y_fpkm_mean, x_fpkm_mean)
            else:
                FC = np.negative(np.divide(x_fpkm_mean, y_fpkm_mean))
            if x_fpkm_mean < 1:
                x_adj = 1
            else:
                x_adj = x_fpkm_mean
            if y_fpkm_mean < 1:
                y_adj = 1
            else:
                y_adj = y_fpkm_mean
            FC_adj = np.log2(y_adj) - np.log2(x_adj)
            self.FC_values.append(str(np.round_(FC,  decimals=1)))
            self.FC_values_adj.append(str(np.round_(FC_adj,  decimals=1)))
            self.FC_dict[y] = str(np.round_(FC_adj,  decimals=1))


#-----------------------------------------------------
# Step X
# Process fpkm data
#-----------------------------------------------------

expression_dict = defaultdict(list)
for line in DEG_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    condition = split_line[1]
    if expression_dict[gene_id]:
        exp_obj = expression_dict[gene_id][0]
    else:
        exp_obj = Expression_obj(gene_id)
    exp_obj.add_DEG(condition)
    # print "\t".join([gene_id, condition])
    # print exp_obj.DEG_conditions
    expression_dict[gene_id].append(exp_obj)

conditions_list = fpkm_lines[0].rstrip().split("\t")
conditions_list = conditions_list[2:]
# print conditions_list
for line in fpkm_lines[1:]:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    fpkm_list = split_line[2:]
    if expression_dict[gene_id]:
        expression_dict[gene_id][0].add_fpkm(conditions_list, fpkm_list)
    else:
        exp_obj = Expression_obj(gene_id)
        exp_obj.add_fpkm(conditions_list, fpkm_list)
        expression_dict[gene_id].append(exp_obj)

for gene_id in expression_dict.keys():
    for obj in expression_dict[gene_id]:
        # print (gene_id + "\t" + "\t".join(obj.fpkm_list))
        obj.add_LFC()
        # print obj.DEG_conditions
        obj.filter_DEG()
        # print (gene_id + "\t" + "\t".join(obj.FC_values_adj))


#-----------------------------------------------------
# Step 2
# Collect information on contig length from the genome
# assembly file. This can be used to determine if a
# gene has 2kb of sequence assembled up and downstream.
# This important for knock out design.
#-----------------------------------------------------

contig_len_dict = defaultdict(list)
contig_id = ""
seq_lines = ""
for line in contig_lines:
    line = line.rstrip()
    if line.startswith(">"):
        last_seq_length = len(seq_lines)
        contig_len_dict[contig_id] = len(seq_lines)
        split_line = line.split(" ")
        contig_id = split_line[0].replace(">", "")
        seq_lines = ""
    else:
        seq_lines += line



#-----------------------------------------------------
# Step 2.5
# Order gff features by input contigs and gene start
#-----------------------------------------------------

contig_list = []
gene_start_dict = defaultdict(list)
features_dict = defaultdict(list)
ordered_gff_lines = []

for line in unordered_gff_lines:
    line = line.strip("\n")
    if line.startswith("#"):
        ordered_gff_lines.append(line)
        continue
    split_line = line.split()
    if split_line[2] == 'gene':
        contig = split_line[0]
        gene_start = split_line[3]
        if contig not in contig_list:
            contig_list.append(contig)
        gene_start_dict[contig].append(int(gene_start))
        key = "_".join([contig, gene_start])
    features_dict[key].append(line)

gff_lines = []
for contig in sorted(contig_list, key = lambda x: (int(x.split('_')[1]))):
    gene_start_list = gene_start_dict[contig]
    for gene_start in sorted(gene_start_list):
        key = "_".join([contig, str(gene_start)])
        ordered_gff_lines.extend(features_dict[key])

#-----------------------------------------------------
# Step 3
# Append co-ordinates from the FoC gene gff, showing
# gene locations.
# Also identify whether there is 2kb sequence data up
# and downstream of the gene allowing design of
# knockouts
#-----------------------------------------------------

gene_id_set = Set([])
gene_id_list = []
FoC_genes_dict = defaultdict(list)
for line in ordered_gff_lines:
    if "gff-version" in line:
        continue
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split("\t")
    if 'mRNA' in split_line[2]:
        gene_features = split_line[8].split(';')
        gene_id = gene_features[0]
        gene_id = gene_id.replace('ID=', '')
        column_list = ["", "", "", ""]
        if gene_id not in gene_id_set:
            gene_id_list.append(gene_id)
            gene_id_set.add(gene_id)
        column_list=itemgetter(0, 3, 4, 6)(split_line)
        for column in column_list:
            FoC_genes_dict[gene_id].append(column)

        contig_id = column_list[0]
        feature_start=int(column_list[1])
        feature_end=int(column_list[2])
        contig_length = contig_len_dict[contig_id]
        if (feature_start - 2000) > 0 and (feature_end +2000) < contig_length:
            FoC_genes_dict[gene_id].append("Flank")
        else:
            FoC_genes_dict[gene_id].append("")


#-----------------------------------------------------
# Step 4
# Build a dictionary of interproscan annotations
# Annotations first need to be filtered to remove
# redundancy. This is done by first loading anntoations
# into a set.
#-----------------------------------------------------

interpro_set =  Set([])
interpro_dict = defaultdict(list)

for line in InterPro_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    interpro_columns = []
    index_list = [0, 4, 5, 11, 12]
    for x in index_list:
        if len(split_line) > x:
            interpro_columns.append(split_line[x])
    set_line = ";".join(interpro_columns)
    if set_line not in interpro_set:
        gene_id = interpro_columns[0]
        interpro_feat = ";".join(interpro_columns[1:])
        interpro_dict[gene_id].append(interpro_feat)
    interpro_set.add(set_line)


#-----------------------------------------------------
# Step 5
# Build a dictionary of Swissprot annotations
#-----------------------------------------------------

swissprot_dict = defaultdict(list)

for line in swissprot_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    swissprot_columns = itemgetter(14, 12, 13)(split_line)

    swissprot_dict[gene_id].extend(swissprot_columns)


#-----------------------------------------------------
# Step 5
# Build a dictionary of Secondary Metabolite annotations
#-----------------------------------------------------

antismash_dict = defaultdict(list)
for line in antismash_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    secmet_func = split_line[2]
    cluster = "AS_" + split_line[3]
    # secmet_columns = itemgetter(3, 4)(split_line)
    antismash_dict[gene_id].extend([secmet_func, cluster])

smurf_dict = defaultdict(list)
for line in smurf_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    secmet_func = split_line[2]
    cluster = "SM_" + split_line[3]
    # secmet_columns = itemgetter(3, 4)(split_line)
    smurf_dict[gene_id].extend([secmet_func, cluster])



#-----------------------------------------------------
# Step 6
# Build a dictionary of vitamin metabolism gene homologs
#-----------------------------------------------------

vitamin_dict = defaultdict(list)
for line in vitamin_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    Fg_homolog = split_line[2]
    kegg_term = split_line[3]
    vit_pathway = split_line[5]
    if vitamin_dict[gene_id]:
        vitamin_dict[gene_id][0]+=";" + str(vit_pathway)
        vitamin_dict[gene_id][1]+=";" + str(kegg_term)
        vitamin_dict[gene_id][2]+=";" + str(Fg_homolog)
    else:
        vitamin_dict[gene_id].extend([vit_pathway, kegg_term, Fg_homolog])

#-----------------------------------------------------
# Step 6
# Build a dictionary of vitamin metabolism gene homologs
#-----------------------------------------------------

TF_dict = defaultdict(list)
for line in TF_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    TF_function = split_line[2]
    TF_dict[gene_id].append(TF_function)


#-----------------------------------------------------
# Step 7
# Build a dictionary of orthogroups
#-----------------------------------------------------

strain_id = 'A3_5'
strain_id = strain_id + "|"

PH1_orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)
all_isolates = ['A3_5', 'PH1']
for line in PH1_orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        content_counts = ":".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes
    for gene_id in split_line[1:]:
        content_str = ",".join(split_line[1:])
        if strain_id in gene_id:
            gene_id = gene_id.replace(strain_id, '')
            PH1_orthogroup_dict[gene_id].extend([orthogroup_id, content_counts, content_str])

GR1_orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)
all_isolates = ['A3_5', 'FusGr1']
for line in GR1_orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        content_counts = ":".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes
    for gene_id in split_line[1:]:
        content_str = ",".join(split_line[1:])
        if strain_id in gene_id:
            gene_id = gene_id.replace(strain_id, '')
            GR1_orthogroup_dict[gene_id].extend([orthogroup_id, content_counts, content_str])


#-----------------------------------------------------
# Step 6
# Print final table of information on query, blast
# results and genes intersecting blast results
#-----------------------------------------------------

header_line = [
"gene_id", "contig", "gene_start", "gene_end", "gene_strand",
"2kb_flank",
"Cluster_ID", "SecMet_function", "SecMet_program", "Secmet_cluster",
"Vitamin_pathway", "KEGG_ID", "Fg_homolog",
"TF",
"Swissprot_organism", "Swissprot_hit", "Swissprot_function",
"Interpro_annotations",
"PH1_orthogroup", "PH1_genes_per_isolate", "PH1_orthogroup_contents",
"GR1_orthogroup", "GR1_genes_per_isolate", "GR1_orthogroup_contents"
]
header_line.extend(conditions_list)
header_line.append("DEG")
header_line.extend(['02793 LFC', 'F55 LFC', '10170 LFC', 'MWT LFC', 'MOL LFC', 'MKO LFC', 'TJ LFC'])
print ("\t".join(header_line))

in_cluster = False
cluster_num = 0
for gene_id in gene_id_list:
    useful_columns=[gene_id]
    useful_columns.extend(FoC_genes_dict[gene_id])
    if antismash_dict[gene_id] or smurf_dict[gene_id]:
        if in_cluster == False: # group antismash and smurf clusters into a single set of clusters
            cluster_num += 1
            cluster_name = "SecMet_cluster_" + str(cluster_num)
            in_cluster=True
        if antismash_dict[gene_id] and smurf_dict[gene_id]:
            # print gene_id
            antismash_cols = antismash_dict[gene_id]
            smurf_cols = smurf_dict[gene_id]
            prog = "antismash;smurf"
            secmet_func = ";".join([antismash_cols[0], smurf_cols[0]])
            cluster = ";".join([antismash_cols[1], smurf_cols[1]])
        elif antismash_dict[gene_id]:
            antismash_cols = antismash_dict[gene_id]
            prog = "antismash"
            secmet_func = antismash_cols[0]
            cluster = antismash_cols[1]
        elif smurf_dict[gene_id]:
            smurf_cols = smurf_dict[gene_id]
            prog = "smurf"
            secmet_func = smurf_cols[0]
            cluster = smurf_cols[1]
        else:
            continue
        useful_columns.extend([cluster_name, secmet_func, prog, cluster])
    else:
        useful_columns.extend(["","","",""])
        in_cluster=False

    key = gene_id.split(".")[0]
    if vitamin_dict[key]:
        vit_cols = vitamin_dict[key]
        useful_columns.extend(vit_cols)
    else:
        useful_columns.extend(["","",""])

    if TF_dict[gene_id]:
        TF_functions = TF_dict[gene_id]
        useful_columns.append(";".join(TF_functions))
    else:
        useful_columns.append("")

    if swissprot_dict[gene_id]:
        useful_columns.extend(swissprot_dict[gene_id])
    else:
        useful_columns.extend(["","",""])

    if interpro_dict[gene_id]:
        interpro_col = "|".join(interpro_dict[gene_id])
        useful_columns.append(interpro_col)
    else:
        useful_columns.append("")

    if PH1_orthogroup_dict[gene_id]:
        useful_columns.extend(PH1_orthogroup_dict[gene_id])
    else:
        useful_columns.extend(["", "", ""])

    if GR1_orthogroup_dict[gene_id]:
        useful_columns.extend(GR1_orthogroup_dict[gene_id])
    else:
        useful_columns.extend(["", "", ""])

    id = re.sub(r"\.t.*", "", gene_id)
    # print id
    for obj in expression_dict[id]:
        useful_columns.extend(obj.fpkm_list)
        useful_columns.append(",".join(obj.DEG_conditions_adj))
        useful_columns.extend(obj.FC_values_adj)

    print ("\t".join(useful_columns))
