import operator
import os
import pickle
import re
import warnings
from argparse import ArgumentParser
from collections import defaultdict
from functools import reduce
from BCBio import GFF
from Bio import BiopythonWarning, SeqIO
warnings.simplefilter('ignore', BiopythonWarning)
SHORT_LEN = 10
prefix = os.path.basename(os.path.splitext(gff3_file)[0])
genome_assembly_file = os.path.abspath('/Users/pierrj/fungap_runs/guy11/guy11_genome_baoetal2017.fasta')
in_seq_handle = open(genome_assembly_file)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, 'fasta'))
in_seq_handle.close()


gff3_file = os.path.abspath('/Users/pierrj/fungap_runs/guy11/fungap_out/augustus_out/augustus.gff3')
gff3_file = os.path.abspath('/Users/pierrj/fungap_runs/guy11/fungap_out/maker_out/SRR8842990/maker_SRR8842990.gff3')
gff3_file = os.path.abspath('/Users/pierrj/fungap_runs/guy11/fungap_out/braker_out/SRR8842990/braker_SRR8842990.gff3')

gff3_file = os.path.abspath('/Users/pierrj/fungap_runs/guy11/fungap_out/braker_out/SRR8842990/braker_makerscript.gff3')




in_handle = open(gff3_file)
for rec in GFF.parse(in_handle, base_dict=seq_dict):
    gene_features = rec.features
    for gene_feature in gene_features:
        mrna_features = gene_feature.sub_features
        print('printing mrna_features')
        print(mrna_features)
        for mrna_feature in mrna_features:
            print('printing mrna_feature')
            print(mrna_feature)
            mrna_sub_features = mrna_feature.sub_features
            print('printing mrna_subfeatures')
            print(mrna_sub_features)
            mrna_sub_features_s = sorted(
                mrna_sub_features, key=lambda x: x.location.start
            )
            seq_cds = []
            coords = []
            mrna_sub_features_s2 = []
            for feature in mrna_sub_features_s:
                print('printing feature')
                print(feature)
                if feature.type != 'CDS':
                    continue
                mrna_sub_features_s2.append(feature)
                seq_cds.append(rec.seq[
                    feature.location.start:
                    feature.location.end])
                coords.append(
                    (feature.location.start, feature.location.end)
                )
            i = 1
            while i < len(coords):
                intron_start = coords[i - 1][1]
                intron_end = coords[i][0]
                intron_len = intron_end - intron_start
                if intron_len < 10:
                    d_bad[(prefix, mrna_feature.id)] = True
                    d_intron[prefix] += 1
                i += 1
            gene_seq = reduce(operator.add, seq_cds)


gff3_file = os.path.abspath('/Users/pierrj/fungap_runs/guy11/fungap_out/maker_out/SRR8842990/maker_SRR8842990.gff3')
in_handle = open(gff3_file)
for rec in GFF.parse(in_handle, base_dict=seq_dict):
    gene_features = rec.features
    print(gene_features)

gff3_file = os.path.abspath('/Users/pierrj/fungap_runs/guy11/fungap_out/braker_out/SRR8842990/braker_SRR8842990.gff3')
in_handle = open(gff3_file)
for rec in GFF.parse(in_handle, base_dict=seq_dict):
    gene_features = rec.features
    for gene_feature in gene_features:
        mrna_features = gene_feature.sub_features
        print(mrna_features)