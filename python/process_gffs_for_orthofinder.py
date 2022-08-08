from BCBio import GFF
import csv
import os
import shutil
import sys
from os import join

genome_list = sys.argv[1]

with open(genome_list) as file:
    genomes = [genomes.strip() for genomes in file]

out_dir = 'all_gffs_fixed'
if os.path.isdir(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)

for genome in genomes:
    accession = genome.split('_')[0] + genome.split('_')[1]
    print(genome)
    in_file = join(genome, 'fungap_out','fungap_out', 'fungap_out', 'fungap_out.gff3')
    in_handle = open(in_file)
    out_file = join(out_dir,genome+"fungap_out.fixed.gff3")
    with open(out_file, "w") as out_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                try:
                    original_feature_id = feature.id[:]
                    feature.id = original_feature_id + '_' + accession
                    feature.qualifiers['ID'] = [original_feature_id + '_' + accession]
                    feature.qualifiers['Name'] = [original_feature_id + '_' + accession]
                    feature.sub_features[0].id = [original_feature_id + '_' + accession+'T0']
                    feature.sub_features[0].qualifiers['ID'] = [original_feature_id + '_' + accession+'T0']
                    feature.sub_features[0].qualifiers['Parent'] = [original_feature_id + '_' + accession]
                    for sub_feature in feature.sub_features[0].sub_features:
                        sub_feature.id = [original_feature_id + '_' + accession+'T0']
                        sub_feature.qualifiers['Parent'] = [original_feature_id + '_' + accession+'T0']
                except KeyError:
                    print(feature)
            GFF.write([rec], out_handle)
    in_handle.close()