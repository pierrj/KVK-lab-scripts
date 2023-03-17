#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
from BCBio import GFF
import os
import shutil
import sys
from os.path import join

genome_list = sys.argv[1]
# lineage_info_file = sys.argv[2]

# get list of genomes
with open(genome_list) as file:
    genomes = [genomes.strip() for genomes in file]

### read in lineage info file to add lineage info to protein names (only needed for rice blast)
## for gladieux proteomes/rice
# lineage_info = {}
# with open(lineage_info_file, newline = '') as file:
#     file_reader = csv.reader(file, delimiter = '\t')
#     for row in file_reader:
#         if row[0] == 'WD-3-1_1':
#             lineage_info['WD-3-1'] = row[2]+'_'+row[3]
#         else:
#             lineage_info[row[0]] = row[2]+'_'+row[3]

out_dir = 'all_gffs_fixed'
if os.path.isdir(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)

## read in all gffs, change the names and output processed gffs
for genome in genomes:
    # originally this for wheat
    accession = genome.split('_')[0] + genome.split('_')[1]
    ## but this for rice
    # isolate = genome
    # lineage = lineage_info[isolate]
    # accession = isolate + '_' + lineage
    ## and this for guy11
    # accession = genome
    in_file = join(genome, 'fungap_out','fungap_out', 'fungap_out.nocomments.gff3')
    in_handle = open(in_file)
    out_file = join(out_dir,genome+"_fungap_out.fixed.gff3")
    with open(out_file, "w") as out_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                # fix names so that the genome name is included in the gene name
                if feature.type != "remark":
                    try:
                        original_feature_id = feature.id[:]
                        gene_number = str(int(original_feature_id.split('_')[1].lstrip("0"))-1)
                        original_feature_id = original_feature_id.split('_')[0] + '_' + gene_number ## remove zeroes...
                        feature.id = original_feature_id + '_' + accession
                        feature.qualifiers['ID'] = [original_feature_id + '_' + accession]
                        feature.qualifiers['Name'] = [original_feature_id + '_' + accession]
                        feature.sub_features[0].id = [original_feature_id + '_' + accession+'T0']
                        feature.sub_features[0].qualifiers['ID'] = [original_feature_id + '_' + accession+'T0']
                        feature.sub_features[0].qualifiers['Parent'] = [original_feature_id + '_' + accession]
                        for sub_feature in feature.sub_features[0].sub_features:
                            sub_feature.id = [original_feature_id + '_' + accession+'T0']
                            sub_feature.qualifiers['Parent'] = [original_feature_id + '_' + accession+'T0']
                            if len(sub_feature.qualifiers['ID']) > 1:
                                print('sub_feature too long')
                            sub_feature.qualifiers['ID'] = [original_feature_id + '_' + accession+'T0'+'.'+sub_feature.qualifiers['ID'][0].split('.')[-1]]
                    except KeyError:
                        print(feature)
            GFF.write([rec], out_handle)
    in_handle.close()