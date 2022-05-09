import sys
import os
from Bio import SeqIO

input_dir = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, 'w') as output:
    for og_fasta in os.listdir(input_dir):
        orthogroup = og_fasta[0:9]
        record_dict = SeqIO.to_dict(SeqIO.parse(input_dir+'/'+og_fasta, "fasta"))
        for record in record_dict:
            record_id = record_dict[record].id
            new_id = record_dict[record].id + '_' + orthogroup
            record_dict[record].id = new_id
            SeqIO.write(record_dict[record], output, 'fasta')