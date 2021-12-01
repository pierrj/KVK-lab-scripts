from Bio import SeqIO
import sys
import os
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

msa_dir = sys.argv[1]
genomes = sys.argv[2]
output = sys.argv[3]

msa_list = os.listdir(msa_dir)

## make dict of sequence lists
genomes_dict = {}

with open(genomes, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        genomes_dict[row[0].split("_")[0]] = SeqRecord(Seq(""), id=row[0])


for msa in msa_list:
    msa_path = msa_dir+ '/' +msa
    for record in SeqIO.parse(msa_path, 'fasta'):
        genome = record.id.split("_")[2]
        genomes_dict[genome].seq += record.seq

with open(output, 'w') as handle:
    SeqIO.write(genomes_dict.values(), handle, 'fasta')