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
from Bio import SeqIO
import sys
import os
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

msa_dir = sys.argv[1]
genomes = sys.argv[2]
output = sys.argv[3]

msa_list = sorted(os.listdir(msa_dir))

## make dict of sequence lists
genomes_dict = {}

with open(genomes, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        genomes_dict[row[0].split("_")[0]] = SeqRecord(Seq(""), id=row[0])

## concatenate all of the sequences within dictionary
for msa in msa_list:
    print(msa)
    msa_path = msa_dir+ '/' +msa
    for record in SeqIO.parse(msa_path, 'fasta'):
        genome = record.id.split("_")[2]
        genomes_dict[genome].seq += record.seq

for genome in genomes_dict.keys():
    print(genome)
    print(genomes_dict[genome].seq.count('-'))

# print concatenated sequences
with open(output, 'w') as handle:
    SeqIO.write(genomes_dict.values(), handle, 'fasta')