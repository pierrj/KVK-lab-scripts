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
import sys
import os
from Bio import SeqIO

input_dir = sys.argv[1]
output_file = sys.argv[2]

# read in all of the orthogroup fasta files and concatenate them all together
with open(output_file, 'w') as output:
    for og_fasta in os.listdir(input_dir):
        orthogroup = og_fasta[0:9]
        record_dict = SeqIO.to_dict(SeqIO.parse(input_dir+'/'+og_fasta, "fasta"))
        for record in record_dict:
            record_id = record_dict[record].id
            new_id = record_dict[record].id + '_' + orthogroup
            record_dict[record].id = new_id
            SeqIO.write(record_dict[record], output, 'fasta')