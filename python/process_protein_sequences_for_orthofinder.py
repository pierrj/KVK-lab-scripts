from Bio import SeqIO
import sys
import os
import shutil
import csv

seq_dir = sys.argv[1]
out_dir = sys.argv[2]
lineage_info_file = sys.argv[3]
output_accessions = sys.argv[4]


lineage_info = {}
with open(lineage_info_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    for row in file_reader:
        if row[0] == 'WD-3-1_1':
            lineage_info['WD-3-1'] = row[2]+'_'+row[3]
        else:
            lineage_info[row[0]] = row[2]+'_'+row[3]

seq_list = os.listdir(seq_dir)

if os.path.isdir(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)

accessions = []

for seq in seq_list:
    print(seq)
    if '_fungap_out_prot.faa' in seq:
        accession = seq.split('_')[0] + seq.split('_')[1]
    elif seq == 'GCA_000002495.2_MG8_protein.faa':
        accession = '70-15'
    elif seq == 'GCA_004355905.1_PgNI_protein.faa':
        accession = 'NI907'
    elif '_protein.fasta' in seq:
        isolate = seq.split('_')[0]
        lineage = lineage_info[isolate]
        accession = seq.split('_')[0] + '_' + lineage
    else:
        print('couldnt process accession')
    print(accession)
    accessions.append(accession)
    out_file = out_dir + '/' + accession + '_proteome_filtered.fasta'
    seq_path = seq_dir + '/' + seq
    record_list = list(SeqIO.parse(seq_path, 'fasta'))
    with open(out_file, 'w') as corrected:
        for i in range(len(record_list)):
            record = record_list[i]
            record.id = 'gene_' + str(i) + '_' + accession ## rename records to have genome name in them
            record.description = ''
            if '*' in record.seq:
                if record.seq[-1] == '*': ## remove stop codon from end of sequences
                    record.seq = record.seq[:-1]
                    SeqIO.write(record, corrected, 'fasta')
            else:
                SeqIO.write(record, corrected, 'fasta')

with open(output_accessions, 'w', newline = '') as output_csv:
    w = csv.writer(output_csv, delimiter = '\t')
    for row in accessions:
        w.writerow([row])