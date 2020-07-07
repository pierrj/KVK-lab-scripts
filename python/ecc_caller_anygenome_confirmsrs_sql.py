import csv
import ipyparallel as ipp
import sys
import sqlite3
from collections import Counter
from itertools import compress

split_read_file = str(sys.argv[1])

outwardfacing_read_file = str(sys.argv[2])

output_name = str(sys.argv[3])

scaffold_number = int(sys.argv[4])

# open putative ecc list and index to speed up confirming eccs
with open(split_read_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    eccloc_list = []
    for row in file_reader:
        ecc_loc = [int(row[0]) - 1, int(row[1]), int(row[2])]
        eccloc_list.append(ecc_loc)

# open opposite facing discordant read file as SQL databases
with open(outwardfacing_read_file) as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    for row in discordant_reader:
        for i in range(scaffold_number):
            conn = sqlite3.connect(r"scaffold"+str(i+1)+"_sql.db")
            c = conn.cursor()
            c.execute('''Drop TABLE if exists server''')
            c.execute('''Create TABLE if not exists server(start, end, name)''')
            to_add = []
            while int(row[0]) == i+1:
                to_add.append((int(float(row[1])), int(float(row[2])), str(row[3]))) ### read names here too
                try:
                    row = next(discordant_reader)
                except StopIteration:
                    break
            c.executemany("INSERT INTO server(start, end, name) VALUES(?,?,?)", to_add)
            conn.commit()
            conn.close()

# does proximity filtering based off an estimated insert size of 400 + 25%

def confirmeccs(ecc):
    conn = sqlite3.connect(r"scaffold"+str(ecc[0]+1)+"_sql.db")
    c = conn.cursor()
    c.execute("SELECT * FROM server WHERE start >= "+str(ecc[1])+" AND start <= " +str(ecc[2])+ " AND end >= "+str(ecc[1])+" AND end <= " + str(ecc[2]) )
    opposite_read_names = []
    opposite_reads = []
    for row in c.fetchall():
        opposite_read_names.append(row[2])
        opposite_reads.append((row[0], row[1], row[2]))
    conn.close()
    counted = Counter(opposite_read_names) 
    filtered_names = set([el for el in opposite_read_names if counted[el] >= 2])
    filtered_reads = [[row[0], row[1], row[2]] for row in opposite_reads if row[2] in filtered_names]
    filtered_reads.sort(key = lambda x: x[2])
    for i in range(0, len(filtered_reads), 2):
        reada = opposite_reads[i]
        readb = opposite_reads[i+1]
        if reada[0] < readb[0]:
            read1 = reada
            read2 = readb
        else:
            read1 = readb
            read2 = reada
        if ecc[1] <= read1[0] <= ecc[2] and ecc[1] <= read1[1] <= ecc[2] and ecc[1] <= read2[0] <= ecc[2] and ecc[1] <= read2[1] <= ecc[2]:
            distance = abs(read1[1] - ecc[1]) + abs(read2[0] - ecc[2])
            if distance <= 500:
                return True
    return False

# open parallelization client
rc = ipp.Client(profile='default', cluster_id = "cluster-id-" + output_name)
dview = rc[:]
dview.block = True
lview = rc.load_balanced_view()
lview.block = True

dview.execute('from collections import Counter')
dview.execute('import sqlite3')

# get true/false list if each ecc is confirmed, then compress only keeps where true is in the list
yesornoeccs = list(lview.map(confirmeccs, eccloc_list))
confirmedeccs = list(compress(eccloc_list, yesornoeccs))

# write confirmed eccs to file
with open('parallel.confirmed', 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    w.writerows(confirmedeccs)