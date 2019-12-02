import csv
import ipyparallel as ipp
from itertools import groupby

with open('merged.splitreads.sorted.reverseread1.G3_1A_bwamem.bed', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    eccloc_list = [[int(row[0][10:12]), int(row[1]), int(row[2])] for row in eccloc_reader]
    
with open('discordantmappedreads.oppositefacing.bed', newline = '') as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    discordant_list = [[int(row[0][10:12]), int(row[1]), int(row[2])] for row in discordant_reader]
    discordant_indexed = [[x[1:] for x in g] for k, g in groupby(discordant_list, key = lambda x: x[0])]

def confirmeccs(ecc):
    for i in range(0, len(discordant_indexed[ecc[0]]), 2):
        read1 = discordant_indexed[ecc[0]][i]
        read2 = discordant_indexed[ecc[0]][i+1]
        if ecc[1] <= read1[0] <= ecc[2] and ecc[1] <= read1[1] <= ecc[2] and ecc[1] <= read2[0] <= ecc[2] and ecc[1] <= read2[1] <= ecc[2]:
            return True
    return False

rc = ipp.Client(profile='default', cluster_id='')
dview = rc[:]
dview.block = True
lview = rc.load_balanced_view()
lview.block = True

mydict = dict(discordant_list = discordant_list)
dview.push(mydict)

yesornoeccs = list(lview.map(confirmeccs, eccloc_list))
from itertools import compress
confirmedeccs = list(compress(eccloc_list, yesornoeccs))

with open('parallel.confirmed', 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    w.writerows(confirmedeccs)
