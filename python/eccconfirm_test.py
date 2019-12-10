import csv
import ipyparallel as ipp
from itertools import groupby
from itertools import compress
import subprocess
import statistics
import collections

with open('/global/home/users/pierrj/testfiles/500.testbedtools.processed', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    eccloc_list = [[int(row[0][10:12]) - 1, int(row[1]), int(row[2])] for row in eccloc_reader]

with open('/global/home/users/pierrj/testfiles/sorted.discordantmappedreads.oppositefacing.bed', newline = '') as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    discordant_indexed = [[] for i in range(56)]
    for row in discordant_reader:
        discordant_indexed[(int(row[0][10:12])-1)].append([int(row[1]), int(row[2])])

rc = ipp.Client(profile='default', cluster_id='')
dview = rc[:]
dview.block = True
lview = rc.load_balanced_view()
lview.block = True

mydict = dict(discordant_indexed = discordant_indexed)
dview.push(mydict)

def confirmeccs(ecc):
    for i in range(0, len(discordant_indexed[ecc[0]]), 2):
        read1 = discordant_indexed[ecc[0]][i]
        read2 = discordant_indexed[ecc[0]][i+1]
        if ecc[1] <= read1[0] <= ecc[2] and ecc[1] <= read1[1] <= ecc[2] and ecc[1] <= read2[0] <= ecc[2] and ecc[1] <= read2[1] <= ecc[2]:
            return True
    return False

def confirmeccs_nomap(eccloc):
    confirmedeccs = []
    for l in range(len(eccloc)):
        ecc = eccloc[l]
        for i in range(0, len(discordant_indexed[ecc[0]]), 2):
            read1 = discordant_indexed[ecc[0]][i]
            read2 = discordant_indexed[ecc[0]][i+1]
            if ecc[1] <= read1[0] <= ecc[2] and ecc[1] <= read1[1] <= ecc[2] and ecc[1] <= read2[0] <= ecc[2] and ecc[1] <= read2[1] <= ecc[2]:
                confirmedeccs.append(ecc)
                break
    return confirmedeccs

yesornoeccs = list(lview.map(confirmeccs, eccloc_list))
confirmedeccs = list(compress(eccloc_list, yesornoeccs))
confirmedeccs_nomap = confirmeccs_nomap(eccloc_list)

if confirmedeccs == confirmedeccs_nomap:
    print('true')
else:
    print('false')