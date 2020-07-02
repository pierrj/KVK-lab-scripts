import csv
import ipyparallel as ipp
import sys

split_read_file = str(sys.argv[1])

outwardfacing_read_file = str(sys.argv[2])

output_name = str(sys.argv[3])

scaffold_number = int(sys.argv[4])

print('successfully load modules')

# open putative ecc list and index to speed up confirming eccs
with open(split_read_file, newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    eccloc_list = []
    for row in file_reader:
        ecc_loc = [int(row[0]) - 1, int(row[1]), int(row[2])]
        eccloc_list.append(ecc_loc)

print('successfully opened split read file')

# open opposite facing discordant read file
with open(outwardfacing_read_file, newline = '') as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    # index discordant read file so that confirmeccs() only looks at discordant reads found on the same chromosome
    discordant_indexed = [[] for i in range(scaffold_number)]
    for row in discordant_reader:
        discordant_indexed[(int(row[0])-1)].append([int(row[1]), int(row[2])])

print('successfully opened files')

def confirmeccs(ecc):
    for i in range(0, len(discordant_indexed[ecc[0]]), 2):
        reada = discordant_indexed[ecc[0]][i]
        readb = discordant_indexed[ecc[0]][i+1]
        # on the fly sorting because the discordant reads are next to each other but not sorted by base position
        if reada[0] < readb[0]:
            read1 = reada
            read2 = readb
        else:
            read1 = readb
            read2 = reada
        if ecc[1] <= read1[0] <= ecc[2] and ecc[1] <= read1[1] <= ecc[2] and ecc[1] <= read2[0] <= ecc[2] and ecc[1] <= read2[1] <= ecc[2]:
            # calculate total distance from ecc start and end
            distance = abs(read1[1] - ecc[1]) + abs(read2[0] - ecc[2])
            # set here for insert size distribution
            if distance <= 500:
                return True
    return False

rc = ipp.Client(profile='pierrj', cluster_id = "cluster-id-memory-test")
dview = rc[:]
dview.block = True
lview = rc.load_balanced_view()
lview.block = True

mydict = dict(discordant_indexed = discordant_indexed)
dview.push(mydict)

print('successfully pushed')

# get true/false list if each ecc is confirmed, then compress only keeps where true is in the list
yesornoeccs = list(lview.map(confirmeccs, eccloc_list))

print('succesfully confirmed')

# confirmedeccs = list(compress(eccloc_list, yesornoeccs))