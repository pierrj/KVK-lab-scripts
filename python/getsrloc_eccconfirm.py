import csv
import ipyparallel as ipp
from itertools import groupby
from itertools import compress
import subprocess


## length filter junctions
with open('samechromosome.exactlytwice.reverseread1.G3_1A_bwamem.bam', newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open('getsrloc_test', 'w', newline = '') as filtered:
        w = csv.writer(filtered, delimiter = '\t')
        for line1 in file_reader:
            line2 = next(file_reader)
            sum = int(line1[3]) - int(line2[3])
            if abs(sum) <= 1000000:
                w.writerow(line1)
                w.writerow(line2)

filename = 'getsrloc_test'
filtered_filename = 'qualityfiltered.' + filename
exactlytwice_filename = 'exactlytwice.' + filtered_filename
actualbam_filename = 'actualbam.' + exactlytwice_filename
bedfile = filename + '.bed'
samtools = '/global/home/groups/consultsw/sl-7.x86_64/modules/samtools/1.8/bin/samtools'
bedtools = '/global/home/groups/consultsw/sl-7.x86_64/modules/bedtools/2.28.0/bin/bedtools'

quality_filter = '''awk '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\\\1", "", $6); if((a !~ /[DMIHS]/ && int(a) > 49 ) || (b !~ /[DMIHS]/ && int(b) > 49)) print $0}' ''' + filename + ''' > ''' + filtered_filename
exactlytwice_filter = '''awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' ''' + filtered_filename + ''' ''' + filtered_filename + ''' > ''' + exactlytwice_filename
make_bam = 'bash -c \"' + samtools + ' view -b -h <(cat <(' + samtools + ' view -H G3_1A_bwamem.bam) ' + exactlytwice_filename +') > ' + actualbam_filename + '\"'
bamtobed_sort = bedtools + ' bamtobed -i ' + actualbam_filename + ' | sort -k 4,4 -k 2,2 > ' + bedfile

subprocess.run(quality_filter, shell= True)
subprocess.run(exactlytwice_filter, shell= True)
subprocess.run(make_bam, shell= True)
subprocess.run(bamtobed_sort, shell= True)

with open('getsrloc_test.bed', newline = '') as file:
    file_reader = csv.reader(file, delimiter = '\t')
    with open('merged.getsrloc_test', 'w', newline = '') as filtered:
        w = csv.writer(filtered, delimiter = '\t')
        for line1 in file_reader:
            line2 = next(file_reader)
            line1[2] = line2[1]
            w.writerow(line1)

with open('merged.getsrloc_test', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    eccloc_list = [[int(row[0][10:12]) - 1, int(row[1]), int(row[2])] for row in eccloc_reader]
    
with open('discordantmappedreads.oppositefacing.bed', newline = '') as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    discordant_indexed = [[] for i in range(56)]
    for row in discordant_reader:
        discordant_indexed[(int(row[0][10:12])-1)].append([int(row[1]), int(row[2])])

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

mydict = dict(discordant_indexed = discordant_indexed)
dview.push(mydict)

yesornoeccs = list(lview.map(confirmeccs, eccloc_list))
confirmedeccs = list(compress(eccloc_list, yesornoeccs))

with open('parallel.confirmed.getsrlocmergetest', 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    w.writerows(confirmedeccs)
