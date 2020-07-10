import ipyparallel as ipp
from itertools import compress

list1 = [ [i, i, i] for i in range(4000000)]
list2 = [ [i, i, i] for i in range(2000000, 6000000)]

print('generated list')

def loop(item):
    lview.results.clear()
    for i in range(len(list2)):
        if list2[i][0] == item[0]:
            return True
    return False

rc = ipp.Client(profile='default', cluster_id = "cluster-id-dummy")
dview = rc[:]
dview.block = True
lview = rc.load_balanced_view()
lview.block = True

print('connected to engines')

mydict = dict(list2 = list2)
dview.push(mydict)

print('pushed list 2')

trueorfalse = list(lview.map(loop, list1))

print('looped through list 2')

# found = list(compress(list1, trueorfalse))

# print('compressed lists')