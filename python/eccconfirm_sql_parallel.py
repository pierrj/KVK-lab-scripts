
# coding: utf-8

# In[1]:


import csv
import ipyparallel as ipp
from itertools import compress
import sqlite3


# In[2]:


def data_entry(data):
    c.execute('''Drop TABLE if exists server''')
    c.execute('''Create TABLE if not exists server(chrom, start, end)''')
    for i in range(len(data)):
        c.execute("INSERT INTO server(chrom, start, end) VALUES(?,?,?)", (test[i]))
    conn.commit()


# In[3]:


def confirmeccs(ecc):
    c.execute("SELECT * FROM server")
    for read1 in c:
        read2 = next(c)
        if read1[0]==ecc[0]==read2[0]:
            if ecc[1] <= read1[1] <= ecc[2] and ecc[1] <= read1[2] <= ecc[2] and ecc[1] <= read2[1] <= ecc[2] and ecc[1] <= read2[2] <= ecc[2]:
                return True
    return False


# In[ ]:


with open('merged.splitreads.sorted.reverseread1.G3_1A_bwamem.bed', newline = '') as eccloc:
    eccloc_reader = csv.reader(eccloc, delimiter = '\t')
    eccloc_list = [[int(row[0][10:12]), int(row[1]), int(row[2])] for row in eccloc_reader]
    
with open('discordantmappedreads.oppositefacing.bed', newline = '') as discordant:
    discordant_reader = csv.reader(discordant, delimiter = '\t')
    discordant_list = [[int(row[0][10:12]), int(row[1]), int(row[2])] for row in discordant_reader]


# In[ ]:


conn = sqlite3.connect(r"eccloclist_sql.db")
c = conn.cursor()
data_entry(discordant_list)


# In[ ]:


rc = ipp.Client(profile='default', cluster_id='')
dview = rc[:]
dview.block = True
lview = rc.load_balanced_view()
lview.block = True


# In[ ]:


dview.execute('import sqlite3')
dview.execute('conn = sqlite3.connect(r"eccloclist_sql.db")')
dview.execute('c = conn.cursor()')


# In[ ]:


yesornoeccs = list(lview.map(confirmeccs, eccloc_list))
confirmedeccs = list(compress(eccloc_list, yesornoeccs))


# In[ ]:


with open('parallel.confirmed', 'w', newline = '') as confirmed:
    w = csv.writer(confirmed, delimiter = '\t')
    w.writerows(confirmedeccs)


# In[ ]:


dview.execute('conn.close()')
conn.close()

