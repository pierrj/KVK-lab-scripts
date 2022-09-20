#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
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

import pandas as pd

table_master = pd.read_csv('shuffled_1/venn_diagram_stats_out.tsv',sep='\t',header=None)
table_master.columns = ['overlap', 'biorep', 'rep1','rep2', 'rep3','rep1+2','rep2+3','rep1+3','rep1+2+3']
table_master = table_master.drop(columns=['rep1','rep2', 'rep3','rep1+2','rep2+3','rep1+3','rep1+2+3'])

rep_columns = []
for i in range(1,101):
    shuffled_table = pd.read_csv('shuffled_'+str(i)+'/venn_diagram_stats_out.tsv',sep='\t',header=None)
    shuffled_table.columns = ['overlap', 'biorep', 'rep1','rep2', 'rep3','rep1+2','rep2+3','rep1+3','rep1+2+3']
    shuffled_table['sum'] = shuffled_table[['rep1','rep2', 'rep3','rep1+2','rep2+3','rep1+3','rep1+2+3']].sum(axis=1)
    shuffled_table.loc[shuffled_table['biorep'] != 'G3_3','percentage'] = shuffled_table['rep1+2+3']/shuffled_table['sum']
    shuffled_table.loc[shuffled_table['biorep'] == 'G3_3','percentage'] = shuffled_table['rep1+2']/shuffled_table['sum']
    rep_columns.append('rep_'+str(i))
    table_master['rep_'+str(i)] = shuffled_table['percentage']

table_master['shuffled_mean'] = table_master[rep_columns].mean(axis=1)
table_master['shuffled_std'] = table_master[rep_columns].std(axis=1)
table_master = table_master.drop(columns=rep_columns)

table_master.to_csv('venn_diagram_stats_out.shuffled_all.tsv', sep='\t')