import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import sys
import pickle

input_table = sys.argv[1]
rf_pkl = sys.argv[2]
output_table = sys.argv[3]

df_genes = pd.read_csv(input_table)

params = [
    'any_te',
    'gene_nearby',
    'gene_gc',
    'flanking_1kb_gc',
    'lengths',
    'tm',
    'signalp',
    'effectorp',
    'H3K27ac',
    'H3K27me3',
    'H3K36me3',
    'cm_expression',
    'ip_expression',
    'eccdna_cov',
    'methylation',
    'go',
    'pfam'
]

df_genes = df_genes[params]

file = open(rf_pkl, 'rb')
model = pickle.load(file)

y_pred = model.predict(df_genes)

df_genes = pd.read_csv(input_table)
df_genes['predicted_lineage_pav'] = y_pred

df_genes.to_csv(output_table, index=False)