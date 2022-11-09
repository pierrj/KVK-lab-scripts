import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import sys
import pickle

input_table = sys.argv[1]
rf_pkl = sys.argv[2]
output_table = sys.argv[3]

df_genes = pd.read_csv(input_table)
df_genes = df_genes.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                            'genome', 'lineage'], axis=1)


file = open(rf_pkl, 'rb')
model = pickle.load(file)

y_pred = model.predict(df_genes)

df_genes = pd.read_csv(input_table)
df_genes['predicted_lineage_pav'] = y_pred