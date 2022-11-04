import pandas as pd
import sys

input_file = sys.argv[1]
output_name = sys.argv[2]

df = pd.read_csv(input_file, sep='\t', header = None)
df.columns = ['approach', 'majority_fraction', 
            'n_estimators', 'min_samples_split',
            'min_samples_leaf', 'max_features',
            'max_depth', 'bootstrap',
            'recall', 'precision', 'ap', 'auc',
            'TP', 'FN', 'FP', 'TN']


df_grouped = df.groupby(['approach', 'majority_fraction', 
            'n_estimators', 'min_samples_split',
            'min_samples_leaf', 'max_features',
            'max_depth', 'bootstrap']).mean().reset_index()

df_grouped.to_csv(output_name,sep='\t', index=False)