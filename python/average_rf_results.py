import pandas as pd
import sys

input_file = sys.argv[1]
output_name = sys.argv[2]

df = pd.read_csv(input_file, sep='\t', header = None)
df.columns = ['approach', 'majority_fraction', 'recall', 'precision', 'ap', 'auc']

df_grouped = df.groupby(["approach", "majority_fraction"]).mean().reset_index()

df_grouped.to_csv(output_name + "average_results.txt",sep='\t', index=False)