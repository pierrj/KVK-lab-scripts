import pandas as pd
import sys

input_file = sys.argv[1]
output_name = sys.argv[2]

df = pd.read_csv(input_file, sep='\t')

df_means = dict(zip(df.columns.to_list(), df.mean()))
df_means = pd.DataFrame(df_means, index=['1'])

df_means.to_csv(output_name,sep='\t', index=False)