import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from rfpimp import importances
from rfpimp import dropcol_importances


df_genes = pd.read_csv('all_genes_distance_and_gc.csv')

## convert true/false to int

df_genes["te_nearby"] = df_genes["te_nearby"].astype(int)
df_genes["gene_nearby"] = df_genes["gene_nearby"].astype(int)
df_genes["lineage_conserved"] = df_genes["lineage_conserved"].astype(int)

df_genes = df_genes[['te_nearby','gene_nearby','lineage_conserved', 'gene_gc', 'flanking_1kb_gc', 'flanking_5kb_gc']]

y = df_genes['lineage_conserved']
x = df_genes.drop('lineage_conserved', axis=1)

x_train, x_test, y_train, y_test = train_test_split(x,y,test_size=0.1,random_state=101)

model = RandomForestClassifier()
model.fit(x_train,y_train)

print('score')
print(model.score(x_train,y_train))
print(model.score(x_test,y_test))

predictions = model.predict(x_test)

TP = len(predictions[(predictions == 1) & (y_test == 1)])
FN = len(predictions[(predictions == 0) & (y_test == 1)])
TN = len(predictions[(predictions == 0) & (y_test == 0)])
FP = len(predictions[(predictions == 1) & (y_test == 0)])

# specificity, how specific is the test? TN/TN+FP
print('specificity')
print(TN/(TN+FP))

# sensitivity, how sensitive is the test? TP/TP+FN
print('sensitivity')
print(TP/(TP+FN))

## PPV, how powerful is a positive? TP/TP+FP
print('PPV')
print(TP/(TP+FP))

## NPV, how powerful is a negative? TN/TN+FN
print('NPV')
print(TN/(TN+FN))

# if we called everything positive
print('if we just called everything positive')
TP = len(y_test[y_test == 1])
FN = 0
TN = 0
FP = len(y_test[y_test == 0])

print('specificity')
print(TN/(TN+FP))
print('sensitivity')
print(TP/(TP+FN))
print('PPV')
print(TP/(TP+FP))
# print('NPV')
# print(TN/(TN+FN))

## default importance
I = pd.DataFrame()
I['Feature'] = x_train.columns
I['Importance'] = model.feature_importances_
I = I.sort_values('Importance', ascending=False)
I = I.set_index('Feature')
I.to_csv('default_importances.txt',sep='\t')

# permutation importances
I = importances(model, x_test, y_test)
I.to_csv('permutation_importances.txt',sep='\t')

# # drop col importances
# I = dropcol_importances(model, x_train, y_train, x_test, y_test)
# I.to_csv('dropcol_importances.txt', sep='\t')