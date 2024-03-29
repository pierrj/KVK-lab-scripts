import matplotlib.pyplot as plt
from statistics import mean
from matplotlib import pyplot
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import plot_confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE
import sys

input_table = sys.argv[1]
test_table = sys.argv[2]
output_string = sys.argv[3]

df_genes = pd.read_csv(input_table)

# remove lineage 4 because there are so few genomes in that lineage
df_genes = df_genes[df_genes['lineage']!=4]

df_genes = df_genes.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                        'genome', 'lineage', 'lineage_conserved', 'proportion',
                        'LTR/Gypsy', 'Unknown', 'DNA',
                        'DNA/TcMar-Fot1', 'LINE/Tad1', 'DNA/Tc-Mar',
                        'LTR/Copia','DNA/MULE-MuDR','DNA/hAT-Ac',
                        'DNA/CMC-EnSpm','LINE/CRE', 'DNA/Kolobok-H', 
                        'LTR/Unknown','LTR/Pao','DNA/TcMar-Pogo',
                        'LINE/R2-NeSL','LINE/Penelope'], axis=1)


y = df_genes['lineage_pav']
X = df_genes.drop('lineage_pav', axis=1)

#Use SMOTE to oversample the minority class
oversample = SMOTE()
over_X, over_y = oversample.fit_resample(X, y)
over_X_train, over_X_test, over_y_train, over_y_test = train_test_split(over_X, over_y, test_size=0.1, stratify=over_y)

#Build SMOTE SRF model
SMOTE_SRF = RandomForestClassifier(random_state=0)
#Create Stratified K-fold cross validation
cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
scoring = ('f1', 'recall', 'precision')
#Evaluate SMOTE SRF model
scores = cross_validate(SMOTE_SRF, over_X, over_y, scoring=scoring, cv=cv)
#Get average evaluation metrics
print('Mean f1: %.3f' % mean(scores['test_f1']))
print('Mean recall: %.3f' % mean(scores['test_recall']))
print('Mean precision: %.3f' % mean(scores['test_precision']))

#Randomly spilt dataset to test and train set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, stratify=y)
#Train SMOTE SRF
SMOTE_SRF.fit(over_X_train, over_y_train)
#SMOTE SRF prediction result
y_pred = SMOTE_SRF.predict(X_test)

fig = plot_confusion_matrix(SMOTE_SRF, X_test, y_test, display_labels=['Will not be lost', 'Will be lost'], cmap='Greens')
plt.title('SMOTE + Standard Random Forest Confusion Matrix')
plt.tight_layout()
plt.savefig(output_string + '_confusion_matrix_internal_test.png',facecolor='w',dpi=200)

TP = len(y_pred[(y_pred == 1) & (y_test == 1)])
FN = len(y_pred[(y_pred == 0) & (y_test == 1)])
TN = len(y_pred[(y_pred == 0) & (y_test == 0)])
FP = len(y_pred[(y_pred == 1) & (y_test == 0)])

# specificity, how specific is the test? TN/TN+FP
print('specificity')
print(TN/(TN+FP))

# sensitivity, how sensitive is the test? TP/TP+FN aka recall
print('sensitivity')
print(TP/(TP+FN))

## PPV, how powerful is a positive? TP/TP+FP aka precision
print('PPV')
print(TP/(TP+FP))

## NPV, how powerful is a negative? TN/TN+FN
print('NPV')
print(TN/(TN+FN))

df_genes_other_host = pd.read_csv(test_table)

# remove lineage 4 because there are so few genomes in that lineage
df_genes_other_host = df_genes_other_host[df_genes_other_host['lineage']!=4]

df_genes_other_host = df_genes_other_host.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                        'genome', 'lineage', 'lineage_conserved', 'proportion',
                        'LTR/Gypsy', 'Unknown', 'DNA',
                        'DNA/TcMar-Fot1', 'LINE/Tad1', 'DNA/Tc-Mar',
                        'LTR/Copia','DNA/MULE-MuDR','DNA/hAT-Ac',
                        'DNA/CMC-EnSpm','LINE/CRE', 'DNA/Kolobok-H', 
                        'LTR/Unknown','LTR/Pao','DNA/TcMar-Pogo',
                        'LINE/R2-NeSL','LINE/Penelope'], axis=1)


y_other_host = df_genes_other_host['lineage_pav']
X_other_host = df_genes_other_host.drop('lineage_pav', axis=1)

y_pred_other_host = SMOTE_SRF.predict(X_other_host)

print('preds: ' + output_string)
TP = len(y_pred_other_host[(y_pred_other_host == 1) & (y_other_host == 1)])
FN = len(y_pred_other_host[(y_pred_other_host == 0) & (y_other_host == 1)])
TN = len(y_pred_other_host[(y_pred_other_host == 0) & (y_other_host == 0)])
FP = len(y_pred_other_host[(y_pred_other_host == 1) & (y_other_host == 0)])

# specificity, how specific is the test? TN/TN+FP
print('specificity')
print(TN/(TN+FP))

# sensitivity, how sensitive is the test? TP/TP+FN aka recall
print('sensitivity')
print(TP/(TP+FN))

## PPV, how powerful is a positive? TP/TP+FP aka precision
print('PPV')
print(TP/(TP+FP))

## NPV, how powerful is a negative? TN/TN+FN
print('NPV')
print(TN/(TN+FN))