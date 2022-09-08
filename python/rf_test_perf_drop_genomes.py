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
output_string = sys.argv[2]

df_genes = pd.read_csv(input_table)

# remove lineage 4 because there are so few genomes in that lineage
df_genes = df_genes[df_genes['lineage']!=4]

# pick 5 genomes for final testing (in a way that scientists will understand better)
genome_test_subset = np.random.choice(np.unique(df_genes.genome), size=5,replace=False)
df_genes_test_subset = df_genes[df_genes.genome.isin(genome_test_subset)]
df_genes = df_genes[~df_genes.genome.isin(genome_test_subset)]

# drop columns
df_genes = df_genes.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                        'genome', 'lineage', 'lineage_conserved', 'proportion'], axis=1)

y = df_genes['lineage_pav']
X = df_genes.drop('lineage_pav', axis=1)

#Use SMOTE to oversample the minority class
oversample = SMOTE()
over_X, over_y = oversample.fit_resample(X, y)
over_X_train, over_X_test, over_y_train, over_y_test = train_test_split(over_X, over_y, test_size=0.1, stratify=over_y)

#Build SMOTE SRF model
SMOTE_SRF = RandomForestClassifier(n_estimators=900, # default is 100
                                min_samples_split=2, # default is 2
                                min_samples_leaf=1, # default is 1
                                max_features=None, # default is sqrt
                                max_depth=60, # default is none
                                bootstrap=True, # default is True...
                                random_state=1)
#Create Stratified K-fold cross validation
cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
scoring = ('f1', 'recall', 'precision')
#Evaluate SMOTE SRF model
scores = cross_validate(SMOTE_SRF, over_X, over_y, scoring=scoring, cv=cv, n_jobs=16,pre_dispatch='n_jobs')

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

for genome in genome_test_subset:
    df_genes_test_subset_per_genome = df_genes_test_subset[df_genes_test_subset['genome'] == genome]
    df_genes_test_subset_per_genome = df_genes_test_subset_per_genome.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                            'genome', 'lineage', 'lineage_conserved', 'proportion'], axis=1)
    truths = df_genes_test_subset_per_genome['lineage_pav']
    genome_X_test = df_genes_test_subset_per_genome.drop('lineage_pav', axis=1)
    preds = SMOTE_SRF.predict(genome_X_test)
    TP = len(y_pred[(preds == 1) & (truths == 1)])
    FN = len(y_pred[(preds == 0) & (truths == 1)])
    TN = len(y_pred[(preds == 0) & (truths == 0)])
    FP = len(y_pred[(preds == 1) & (truths == 0)])
    print(genome)
    print('specificity')
    print(TN/(TN+FP))
    print('sensitivity')
    print(TP/(TP+FN))
    print('PPV')
    print(TP/(TP+FP))
    print('NPV')
    print(TN/(TN+FN))