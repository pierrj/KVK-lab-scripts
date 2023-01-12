#MIT License
#
#Copyright (c) 2023 Pierre Michel Joubert
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
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import f1_score
import sys
import matplotlib.pyplot as plt

# input parameters
input_df = sys.argv[1]
majority_fraction = float(sys.argv[2])
approach = sys.argv[3]
n_estimators = int(sys.argv[4])
min_samples_split = int(sys.argv[5])
min_samples_leaf = int(sys.argv[6])
max_features = sys.argv[7]
max_depth = sys.argv[8]
bootstrap = eval(sys.argv[9])

# make sure the parameters are properly formated
def none_or_str(value):
    if value == 'None':
        return None
    return value

max_features = none_or_str(max_features)


def none_or_int(value):
    if value == 'None':
        return None
    return int(value)

max_depth = none_or_int(max_depth)

args_dict = {
    "n_estimators": n_estimators,
    "min_samples_split": min_samples_split,
    "min_samples_leaf": min_samples_leaf,
    "max_features": max_features,
    "max_depth": max_depth,
    "bootstrap": bootstrap
}

# make my own train test split command that keeps random genomes out of the data rather than random genes
def train_test_split_mine_downsample(majority_fraction):
    df_genes = pd.read_csv(input_df)
    df_genes = df_genes[df_genes['lineage']!=4]
    ## pick 4 genomes per lineage as testing data
    genome_test_subset = []
    for lineage in np.unique(df_genes.lineage):
        for genome in np.random.choice(df_genes[df_genes.lineage == lineage].genome, size=4,replace=False):
            genome_test_subset.append(genome)
    df_genes_test_subset = df_genes[df_genes.genome.isin(genome_test_subset)]
    df_genes = df_genes[~df_genes.genome.isin(genome_test_subset)]
    # downsample non-pav genes based off majority fraction parameter
    if majority_fraction != 1.0:
        pav_true_subset = df_genes[df_genes['lineage_pav']==True].id
        pav_false_subset_downsampled = np.random.choice(df_genes[df_genes['lineage_pav'] == False].id, size=int(len(df_genes.index)*majority_fraction),replace=False)
        df_genes_downsampled = df_genes[(df_genes.id.isin(pav_false_subset_downsampled)) | (df_genes.id.isin(pav_true_subset))]
    else:
        df_genes_downsampled = df_genes
    # drop extra columns
    df_genes_downsampled = df_genes_downsampled.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                            'genome', 'lineage', 'lineage_conserved', 'proportion'], axis=1)
    df_genes_test_subset = df_genes_test_subset.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                            'genome', 'lineage', 'lineage_conserved', 'proportion'], axis=1)
    y_train = df_genes_downsampled['lineage_pav']
    X_train = df_genes_downsampled.drop('lineage_pav', axis=1)
    y_test = df_genes_test_subset['lineage_pav']
    X_test = df_genes_test_subset.drop('lineage_pav', axis=1)
    return(y_train,X_train,y_test,X_test)

y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)

# different models based off input selection
if approach == "SMOTE":
    oversample = SMOTE()
    over_X_train, over_y_train = oversample.fit_resample(X_train, y_train)
    X_train = over_X_train
    y_train = over_y_train
if approach == "BRFC":
    model = BalancedRandomForestClassifier(**args_dict)
elif approach == "RF_balanced":
    model = RandomForestClassifier(class_weight="balanced", **args_dict)
elif approach == "RF_balanced_subsample":
    model = RandomForestClassifier(class_weight="balanced_subsample", **args_dict)
elif approach == "RF":
    model = RandomForestClassifier(**args_dict)
elif approach == "SMOTE":
    model = RandomForestClassifier(**args_dict)

# fit model
model.fit(X_train, y_train)

# get baseline f1
y_pred = model.predict(X_test)
baseline = f1_score(y_test, y_pred)
permuted_diffs = []
for column in X_test.columns:
    # get f1 score when X_test column is permuted, output f1 and then reset column
    save = X_test[column].copy()
    X_test[column] = np.random.permutation(X_test[column])
    y_pred = model.predict(X_test)
    permuted_f1 = f1_score(y_test, y_pred)
    diff = baseline-permuted_f1
    permuted_diffs.append(diff)
    X_test[column] = save

## print variable names and importances
print('\t'.join(map(str, X_test.columns.to_list())))

print('\t'.join(map(str, permuted_diffs)))