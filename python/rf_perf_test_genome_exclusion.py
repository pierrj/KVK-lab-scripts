import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
import sys

input_df = sys.argv[1]
majority_fraction = float(sys.argv[2])
approach = sys.argv[3]
n_estimators = int(sys.argv[4])
min_samples_split = int(sys.argv[5])
min_samples_leaf = int(sys.argv[6])
max_features = sys.argv[7]
max_depth = sys.argv[8]
bootstrap = eval(sys.argv[9])


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

def reports(model, X_test, y_test):
    y_pred = model.predict(X_test)
    TP = len(y_pred[(y_pred == 1) & (y_test == 1)])
    FN = len(y_pred[(y_pred == 0) & (y_test == 1)])
    FP = len(y_pred[(y_pred == 1) & (y_test == 0)])
    # sensitivity, how sensitive is the test? TP/TP+FN aka recall
    recall = TP/(TP+FN)
    ## PPV, how powerful is a positive? TP/TP+FP aka precision
    precision = TP/(TP+FP)
    ap = average_precision_score(y_test, model.predict_proba(X_test)[:,1])
    auc = roc_auc_score(y_test, model.predict_proba(X_test)[:,1])
    return([recall, precision, ap, auc])

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
    pav_true_subset = df_genes[df_genes['lineage_pav']==True].id
    if majority_fraction != 1.0:
        pav_false_subset_downsampled = np.random.choice(df_genes[df_genes['lineage_pav'] == False].id, size=int(len(df_genes.index)*majority_fraction),replace=False)
        df_genes_downsampled = df_genes[(df_genes.id.isin(pav_false_subset_downsampled)) | (df_genes.id.isin(pav_true_subset))]
    else:
        df_genes_downsampled = df_genes
    # drop columns
    df_genes_downsampled = df_genes_downsampled.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                            'genome', 'lineage', 'lineage_conserved', 'proportion'], axis=1)
    df_genes_test_subset = df_genes_test_subset.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                            'genome', 'lineage', 'lineage_conserved', 'proportion'], axis=1)
    y_train = df_genes_downsampled['lineage_pav']
    X_train = df_genes_downsampled.drop('lineage_pav', axis=1)
    y_test = df_genes_test_subset['lineage_pav']
    X_test = df_genes_test_subset.drop('lineage_pav', axis=1)
    return(y_train,X_train,y_test,X_test)

if approach not in ["RF", "SMOTE", "BRFC", "RF_balanced", "RF_balanced_subsample"]:
    raise NameError('Approach not in approved list')

column_names = ["recall", "precision", "ap", "auc"]
df_results = pd.DataFrame(columns = column_names)

y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)
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
model.fit(X_train, y_train)
results = reports(model, X_test, y_test)

# n_estimators = int(sys.argv[4])
# min_samples_split = int(sys.argv[5])
# min_samples_leaf = int(sys.argv[6])
# max_features = sys.argv[7]
# max_depth = sys.argv[8]
# bootstrap = eval(sys.argv[9])

print(approach + '\t' + 
            str(majority_fraction) + '\t' +
            str(n_estimators) + '\t' +
            str(min_samples_split) + '\t' +
            str(min_samples_leaf) + '\t' +
            str(max_features) + '\t' +
            str(max_depth) + '\t' +
            str(bootstrap) + '\t' +
            str(results[0]) + '\t' + 
            str(results[1]) + '\t' + 
            str(results[2]) + '\t' + 
            str(results[3]))