import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
import sys

input_df = sys.argv[1]
reps = int(sys.argv[2])
majority_fraction = float(sys.argv[3])
approach = sys.argv[4]

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
    gene_test_subset = np.random.choice(df_genes.id, size=int(len(df_genes.index)*0.1),replace=False)
    df_genes_test_subset = df_genes[df_genes.id.isin(gene_test_subset)]
    df_genes = df_genes[~df_genes.id.isin(gene_test_subset)]
    pav_true_subset = df_genes[df_genes['lineage_pav']==True].id
    pav_false_subset_downsampled = np.random.choice(df_genes[df_genes['lineage_pav'] == False].id, size=int(len(df_genes.index)*majority_fraction),replace=False)
    df_genes_downsampled = df_genes[(df_genes.id.isin(pav_false_subset_downsampled)) | (df_genes.id.isin(pav_true_subset))]
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

if majority_fraction > 0.95:
    raise ValueError('Majority fraction has to be less than 0.95 otherwise it isnt compatible with downsampling')

column_names = ["recall", "precision", "ap", "auc"]
df_results = pd.DataFrame(columns = column_names)

for replicate in range(reps):
    y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)
    if approach == "SMOTE":
        oversample = SMOTE()
        over_X_train, over_y_train = oversample.fit_resample(X_train, y_train)
        X_train = over_X_train
        y_train = over_y_train
    if approach == "BRFC":
        model = BalancedRandomForestClassifier()
    elif approach == "RF_balanced":
        model = RandomForestClassifier(class_weight="balanced")
    elif approach == "RF_balanced_subsample":
        model = RandomForestClassifier(class_weight="balanced_subsample")
    elif approach == "RF":
        model = RandomForestClassifier()
    elif approach == "SMOTE":
        model = RandomForestClassifier()
    model.fit(X_train, y_train)
    row = reports(model, X_test, y_test)
    df_results.loc[len(df_results.index)] = row

averages = df_results.mean().tolist()
print(approach + '\t' + str(majority_fraction) + '\t' + str(averages[0]) + '\t' + str(averages[1]) + '\t' + str(averages[2]) + '\t' + str(averages[3]))