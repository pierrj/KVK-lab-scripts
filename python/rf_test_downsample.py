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
from sklearn.metrics import classification_report
from imblearn.over_sampling import SMOTE
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
import sys
import pickle
from sklearn.metrics import PrecisionRecallDisplay

def reports_short(model, X_test, y_test, majority_fraction, replicate):
    y_pred = model.predict(X_test)
    TP = len(y_pred[(y_pred == 1) & (y_test == 1)])
    FN = len(y_pred[(y_pred == 0) & (y_test == 1)])
    TN = len(y_pred[(y_pred == 0) & (y_test == 0)])
    FP = len(y_pred[(y_pred == 1) & (y_test == 0)])
    # sensitivity, how sensitive is the test? TP/TP+FN aka recall
    recall = TP/(TP+FN)
    ## PPV, how powerful is a positive? TP/TP+FP aka precision
    precision = TP/(TP+FP)
    ap = average_precision_score(y_test, model.predict_proba(X_test)[:,1])
    auc = roc_auc_score(y_test, model.predict_proba(X_test)[:,1])
    return([majority_fraction, replicate, recall, precision, ap, auc])

def train_test_split_mine_downsample(majority_fraction):
    df_genes = pd.read_csv('gene_info_w_proportions.csv')
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

column_names = ["majority_fraction", "replicate", "recall", "precision", "ap", "auc"]
df_srf = pd.DataFrame(columns = column_names)
df_brfc = pd.DataFrame(columns = column_names)
df_smotesrf = pd.DataFrame(columns = column_names)
df_srf_balanced = pd.DataFrame(columns = column_names)
df_srf_balanced_subsample = pd.DataFrame(columns = column_names)

reps = 3
fractions = [0.05, 0.1, 0.25,0.5,0.75,0.95] # max is 0.99

for majority_fraction in fractions:
    for replicate in range(reps):
        y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)
        SRF = RandomForestClassifier()
        SRF.fit(X_train, y_train)
        row = reports_short(SRF, X_test, y_test, majority_fraction, replicate)
        df_srf.loc[len(df_srf.index)] = row
    for replicate in range(reps):
        y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)
        BRFC = BalancedRandomForestClassifier()
        BRFC.fit(X_train, y_train)
        row = reports_short(BRFC, X_test, y_test, majority_fraction, replicate)
        df_brfc.loc[len(df_brfc.index)] = row
    for replicate in range(reps):
        y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)
        oversample = SMOTE()
        over_X_train, over_y_train = oversample.fit_resample(X_train, y_train)
        SMOTE_SRF = RandomForestClassifier()
        SMOTE_SRF.fit(over_X_train, over_y_train)
        row = reports_short(SMOTE_SRF, X_test, y_test, majority_fraction, replicate)
        df_smotesrf.loc[len(df_smotesrf.index)] = row
    for replicate in range(reps):
        y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)
        SRF_balanced = RandomForestClassifier(class_weight="balanced")
        SRF_balanced.fit(X_train, y_train)
        row = reports_short(SRF_balanced, X_test, y_test, majority_fraction, replicate)
        df_srf_balanced.loc[len(df_srf_balanced.index)] = row
    for replicate in range(reps):
        y_train,X_train,y_test,X_test = train_test_split_mine_downsample(majority_fraction)
        SRF_balanced_subsample = RandomForestClassifier(class_weight="balanced_subsample")
        SRF_balanced_subsample.fit(X_train, y_train)
        row = reports_short(SRF_balanced_subsample, X_test, y_test, majority_fraction, replicate)
        df_srf_balanced_subsample.loc[len(df_srf_balanced_subsample.index)] = row

df_srf.to_csv('srf_downsample_test_all_tes.csv')
df_brfc.to_csv('brfc_downsample_test_all_tes.csv')
df_smotesrf.to_csv('smotesrf_downsample_test_all_tes.csv')
df_srf_balanced.to_csv('srf_balanced_downsample_test_all_tes.csv')
df_srf_balanced_subsample.to_csv('srf_balanced_subsample_downsample_test_all_tes.csv')