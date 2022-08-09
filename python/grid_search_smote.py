from statistics import mean
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SMOTE

df_genes = pd.read_csv('gene_info_w_proportions.csv')

# remove lineage 4
df_genes = df_genes[df_genes['lineage']!=4]

df_genes = df_genes.drop(['id', 'scaffold', 'start', 'end', 'orientation', 'orthogroups', 'enough_space_te', 'enough_space_gene',
                        'genome', 'lineage', 'lineage_conserved', 'proportion'], axis=1)

y = df_genes['lineage_pav']
X = df_genes.drop('lineage_pav', axis=1)

#Use SMOTE to oversample the minority class
oversample = SMOTE()
over_X, over_y = oversample.fit_resample(X, y)
over_X_train, over_X_test, over_y_train, over_y_test = train_test_split(over_X, over_y, test_size=0.1, stratify=over_y)

n_estimators = [int(x) for x in np.linspace(start = 800, stop = 1200, num = 9)]
print(n_estimators)
max_features = ["sqrt", "log2", None]
max_depth = [int(x) for x in np.linspace(30, 70, num = 5)]
print(max_depth)
min_samples_split = [2]
min_samples_leaf = [1]
bootstrap = [True,False]
# Create the random grid
param_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap}

# Use the random grid to search for best hyperparameters
# First create the base model to tune
rf = RandomForestClassifier()
# Random search of parameters, using 3 fold cross validation, 
# search across 100 different combinations, and use all available cores
rf_random = GridSearchCV(estimator = rf, param_grid = param_grid, cv = 3, verbose=2, n_jobs = 16, pre_dispatch='n_jobs')
# Fit the random search model
rf_random.fit(over_X_train, over_y_train)

print(rf_random.best_params_)