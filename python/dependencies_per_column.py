import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import f1_score
from sklearn.metrics import r2_score
from pandas.api.types import is_bool_dtype
import sys

input_xtrain = sys.argv[1]
input_xtest = sys.argv[2]
input_dep_feature = sys.argv[3]
majority_fraction = float(sys.argv[4])
approach = sys.argv[5]
n_estimators = int(sys.argv[6])
min_samples_split = int(sys.argv[7])
min_samples_leaf = int(sys.argv[8])
max_features = sys.argv[9]
max_depth = sys.argv[10]
bootstrap = eval(sys.argv[11])

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

## make sure colnames are okay
X_train = pd.read_csv(input_xtrain)
X_test = pd.read_csv(input_xtest)

all_cols = [col for col in X_train]
boolcols = [col for col in X_train if is_bool_dtype(X_train[col])]

dep_feature = input_dep_feature
row = [dep_feature]
X_dep_train, y_dep_train = X_train.drop(dep_feature, axis=1), X_train[dep_feature]
X_dep_test, y_dep_test = X_test.drop(dep_feature, axis=1), X_test[dep_feature]
if dep_feature in boolcols:
    model_dep = RandomForestClassifier(**args_dict)
else:
    model_dep = RandomForestRegressor(**args_dict)
model_dep.fit(X_dep_train,y_dep_train)
y_dep_pred = model_dep.predict(X_dep_test)
if dep_feature in boolcols:
    baseline = f1_score(y_dep_test, y_dep_pred)
else:
    baseline = r2_score(y_dep_test, y_dep_pred)
row.append(baseline)
for perm_feature in X_train.columns:
    if perm_feature == dep_feature:
        row.append('x')
        continue
    save = X_dep_test[perm_feature].copy()
    X_dep_test[perm_feature] = np.random.permutation(X_dep_test[perm_feature])
    y_dep_pred_permuted = model_dep.predict(X_dep_test)
    if dep_feature in boolcols:
        permuted_score = f1_score(y_dep_test, y_dep_pred_permuted)
    else:
        permuted_score = r2_score(y_dep_test, y_dep_pred_permuted)
    diff = baseline-permuted_score
    row.append(diff)
    X_dep_test[perm_feature] = save

print('\t'.join(map(str, row)))