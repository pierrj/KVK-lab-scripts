import numpy as np
from itertools import product
from sklearn.model_selection import ParameterGrid

approach = ["RF", "SMOTE"]
majority_fraction = [0.25, 0.5, 1.0]
n_estimators = [int(x) for x in np.linspace(start = 1000, stop = 4000, num = 5)]
min_samples_split = [2, 5]
min_samples_leaf = [1]
max_features = ["sqrt", "log2", None]
max_depth = [int(x) for x in np.linspace(50, 200, num = 6)]
max_depth.append(None)
bootstrap = [True,False]
grid = product(approach,
            majority_fraction,
            n_estimators,
            min_samples_split,
            min_samples_leaf,
            max_features,
            max_depth,
            bootstrap)

            
for i in grid:
    print(str(i[0]) + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i[5]) + '\t' + str(i[6]) + '\t' + str(i[7]))
