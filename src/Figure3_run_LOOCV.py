import pandas as pd
import lifelines 
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index 
from sksurv.ensemble import RandomSurvivalForest
from sklearn.inspection import permutation_importance
from collections import defaultdict
import numpy as np
import scipy.stats as stat
import os, time, sys, warnings, random
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
warnings.simplefilter(action='ignore', category=FutureWarning)
home_dir = os.path.expanduser('~')
exec(open('%s/TMB_profile/toolbox/load_hierarchy.py'%home_dir).read())
exec(open('%s/TMB_profile/toolbox/load_ICI.py'%home_dir).read())
exec(open('%s/TMB_profile/toolbox/calculate_TMB_profiles.py'%home_dir).read())
exec(open('%s/TMB_profile/toolbox/measure_performance.py'%home_dir).read())

## initalize
hierarchy = 'NEST' 
train_dataset = 'Samstein'
test_types = ['TMB', 'TMB_profile', 'gene']
cox_penalizer = np.arange(0.1, 1.1, 0.1)
max_depths = [3,5,None]


## load datasets
data_df = load_datasets([train_dataset])


## Calculate TMB profiles
TMB_profiles = {}
tmp = Calculate_TMB_profiles(input_df=data_df[train_dataset], hierarchy=hierarchy).Simple_Sum()
TMB_profiles[train_dataset] = tmp.sort_values(by='name')
del tmp


# input dataframe
X_train = TMB_profiles[train_dataset]
X_train_genes = data_df[train_dataset]


##==========================================================================
## LOOCV predictions
# input features
train_input = defaultdict(list)
X_train = X_train.sort_values(by='name')
features = X_train['name'].tolist()
for idx, feature in enumerate(features):
    train_input[feature] = X_train.T.values[1:][:,idx]
samples = X_train.columns[1:]

# response label
mdf, pdf = load_ICI_data(train_dataset, 'mutation')
OS_months = [pdf.loc[pdf['patient']==sample,:]['OS_MONTHS'].tolist()[0] for sample in samples]
OS_status = [pdf.loc[pdf['patient']==sample,:]['OS_STATUS'].tolist()[0] for sample in samples]
OS_status = [1 if status=='1:DECEASED' else 0 for status in OS_status]
train_input['month'] = OS_months; train_input['status'] = OS_status
train_input = pd.DataFrame(train_input)
for col in train_input.columns:
    if train_input[col].sum() == 0:
        train_input = train_input.drop(col, axis=1)


# predict OS
for test_idx in range(train_input.shape[0]):
    train_idx = [idx for idx in range(train_input.shape[0]) if idx != test_idx]
    X_train_input, X_test_input = train_input.iloc[train_idx], train_input.iloc[[test_idx]]

    # prepare data
    X_train_ = X_train_input.copy(); X_train_ = X_train_.drop(['month', 'status'], axis=1)
    X_test_ = X_test_input.copy(); X_test_ = X_test_.drop(['month', 'status'], axis=1)
    y_train, y_test = [], []
    for i in range(X_train_.shape[0]):
        a = tuple([X_train_input.iloc[i]['status'].astype(np.bool_), X_train_input.iloc[i]['month']])
        y_train.append(a)
    for i in range(X_test_.shape[0]):
        a = tuple([X_test_input.iloc[i]['status'].astype(np.bool_), X_test_input.iloc[i]['month']])
        y_test.append(a)

    y_train = np.array(y_train, dtype=[('status','?'), ('month','<f8')])
    y_test = np.array(y_test, dtype=[('status','?'), ('month','<f8')])

    # nested CV
    best_CV_score, best_param, CV_scores = OS_select_param('RandomSurvivalForest', max_depths, X_train_input, 'month', 'status')
    rsf = RandomSurvivalForest(n_estimators=500, max_depth=best_param)
    rsf.fit(X_train_, y_train)
    pred = rsf.predict(X_test_)[0]


