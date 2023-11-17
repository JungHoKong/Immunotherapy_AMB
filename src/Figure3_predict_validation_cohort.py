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
test_dataset = 'Hellmann'
survival_models = ['Cox', 'RandomSurvivalForest']
test_types = ['TMB', 'TMB_profile', 'gene']
cox_penalizer = np.arange(0.1, 1.1, 0.1)
max_depths = [3,5,None]


## load datasets
datasets = [train_dataset, test_dataset]
data_df = load_datasets(datasets)


## Calculate TMB profiles
TMB_profiles = {}
for dataset in datasets:
    tmp = Calculate_TMB_profiles(input_df=data_df[dataset], hierarchy=hierarchy).Simple_Sum()
    TMB_profiles[dataset] = tmp.sort_values(by='name')
    del tmp


# input dataframe
X_train = TMB_profiles[train_dataset]
X_train_genes = data_df[train_dataset]

# TMB
X_train_TMB = defaultdict(list)
TMB_df = {}
for dataset in datasets:
    TMB = load_ICI_TMB(dataset)
    common_samples = sorted(list(set(TMB['patient'].tolist()) & set(TMB_profiles[dataset].columns)))
    TMB = TMB.loc[TMB['patient'].isin(common_samples),:].sort_values(by='patient')['TMB'].tolist()
    # TMB_df
    tmp = {common_samples[i]:TMB[i] for i in range(len(TMB))}
    TMB_df[dataset] = tmp; del tmp
for key, value in TMB_df[train_dataset].items():
    X_train_TMB[key] = [value]
X_train_TMB = pd.DataFrame(data=X_train_TMB, columns=X_train.columns[1:])



##=======================================================
## predictions
output = defaultdict(list)
output2 = defaultdict(list)

# test input dataframe
X_test = TMB_profiles[test_dataset]
X_test_genes = data_df[test_dataset]
X_test_TMB = defaultdict(list)
for key, value in TMB_df[test_dataset].items():
    X_test_TMB[key] = [value]
X_test_TMB = pd.DataFrame(data=X_test_TMB, columns=X_test.columns[1:])


# test types
for test_type in test_types:
    train_input = defaultdict(list)
    test_input = defaultdict(list)

    if 'TMB_profile' == test_type:
        X_train = X_train.sort_values(by='name'); X_test = X_test.sort_values(by='name')
        features = X_train['name'].tolist()
        for idx, feature in enumerate(features):
            train_input[feature] = X_train.T.values[1:][:,idx]
            test_input[feature] = X_test.T.values[1:][:,idx]
    if 'TMB' == test_type:
        train_input['TMB'] = X_train_TMB.values.ravel()
        test_input['TMB'] = X_test_TMB.values.ravel()
    if 'gene' == test_type:
        X_train_genes = X_train_genes.sort_values(by='genes'); X_test_genes = X_test_genes.sort_values(by='genes')
        geneList = X_train_genes['genes'].tolist()
        for g_idx, gene in enumerate(geneList):
            train_input[gene] = X_train_genes.T.values[1:][:,g_idx]
            test_input[gene] = X_test_genes.T.values[1:][:,g_idx]
            
    # training response label
    mdf, pdf = load_ICI_data(train_dataset, 'mutation')
    OS_months = [pdf.loc[pdf['patient']==sample,:]['OS_MONTHS'].tolist()[0] for sample in X_train.columns[1:]]
    OS_status = [pdf.loc[pdf['patient']==sample,:]['OS_STATUS'].tolist()[0] for sample in X_train.columns[1:]]
    OS_status = [1 if status=='1:DECEASED' else 0 for status in OS_status]
    train_input['month'] = OS_months; train_input['status'] = OS_status
    train_input = pd.DataFrame(train_input)
    test_input = pd.DataFrame(test_input)
    for col in train_input.columns:
        if train_input[col].sum() == 0:
            train_input = train_input.drop(col, axis=1)
            test_input = test_input.drop(col, axis=1)

    # test response label
    mdf, pdf = load_ICI_data(test_dataset, 'mutation')
    y_test = []
    for sample in X_test.columns[1:]:
        y_test.append(pdf.loc[pdf['patient']==sample,:]['binary_response'].tolist()[0])
    

    # predict response
    for survival_model in survival_models:
        train_X, test_X = train_input, test_input

        if survival_model == 'Cox':
            if ('TMB_profile' == test_type) or ('gene' == test_type): continue
            # nested CV
            best_CV_score, best_param, CV_scores = OS_select_param('COX', cox_penalizer, train_X, 'month', 'status')
            cph = CoxPHFitter(penalizer=best_param)
            cph.fit(train_X, duration_col = 'month', event_col = 'status')
            pred = cph.predict_partial_hazard(test_X)
            output = test_performance(output, y_test, pred, test_dataset, test_type, survival_model, survival_model_details='cox_penalizer_%s'%best_param)


        elif survival_model == 'RandomSurvivalForest':
            if 'TMB' == test_type: continue
            # prepare data
            X_train_ = train_X.copy()
            X_train_ = X_train_.drop(['month', 'status'], axis=1)
            X_test_ = test_X.copy()
            y_train = []
            for i in range(X_train_.shape[0]):
                a = tuple([train_X.iloc[i]['status'].astype(np.bool_), train_X.iloc[i]['month']])
                y_train.append(a)
            # nested CV
            best_CV_score, best_param, CV_scores = OS_select_param('RandomSurvivalForest', max_depths, train_X, 'month', 'status')
            rsf = RandomSurvivalForest(n_estimators=500, max_depth=best_param)
            y_train = np.array(y_train, dtype=[('status','?'), ('month','<f8')])
            rsf.fit(X_train_, y_train)
            pred = rsf.predict(X_test_)
            output = test_performance(output, y_test, pred, test_dataset, test_type, survival_model, survival_model_details='max_depth_%s'%best_param)
            
        # output2
        for sample, obs, prd in zip(X_test.columns[1:], y_test, pred):
            for key, value in zip(['train_dataset', 'test_dataset', 'sample', 'test_type', 'details', 'obs', 'pred'], [train_dataset, test_dataset, sample, test_type, survival_model, obs, prd]):
                output2[key].append(value)
output = pd.DataFrame(output)
output2 = pd.DataFrame(output2)
