import pandas as pd
import pandas
import lifelines 
from sksurv.ensemble import RandomSurvivalForest
from collections import defaultdict
import numpy as np
import scipy.stats as stat
import os, time, sys, warnings, random
warnings.simplefilter(action='ignore', category=FutureWarning)
sys.path.append('./')
from load_hierarchy import *
from calculate_AMB_profiles import *

# AMB
def train_AMB(X_train, survival_df, assembly_df='NeST', max_depth=3, random_shuffle=False, random_seed=42):
    '''
    Train a random survival forest model using assembly-level mutation burden.


    --------------------------------
    Inputs
    --------------------------------
    X_train : dataframe
            Binarized mutations. 
            Gene names in the first column. 
            Sample IDs starting from the second column.
    

    survival_df : dataframe
            Patient's overall / progression-free survival. 
            Use 'patient', 'months' and 'status' for columns indicating patient ID, survival months and status, respectively.

            
    assembly_df : str or dataframe, default='NeST'
            If assembly_df = 'NeST', use NeST hierarchy (PMID: 34591613).
            
            When using custom dataframe, please see below:
            Dataframe containing assembly name and the genes within in the assembly. 
            Use 'name' and 'gene_id' for columns indicating assembly name and genes, respectively.
            In 'gene_id' columns, separate genes with a space. 


    max_depth : int, default=3
            max depth hyperparameter to use to train a random survival forest.

    
    random_shuffle : bool, default=False
            Randomly shuffle relationship between an assembly and its corresponding genes.
    
    
    random_seed : int, default=42
            Random number generator to reproduce results. Ignored if 'random_shuffle=False'.


    --------------------------------
    Outputs
    --------------------------------
    rsf : Returns a trained model
    X_train_ : input dataframe, containing AMB scores, that is used to train the model
    '''

    # Check input types
    assert isinstance(X_train, pandas.core.frame.DataFrame), 'Expected X_train to be a pandas.dataframe'
    assert isinstance(survival_df, pandas.core.frame.DataFrame), 'Expected survival_df to be a pandas.dataframe'


    # Calculate AMB proflies
    AMB_profiles = Calculate_AMB_profiles(X_train, assembly_df, random_shuffle=random_shuffle, random_seed=random_seed).AMB_score()


    # input features
    train_input = defaultdict(list)
    samples = AMB_profiles.columns[1:]

    # Create an input for the random forest model
    train_input = defaultdict(list)
    months = [survival_df.loc[survival_df['patient']==sample,:]['months'].astype(float).tolist()[0] for sample in samples]
    status = [survival_df.loc[survival_df['patient']==sample,:]['status'].astype(float).tolist()[0] for sample in samples]

    features = AMB_profiles[AMB_profiles.columns[0]].tolist()
    for idx, feature in enumerate(features):
        train_input[feature] = AMB_profiles.T.values[1:][:,idx]
    train_input['months'] = months
    train_input['status'] = status
    
    train_input = pd.DataFrame(train_input)
    for col in train_input.columns:
        if train_input[col].sum() == 0:
            train_input = train_input.drop(col, axis=1)


    # X_train_
    X_train_ = train_input.copy()
    X_train_ = X_train_.drop(['months', 'status'], axis=1)

    # y_train
    y_train = []
    for i in range(X_train_.shape[0]):
        a = tuple([train_input.iloc[i]['status'].astype(np.bool_), train_input.iloc[i]['months']])
        y_train.append(a)
    y_train = np.array(y_train, dtype=[('status','?'), ('months','<f8')])


    # predict
    rsf = RandomSurvivalForest(n_estimators=500, max_depth=max_depth)
    rsf.fit(X_train_, y_train)
    
    # out
    return rsf, X_train_
