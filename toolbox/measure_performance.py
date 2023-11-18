from sklearn.metrics import *
from sklearn.model_selection import *
from sklearn.preprocessing import *
import scipy.stats as stat
import numpy as np
from sklearn.linear_model import *
from sklearn.mixture import *
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.svm import SVR, SVC
from sklearn.ensemble import *
from sklearn.model_selection import GridSearchCV
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sksurv.ensemble import RandomSurvivalForest


## function
def test_performance(output, y_test, pred, test_dataset, test_type, survival_model, survival_model_details):
    fpr, tpr = [], []
    fpr, tpr, _ = roc_curve(y_test, pred, pos_label=0)
    AUC = auc(fpr, tpr)
    print('%s, AUC = %.4f'%(test_type, AUC))
    for key, value in zip(['test_type', 'survival_model', 'details', 'test_dataset', 'metric', 'score', 'fpr', 'tpr'], [test_type, survival_model, survival_model_details, test_dataset, 'AUC', AUC, ';'.join(map(str,fpr)), ';'.join(map(str,tpr))]):
        output[key].append(value)
    return output



def OS_select_param(model, parameters, X, duration_col, event_col, nFold=5, max_features='sqrt'):
    kf = KFold(n_splits=nFold)
    best_score = [0, None]
    CV_scores = defaultdict(list)
    try:
        for param in parameters:
            validation_scores = []
            # parameter sweep
            if model.upper() == 'COX':
                for train, test in kf.split(X):
                    X_train, X_test = X.iloc[train], X.iloc[test]
                    cph = CoxPHFitter(penalizer=param)
                    cph.fit(X_train, duration_col=duration_col, event_col=event_col)
                    score = cph.score(X_test, scoring_method='concordance_index')
                    validation_scores.append(score)
            if model == 'RandomSurvivalForest':
                X_ = X.copy()
                X_ = X_.drop([duration_col, event_col], axis=1)
                y = []
                for i in range(X.shape[0]):
                    a = tuple([X.iloc[i][event_col].astype(np.bool_), X.iloc[i][duration_col]])
                    y.append(a)
                y = np.array(y, dtype=[('status','?'), ('month','<f8')])

                for train, test in kf.split(X_):
                    X_train, X_test = X_.iloc[train], X_.iloc[test]
                    y_train, y_test = y[train], y[test]
                    rsh = RandomSurvivalForest(n_estimators=100, max_depth=param, max_features=max_features)
                    rsh.fit(X_train, y_train)
                    score = rsh.score(X_test, y_test)
                    validation_scores.append(score)
            # evaluate
            CV_scores[param] = validation_scores
            if np.mean(validation_scores) >= best_score[0]:
                best_score = [np.mean(validation_scores), param]
    except: pass
    return best_score[0], best_score[1], CV_scores
