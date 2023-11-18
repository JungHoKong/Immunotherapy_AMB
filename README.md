# AMB
AMB is a visible machine learning framework that predicts patient's response to immune checkpoint inhibitors. This framework computes mutational burden at different levels of cancer protein assemblies and extract biomarkers of immunotherapy response.

# Publication
"Prediction of immunotherapy response using mutations to cancer protein assemblies", Kong et al.

# Requirements
This work was mainly performed using python 3.9.13. Key libraries that were used in this study is listed below:

- scikit-learn
- scikit-survival
- networkx
- lifelines
- plotly
- mageck

A full list of libraries used can be found in "py3.yml", "plotly_env.yml" and "mageckenv.yml".

# src
`Figure3_run_LOOCV.py` Perform Leave-One-Out Cross-Validation in Samstein cohort.
`Figure3_predict_validation_cohort.py` Make predictions in validation cohort.
`Figure6_run_MAGeCK.py` Run MAGeCK.

# toolbox
`calculate_AMB_profiles.py` Calculate AMB levels for each sample.
`load_hierarchy.py` Load data on protein assemblies.
`load_ICI.py` Load immunotherapy dataset.
