# AMB
Assembly-level Mutation Burden (AMB) model is a visible machine learning framework that predicts patient's response to immune checkpoint inhibitors. This framework computes mutational burden at different levels of cancer protein assemblies and extract biomarkers of immunotherapy response.

This code is used for the following work:
"Prediction of immunotherapy response using mutations to cancer protein assemblies", Kong et al.



# Requirements
This work was mainly performed using `python 3.9.13` and `Anaconda`. Key libraries that were used in this study are listed below:

- scikit-learn
- scikit-survival
- lifelines

Libraries above can be install using `pip install` or `conda install`. Installation of the libraries would take less than 5 minutes.

A full list of libraries used can be found in "py3.yml" in the 'env' folder.

Setting up the virtual environment can be done by running the following command line:

'conda env create -f env/py3.yml'



# Code description (within **src** folder)
1. *load_hierarchy.py* : loads assemblies and their associated genes from NeST hierarhcy (PMID: 34591613).
2. *calculate_AMB_profiles.py* : calculates AMB level per assembly. 
3. *train.py* : train on patient's survival data using AMB levels.



# Required inputs to train the AMB model (parameter names used to train the model are inside the parenthesis)
1. *Binarized mutation dataframe (input_df)* : A pandas DataFrame where the 1st column is gene IDs and the 2nd column and onwards are sample IDs.
2. *Assembly information (assembly_df)*: String ('NeST') or a pandas DataFrame. As a default, it will use NeST hierarchy. Custom pathway can be used to train an AMB model. When using custom dataframe, please see below:

'Dataframe containing assembly name and the genes within in the assembly. 
Use 'name' and 'gene_id' for columns indicating assembly name and genes, respectively.
In 'gene_id' columns, include genes within an assembly, separated by a space between gene names (e.g. "TP53 FGFR EGFR").'

3. *Maximum depth of each tree in the random survival forest model (max_depth)* : Integer value or *None* to use maximum depth (default=3). 


**train_AMB** will return (1) *a trained random survival forest model* and (2) *AMB scores (pandas dataframe) used to train the model*.
