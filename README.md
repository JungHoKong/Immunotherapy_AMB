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


## Setting up the virtual environment
Setting up the virtual environment can be done by running the following command line:

`conda env create -f env/py3.yml`



# Code description (within **src** folder)
1. *load_hierarchy.py* : loads assemblies and their associated genes from NeST hierarhcy (PMID: 34591613).
2. *calculate_AMB_profiles.py* : calculates AMB level per assembly. 
3. *train.py* : train on patient's survival data using AMB levels.



# Required inputs to train the AMB model (*train.py*)
Parameter names used in the *train.py* code are inside the parenthesis
1. *Binarized mutation dataframe (input_df)* : A pandas DataFrame where the 1st column is gene IDs. The 2nd column and onwards should be sample IDs.
2. *Survival data (survival_df)* : A pandas DataFrame that includes patient's overall or progression-free survival. The columns of the dataframe should include **sample**, **month** and **status**.  
3. *Assembly information (assembly_df)*: String ('NeST') or a pandas DataFrame. As a default, it will use NeST hierarchy. Custom pathway can be used to train an AMB model. When using custom dataframe, use 'name' and 'gene_id' for columns indicating assembly name and genes, respectively.
In 'gene_id' columns, include genes within an assembly and separate gene names by spaces (e.g. "TP53 FGFR EGFR").
4. *Maximum depth of each tree in the random survival forest model (max_depth)* : Integer value (default=3) or *None* (maximum depth). 


**train_AMB** will return the following:
1. A trained random survival forest model
2. AMB scores (pandas dataframe) used to train the model.

