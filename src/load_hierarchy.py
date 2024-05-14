## load hierarchies
import pandas as pd
from collections import defaultdict
import os, time, sys
import networkx as nx
import numpy as np
import scipy.stats as stat
home_path = os.path.expanduser('~')



## main functions
def load_hierarchy_genes(hierarchy_name='NEST'):
    if hierarchy_name.upper() == 'NEST':
        df = load_NEST().Term2Genes()
    return df


## NeST
class load_NEST:
    # init
    def __init__(self):
        self.fi_dir = '../data/NeST'
    
    # load genes
    def Term2Genes(self):
        df = pd.read_csv('%s/NeST_node.csv'%self.fi_dir)
        df = pd.DataFrame(data=df, columns=['name', 'Genes', 'Annotation']); df = df.rename(columns={'Genes':'gene_id'})
        return df


