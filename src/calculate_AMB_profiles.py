import numpy as np
import networkx as nx
import scipy.stats as stat
import time, os, sys, random
from collections import defaultdict
sys.path.append('./')
from load_hierarchy import *


class Calculate_AMB_profiles:
    def __init__(self, input_df, assembly_df, random_shuffle=False, random_seed=42):
        '''
        Calculates AMB levels.
        
        Inputs
        --------------------------------
        input_df : dataframe
                Binarized mutations. 
                Gene names in the first column. 
                Sample IDs starting from the second column.
                
                
        assembly_df : str or dataframe
                If assembly_df = 'NeST', use NeST hierarchy (PMID: 34591613).
                
                When using custom dataframe, please see below:
                Dataframe containing assembly name and the genes contained in the assembly. 
                Use 'name' and 'gene_id' for columns indicating assembly name and genes, respectively.
                In 'gene_id' columns, separate genes with a space. 
        
        
        random_shuffle : bool
                Randomly shuffle relationship between an assembly and its corresponding genes.
        
        
        random_seed : int
                Random number generator to reproduce results. Ignored if 'random_shuffle=False'.
        '''
        self.input_df = input_df
        
        ## assembly info
        if type(assembly_df) == str:
            if assembly_df.upper() == 'NEST':
                self.term2gene = load_NEST().Term2Genes()
        else:
            self.term2gene = assembly_df
            
            
        ## shuffle assembly-gene relationship
        if random_shuffle == True:
            random.seed(random_seed)
            tmp = []
            for term, geneList in zip(self.term2gene['name'].tolist(), self.term2gene['gene_id'].tolist()):
                geneList = geneList.split(' ')
                resampled_geneList = random.sample(self.geneList, len(geneList))
                tmp.append(' '.join(map(str, resampled_geneList)))
            self.term2gene['gene_id'] = tmp

            
            
    def _sum(self, seq_df):
        out = defaultdict(list)
        samples = seq_df.columns[1:]
        for i in range(self.term2gene.shape[0]):
            term, geneList = self.term2gene.iloc[i]['name'], self.term2gene.iloc[i]['gene_id'].split(' ')
            out['name'].append(term)
            tmp = seq_df.loc[seq_df[seq_df.columns[0]].isin(geneList),:]
            for sample in samples:
                Score = tmp[sample].sum()
                out[sample].append(Score)
        out = pd.DataFrame(data=out, columns=np.append(['name'], samples))
        return out
    
    
    
    def AMB_score(self):
        return self._sum(seq_df=self.input_df)
