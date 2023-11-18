import numpy as np
import networkx as nx
import scipy.stats as stat
import time, os, sys, random
from collections import defaultdict
home_dir = os.path.expanduser('~')
exec(open('%s/AMB_ICI_prediction/toolbox/load_hierarchy.py'%home_dir).read())

class Calculate_AMB_profiles:
    def __init__(self, input_df, hierarchy, transpose_input=False, random_shuffle=False, random_seed=42):
        if transpose_input == False:
            self.input_df = input_df
        else:
            tmp = defaultdict(list)
            samples = input_df[input_df.columns[0]].values
            genes = input_df.columns[1:]
            tmp['genes'] = genes
            for i in range(len(samples)):
                sample = samples[i]
                values = input_df.values[i][1:]
                tmp[sample] = values
            self.input_df = pd.DataFrame(data=tmp, columns=np.append(['genes'], samples))
            print('transpose')
        self.hierarchy = hierarchy
        ## load data
        if hierarchy.upper() == 'NEST':
            self.G, self.depth_df, self.geneList = load_NEST().hierarchy()
            self.term2gene = load_NEST().Term2Genes()
            if random_shuffle == True:
                random.seed(random_seed)
                tmp = []
                for term, geneList in zip(self.term2gene['name'].tolist(), self.term2gene['gene_id'].tolist()):
                    geneList = geneList.split(' ')
                    resampled_geneList = random.sample(self.geneList, len(geneList))
                    tmp.append(' '.join(map(str, resampled_geneList)))
                self.term2gene['gene_id'] = tmp

    def _sum(self, seq_df):
        ontotype = defaultdict(list)
        samples = seq_df.columns[1:]
        for i in range(self.term2gene.shape[0]):
            term, geneList = self.term2gene.iloc[i]['name'], self.term2gene.iloc[i]['gene_id'].split(' ')
            ontotype['name'].append(term)
            tmp = seq_df.loc[seq_df[seq_df.columns[0]].isin(geneList),:]
            for sample in samples:
                ontoScore = tmp[sample].sum()
                ontotype[sample].append(ontoScore)
        ontotype = pd.DataFrame(data=ontotype, columns=np.append(['name'], samples))
        return ontotype
    
    def Simple_Sum(self):
        return self._sum(seq_df=self.input_df)
