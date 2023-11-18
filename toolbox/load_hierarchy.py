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
    elif hierarchy_name.upper() == 'REACTOME':
        df = load_reactome().Term2Genes()
    return df



###############################################################################################
## NeST
class load_NEST:
    # init
    def __init__(self):
        self.fi_dir = '%s/AMB_ICI_prediction/data/NeST'%home_path
    
    
    # load genes
    def Term2Genes(self):
        df = pd.read_csv('%s/NeST_node.csv'%self.fi_dir)
        df = pd.DataFrame(data=df, columns=['name', 'Genes', 'Annotation']); df = df.rename(columns={'Genes':'gene_id'})
        return df


    ## load hierarchies
    def hierarchy(self):
        # init
        geneList = []
        output = defaultdict(list)
        term2gene = self.Term2Genes()
        
        # load data
        df = pd.read_csv('%s/NeST_edge.sif'%self.fi_dir, sep='\t', header=None)
        df2 = pd.read_csv('%s/NeST_node.csv'%self.fi_dir)
        # directed graph
        G = nx.DiGraph()
        # add edges between terms
        for i in range(df.shape[0]):
            tmp = df.iloc[i].values
            T1, T2 = tmp[0], tmp[2]
            G.add_edge(T1, T2)
        # add edges between terms and genes
        for term, genes in zip(df2['name'].tolist(), df2['Genes'].tolist()):
            genes = genes.split(' ')
            geneList = list(set(geneList).union(set(genes)))
            for gene in genes:
                G.add_edge(term, gene)
        # calculate shortest path from a term to genes (leaves)
        if not 'NEST_max_depth.txt' in os.listdir(self.fi_dir):
            print('calculating shortest path from a term to genes, %s'%time.ctime())
            depth_dic = defaultdict(list)
            for term_i, term in enumerate(df2['name'].tolist()):
                term_i += 1
                term_genes = term2gene.loc[term2gene['name']==term,:]['gene_id'].tolist()[0].split(' ')
                if (term_i == 1) or (term_i % 2 == 0):
                    print('     running %s (%s), %s / %s, %s'%(term, len(term_genes), term_i, df2.shape[0], time.ctime()))
                tmp_depths = []
                for gene in list(set(geneList) & set(term_genes)): # calculate distances only for the genes included in the term
                    try:
                        for path in nx.all_simple_paths(G, term, gene):
                            tmp_depths.append(len(path))
                    except: pass
                depth = int(np.max(tmp_depths))
                depth_dic[depth].append(term) 
                # output
                output['term'].append(term)
                output['depth'].append(depth)
            output = pd.DataFrame(output)
            output.to_csv('%s/NEST_max_depth.txt'%self.fi_dir, sep='\t', index=False)
            print('calculation finished, %s'%time.ctime())
        else:
            output = pd.read_csv('%s/NEST_max_depth.txt'%self.fi_dir, sep='\t')
        return G, output, geneList



###############################################################################################

## reactome
class load_reactome:
    def __init__(self):
        self.fi_dir = '%s/AMB_ICI_prediction/data/reactome'%home_path
    
    # load genes
    def Term2Genes(self):
        if not 'reactome_node.txt' in os.listdir(self.fi_dir):
            output = defaultdict(list)
            G = nx.DiGraph()
            df = pd.read_csv('%s/reactome.human.ddot'%self.fi_dir, header=None, sep='\t')
            geneList = sorted(df.loc[df[2]=='gene',:][1].unique().tolist())
            for vList in df.values:
                t1, t2 = vList[0], vList[1]
                G.add_edge(t1, t2)

            # term genes
            for node in G.nodes():
                if not node in geneList:
                    tmp = []
                    for gene in geneList:
                        try:
                            shortest_paths = list(nx.shortest_simple_paths(G, node, gene))
                            tmp.append(gene)
                        except: pass
                    output['name'].append(node)
                    output['Genes'].append(' '.join(map(str, tmp)))
            output = pd.DataFrame(output)
            output.to_csv('%s/reactome_node.txt'%(self.fi_dir), sep='\t', index=False)
            df = output
        else:
            df = pd.read_csv('%s/reactome_node.txt'%self.fi_dir, sep='\t')
        df = pd.DataFrame(data=df, columns=['name', 'Genes']); df = df.rename(columns={'Genes':'gene_id'})
        df = df.dropna()
        return df
    

    ## load hierarchies
    def hierarchy(self):
        # init
        geneList = []
        output = defaultdict(list)
        term2gene = self.Term2Genes()
        
        # load data
        df = pd.read_csv('%s/reactome.human.ddot'%self.fi_dir, header=None, sep='\t')
        # directed graph
        G = nx.DiGraph()
        # add edges between terms
        for i in range(df.shape[0]):
            tmp = df.iloc[i].values
            if tmp[2] != 'gene':
                T1, T2 = tmp[0], tmp[1]
                G.add_edge(T1, T2)
        # add edges between terms and genes
        for term, genes in zip(term2gene['name'].tolist(), term2gene['gene_id'].tolist()):
            genes = genes.split(' ')
            geneList = list(set(geneList).union(set(genes)))
            for gene in genes:
                G.add_edge(term, gene)
        # calculate shortest path from a term to genes (leaves)
        fiName = 'reactome_max_depth.txt'
        if not fiName in os.listdir(self.fi_dir):
            print('calculating shortest path from a term to genes, %s'%time.ctime())
            depth_dic = defaultdict(list)
            for term_i, term in enumerate(term2gene['name'].tolist()):
                term_i += 1
                term_genes = term2gene.loc[term2gene['name']==term,:]['gene_id'].tolist()[0].split(' ')
                if (term_i == 1) or (term_i == 2) or (term_i % 10 == 0):
                    print('     running %s (%s), %s / %s, %s'%(term, len(term_genes), term_i, term2gene.shape[0], time.ctime()))
                tmp_depths = []
                for gene in list(set(geneList) & set(term_genes)): # calculate distances only for the genes included in the term
                    try:
                        for path in nx.all_simple_paths(G, term, gene):
                            tmp_depths.append(len(path))
                    except: pass
                depth = int(np.max(tmp_depths))
                depth_dic[depth].append(term) 
                # output
                output['term'].append(term)
                output['depth'].append(depth)
            output = pd.DataFrame(output)
            output.to_csv('%s/%s'%(self.fi_dir, fiName), sep='\t', index=False)
            print('calculation finished, %s'%time.ctime())
        else:
            output = pd.read_csv('%s/%s'%(self.fi_dir, fiName), sep='\t')
        return G, output, geneList


