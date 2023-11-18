## load immunotherapy treated patient data
import pandas as pd
from collections import defaultdict
import numpy as np
import os, time, sys
import scipy.stats as stat
home_dir = os.path.expanduser('~')



## load_data
def load_ICI_data(cohort, seq_type, response='binary_response'):
    '''
    cohort : 'Samstein', 'Hellmann'
    seq_type : 'mutation'
    '''
    if cohort == 'Samstein':
        seq_output = load_Samstein_seq(seq_type)
        ici_output = load_Samstein_ICI_response()
        ici_output = ici_output.rename(columns={'SAMPLE_ID':'patient'})
    elif cohort == 'Hellmann':
        seq_output = load_Hellmann_seq(seq_type)
        ici_output = load_Hellmann_ICI_response()
        ici_output = ici_output.rename(columns={'PATIENT_ID':'patient', 'response':'binary_response'})
    # sort values
    patients = list(set(seq_output.columns) & set(ici_output['patient']))
    seq_output = pd.DataFrame(data=seq_output, columns=np.append(['gene_id'], patients))
    ici_output = ici_output.loc[ici_output['patient'].isin(patients),:]
    ici_output = ici_output.set_index('patient')
    ici_output = ici_output.loc[patients]
    ici_output['patient'] = ici_output.index
    ici_output.reset_index(drop=True, inplace=True)
    return seq_output, ici_output



def load_ICI_TMB(cohort, binarize_TMB=True):
    '''
    column1: patient
    column2: TMB
    '''
    if cohort.lower() == 'samstein':
        df = pd.read_csv('%s/AMB_ICI_prediction/data/tmb_mskcc_2018/data_clinical_sample.txt'%home_dir, sep='\t', skiprows=4)
        df = df.rename(columns={'SAMPLE_ID':'patient', 'TMB_NONSYNONYMOUS':'TMB'})
        if binarize_TMB == True:
            df2 = pd.DataFrame()
            cancer_types = df['CANCER_TYPE'].unique()
            for cancer_type in cancer_types:
                tmp = df.loc[df['CANCER_TYPE']==cancer_type,:].sort_values(by='TMB', ascending=False)
                cutoff_index = int(np.floor(tmp.shape[0] * 0.2))
                tmpList = [1]*cutoff_index + [0]*(tmp.shape[0]-cutoff_index)
                tmp['high_TMB'] = tmpList
                df2 = pd.concat([df2, tmp])
            df = df2

    elif 'hellmann' in cohort.lower():
        df = pd.read_csv('%s/AMB_ICI_prediction/data/Hellmann_etal/nsclc_mskcc_2018/data_clinical_sample.txt'%(home_dir), sep='\t', header=4)
        df = df.dropna(subset=['TMB_NONSYNONYMOUS'])
        df.reset_index(drop=True, inplace=True)
        df = df.rename(columns={'PATIENT_ID':'patient', 'TMB_NONSYNONYMOUS':'TMB'})
        if binarize_TMB == True:
            df = df.sort_values(by='TMB', ascending=False)
            cutoff_index = int(np.floor(df.shape[0] * 0.2))
            tmpList = [1] * cutoff_index + [0] * (df.shape[0]-cutoff_index)
            df['high_TMB'] = tmpList
    return df





##----------
# load mskcc seq 
def load_Samstein_seq(seq_type):
    output = defaultdict(list)
    if seq_type.lower() == 'mutation':
        fi_dir = '%s/AMB_ICI_prediction/data/Samstein'%home_dir
        if not 'data_mutations_jkong.txt' in os.listdir(fi_dir):
            # load data
            df = pd.read_csv('%s/data_mutations.txt'%fi_dir, sep='\t', low_memory=False)
            mut_types = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'Splice_Site', 'Splice_Region', 'In_Frame_Ins', 'In_Frame_Del']
            df = df.loc[df['Variant_Classification'].isin(mut_types),:]
            # genes
            genes = sorted(df['Hugo_Symbol'].unique())
            output['gene_id'] = genes
            patients = df['Tumor_Sample_Barcode'].unique()                
            for patient in patients: 
                tmp = []
                for gene in genes:
                    if gene in df.loc[df['Tumor_Sample_Barcode']==patient,:]['Hugo_Symbol'].tolist():
                        tmp.append(1)
                    else:
                        tmp.append(0)
                output[patient] = tmp
            output = pd.DataFrame(data=output, columns=np.append(['gene_id'], patients))
            output.to_csv('%s/data_mutations_jkong.txt'%fi_dir, sep='\t', index=False)
        else:
            output = pd.read_csv('%s/data_mutations_jkong.txt'%fi_dir, sep='\t')
    return output


def load_Samstein_ICI_response():
    output = defaultdict(list)
    fi_dir = '%s/AMB_ICI_prediction/data/Samstein'%home_dir
    df1 = pd.read_csv('%s/data_clinical_sample.txt'%fi_dir, sep='\t', header=4)
    df2 = pd.read_csv('%s/data_clinical_patient.txt'%fi_dir, sep='\t', header=4)
    df = pd.merge(df1, df2, on='PATIENT_ID', how='inner')

    return df


##----------
# load Hellmann (lsclc mskcc) 
def load_Hellmann_seq(seq_type):
    output = defaultdict(list)
    if seq_type.lower() == 'mutation':
        fi_dir = '%s/AMB_ICI_prediction/data/Hellmann_etal/nsclc_mskcc_2018'%home_dir
        if not 'data_mutations_jkong.txt' in os.listdir(fi_dir):
            # load data
            df = pd.read_csv('%s/data_mutations_extended.txt'%fi_dir, sep='\t', low_memory=False)
            mut_types = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'Splice_Site', 'Splice_Region', 'In_Frame_Ins', 'In_Frame_Del']
            df = df.loc[df['Variant_Classification'].isin(mut_types),:]
            # sample ID 
            pdf = pd.read_csv('%s/data_clinical_sample.txt'%fi_dir, sep='\t', header=4)
            col_name = {pdf['SAMPLE_ID'].tolist()[i] : pdf['PATIENT_ID'].tolist()[i] for i in range(pdf.shape[0])}

            # genes
            genes = sorted(df['Hugo_Symbol'].unique())
            output['gene_id'] = genes
            patients = df['Tumor_Sample_Barcode'].unique()                
            for patient in patients: 
                tmp = []
                for gene in genes:
                    if gene in df.loc[df['Tumor_Sample_Barcode']==patient,:]['Hugo_Symbol'].tolist():
                        tmp.append(1)
                    else:
                        tmp.append(0)
                output[patient] = tmp
            output = pd.DataFrame(data=output, columns=np.append(['gene_id'], patients))
            output = output.rename(columns=col_name)
            output.to_csv('%s/data_mutations_jkong.txt'%fi_dir, sep='\t', index=False)
        else:
            output = pd.read_csv('%s/data_mutations_jkong.txt'%fi_dir, sep='\t')
    return output


def load_Hellmann_ICI_response():
    output = defaultdict(list)
    fi_dir = '%s/AMB_ICI_prediction/data/Hellmann_etal/nsclc_mskcc_2018'%home_dir
    df = pd.read_csv('%s/data_clinical_patient.txt'%fi_dir, sep='\t', header=4)
    binary_response = []
    for response in df['BEST_OVERALL_RESPONSE'].tolist():
        if response in ['CR', 'PR']:
            binary_response.append(1)
        elif response in ['SD', 'PD']:
            binary_response.append(0)
        else:
            binary_response.append('na')
    df['response'] = binary_response
    df = df.loc[df['response']!='na',:].reset_index(drop=True)
    return df


##---------
# load datasets
def load_datasets(datasets):
    data_df = {}
    for dataset in datasets:
        mdf, pdf = load_ICI_data(dataset, 'mutation')
        mdf = mdf.rename(columns={'gene_id':'genes'})
        if dataset == 'Samstein':
            mdf = pd.DataFrame(data=mdf, columns=np.append(['genes'], list(set(mdf.columns) & set(pdf.loc[pdf['CANCER_TYPE'].isin(['Non-Small Cell Lung Cancer', 'Bladder Cancer']),:]['patient'].tolist()))))
            OS_months = [pdf.loc[pdf['patient']==sample,:]['OS_MONTHS'].tolist()[0] for sample in mdf.columns[1:]]
            OS_status = [pdf.loc[pdf['patient']==sample,:]['OS_STATUS'].tolist()[0] for sample in mdf.columns[1:]]
            OS_status = [1 if status=='1:DECEASED' else 0 for status in OS_status]
        data_df[dataset] = mdf 

    # match genes
    if (len(datasets) > 1) and ('Samstein' in datasets):
        genes = data_df['Samstein']['genes'].tolist()
        for dataset in datasets:
            if dataset == 'Samstein': continue
        data_df[dataset] = data_df[dataset].loc[data_df[dataset]['genes'].isin(genes),:]
        tmp_genes = [gene for gene in genes if not gene in data_df[dataset]['genes'].tolist()]
        tmp_df = defaultdict(list)
        for gene in tmp_genes:
            tmp_df['genes'].append(gene)
            for col in data_df[dataset].columns:
                if not col == 'genes':
                    tmp_df[col].append(0)
        tmp_df = pd.DataFrame(tmp_df)
        data_df[dataset] = data_df[dataset].append(tmp_df)
    return data_df
