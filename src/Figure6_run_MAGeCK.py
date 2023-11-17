import pandas as pd
import os, time, sys

home_dir = os.path.expanduser('~')
exec(open('../toolbox/load_ICI.py').read())

# run MAGeCK
output_dir = '~/TMB_profile/data/in_vivo_mice_CRISPR_screening'


# result 
tmp = pd.read_csv('../data/in_vivo_mice_CRISPR_screening/raw_counts.txt', sep='\t')
control_sgrnas = tmp.loc[tmp['Gene']=='Control',:]['sgRNA'].unique().tolist()

fo = open('../data/in_vivo_mice_CRISPR_screening/control_sgrna.txt', 'w')
for control_sgrna in control_sgrnas:
    fo.write('%s\n'%control_sgrna)
fo.close()

os.chdir(output_dir)
os.system('mageck test -k %s/raw_counts.txt -c BW-I1,BW-I2,BW-I3,BW-I4,BW-I5,BW-I6,BW-I7,BW-I8,BW-I9,BW-I10,BW-I11,BW-I12 -t BW-P1,BW-P2,BW-P3,BW-P4,BW-P5,BW-P6,BW-P7,BW-P8,BW-P9,BW-P10,BW-P11,BW-P12 -n BW_I_vs_BW_P_median --norm-method total --pdf-report --control-sgrna control_sgrna.txt --normcounts-to-file --gene-lfc-method secondbest'%(output_dir))

