#!/usr/bin/env python3

import sys
import pandas as pd
import argparse

## Read command line args

parser = argparse.ArgumentParser(description='Process input and output files.')
parser.add_argument('--joint', help='Joint table path')
parser.add_argument('--patient', help='Patient identifier')
args = parser.parse_args()

#Import input data
joint_table = args.joint
input_data = pd.read_csv(joint_table, sep = '\t')
patient_id = args.patient


#input data preprocessing
def input_processing(data, patient_id):
    new_columns = {'normal_cn': 2}
    df['mutation_id'] = df.apply(lambda row: f"{patient_id}:{row['chr'][3:]}:{row['from']}:{row['alt']}", axis=1)
    df = data.assign(**new_columns)
    df[['major_cn', 'minor_cn']] = df['karyotype'].str.split(':', expand=True)
    df = df.drop(columns=['karyotype'])
    # df['normal_cn'] = df['normal_cn'].astype(int)
    # df['major_cn'] = df['major_cn'].astype(int)
    # df['minor_cn'] = df['minor_cn'].astype(int)
    df = df.rename(columns={'Indiv': 'sample_id', 'NR': 'ref_counts', 'NV': 'alt_counts', 'purity': 'tumour_content'})
    df = df.loc[:, ['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn','minor_cn','tumour_content']]
    #df = df.style.hide(axis='index') #omit row indexes
    return df