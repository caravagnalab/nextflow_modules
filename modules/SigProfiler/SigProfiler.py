#!/usr/bin/env python3

import pandas as pd
import os
import shutil
from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

#create directory
#path = '/orfeo/scratch/cdslab/kdavydzenka/SIGPROFILER/'
#os.mkdir(path)

#import input data
input_data = pd.read_csv("mut_joint_table.txt", sep = '\t')

#input data preprocessing
def input_processing(data):
    new_columns = {'Project': "CLL", 'Genome': 'GRCh37', 'Type': "SOMATIC", 'mut_type': "SNP"}
    df = data.assign(**new_columns)
    df['chr'] = df['chr'].astype(str).str[3:]
    df = df.rename(columns={'sample_id': 'Sample', 'patient_id': 'ID', 'chr': 'chrom', 'from': 'pos_start', 'from': 'pos_start',
                       'to': 'pos_end'})
    df = df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
    #df = df.style.hide(axis='index') #row merge

    return df

input_data = input_processing(input_data)

#saving input matrix to txt
input_data.to_csv('CLL/input_data.txt', sep='\t', index=False, header=True)

#mutations matrix generation
input_matrix = matGen.SigProfilerMatrixGeneratorFunc(
        project = "CLL", 
        reference_genome = "GRCh37", 
        path_to_input_file = "/orfeo/scratch/cdslab/kdavydzenka/CLL/input")

#chose the data type that you would like to import: "vcf" or "matrix"
#data = sig.importdata("matrix")


# Perform model fitting
sig.sigProfilerExtractor("matrix", "results", input_data, 
                         minimum_signatures = 1, 
                         maximum_signatures = 10, 
                         nmf_replicates = 100)


#save the output results
source_dir = "/orfeo/scratch/cdslab/kdavydzenka/results"
dest_dir = "/orfeo/scratch/cdslab/kdavydzenka/SIGPROFILER/results"
shutil.copytree(source_dir, dest_dir)
