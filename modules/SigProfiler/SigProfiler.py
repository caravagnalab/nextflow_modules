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
joint_table = "/u/cdslab/kdavydzenka/mutationsTable.tsv"
input_data = pd.read_csv(joint_table, sep = '\t')
input_path = "/u/cdslab/kdavydzenka/CLL/input_multisample/"

#input data preprocessing
def input_processing(data):
    new_columns = {'Project': "CLL", 'Genome': 'GRCh37', 'Type': "SOMATIC", 'mut_type': "SNP"}
    df = data.assign(**new_columns)
    df['chr'] = df['chr'].astype(str).str[3:]
    df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
    df["ID"] = df["Sample"]
    df = df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
    #df = df.style.hide(axis='index') #omit row indexes
    return df

input_data = input_processing(input_data)

#saving input matrix to txt
input_data.to_csv('CLL/input_multisample/input_data.txt', sep='\t', index=False, header=True)

#mutation's counts matrix generation
input_matrix = matGen.SigProfilerMatrixGeneratorFunc(
        project = "CLL", 
        reference_genome = "GRCh37", 
        path_to_input_files = input_path)

#chose the data type that you would like to import: "vcf" or "matrix"
#data = sig.importdata("matrix")

output_path = "output/SBS/CLL.SBS96.all"
# Perform model fitting
sig.sigProfilerExtractor(input_type = "matrix", 
                         output = "results", 
                         input_data = input_path+output_path,  #path to the file
                         exome = False,
                         minimum_signatures = 1,  #default value
                         maximum_signatures = 10, #default value 25
                         nmf_replicates = 100, #the number of iteration to be performed to extract each number signature
                         resample = True,
                         matrix_normalization = "gmm", #default
                         seeds= "random",
                         nmf_init = "random", #he initialization algorithm for W and H matrix of NMF
                         min_nmf_iterations = 10000, #minimum number of iterations to be completed before NMF converges (default 10000)
                         max_nmf_iterations = 1000000, #default
                         nmf_test_conv = 10000, #defines the number of iterations to done between checking next convergence
                         nmf_tolerance = 1e-15, #tolerance to achieve to converge
                         cpu = -1,
                         gpu = False,
                         cosmic_version = 3.4,
                         make_decomposition_plots = True, #denovo to Cosmic sigantures decompostion plots will be created as a part the results
                         collapse_to_SBS96 = True, #SBS288 and SBS1536 Denovo signatures will be mapped to SBS96 reference signatures
                         get_all_signature_matrices = True,
                         export_probabilities = True)


#save the output results
source_dir = "/orfeo/scratch/cdslab/kdavydzenka/results"
dest_dir = "/orfeo/scratch/cdslab/kdavydzenka/SIGPROFILER/results"
shutil.copytree(source_dir, dest_dir)
