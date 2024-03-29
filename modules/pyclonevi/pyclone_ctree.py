#!/usr/bin/env python3

import sys
import pandas as pd

## Read pyclone best fit table
best_fit_file=sys.argv[2]
df_output = pd.read_csv(best_fit_file, sep='\t')

## Read pyclone input table

df_input = pd.read_csv(sys.argv[1], sep = '\t')

## Caluclate number of mutations per cluster and add to the subset orginal df_output dataframe

df_output_small = df_output[['mutation_id','sample_id','cluster_id','cellular_prevalence']]
df_output_small['nMuts'] = df_output_small.groupby('cluster_id')['mutation_id'].transform('nunique')

## Find the clonal cluster and add the colum 'is.clonal' to identify it

samples = df_output_small['sample_id'].unique()
top_clusters = 0
for s in samples:
  ind = df_output_small[df_output_small['sample_id']==s]['cellular_prevalence'].idxmax()
  top_cluster = df_output_small.loc[ind,'cluster_id']
  if (top_clusters!=top_cluster):
    top_clusters=top_cluster

i = df_output_small[df_output_small['cluster_id']==top_clusters].index

df_output_small['is.clonal']="F"
df_output_small.loc[i,'is.clonal']="T"

## Merge the input information about driver genes to the clusters
df_merged =pd.merge(df_output_small, df_input[['mutation_id','sample_id','patientID','variantID','is.driver']], on = ["mutation_id","sample_id"], how="inner")
df_final = df_merged.drop_duplicates(subset=['sample_id', 'variantID','cluster_id','is.clonal','is.driver'])

## Find clusters with driver genes and find cluster wihtout driver genes
## For the ones without drivers add the 'NA' to variantID and False to is.driver
driver_genes_df = df_final[df_final['is.driver'] == True]

non_driver_genes_df = df_final[df_final['is.driver'] == False].drop_duplicates(subset=['sample_id','cluster_id']).copy()
non_driver_genes_df['variantID'] = 'NA'
non_driver_genes_df['is.driver'] = False
result_df = pd.concat([driver_genes_df, non_driver_genes_df], ignore_index=True)

# Identify clusters for which is.driver is both True and False
conflicting_clusters = result_df.groupby(['sample_id', 'cluster_id'])['is.driver'].nunique() == 2
conflicting_clusters = conflicting_clusters[conflicting_clusters].index

# Create a new DataFrame with only True values for conflicting clusters
true_values_df = result_df[(result_df['is.driver'] == True) & result_df.set_index(['sample_id', 'cluster_id']).index.isin(conflicting_clusters)]

# Create a new DataFrame with all values for non-conflicting clusters
non_conflicting_df = result_df[~result_df.set_index(['sample_id', 'cluster_id']).index.isin(conflicting_clusters)]

# Concatenate the DataFrames
result_df = pd.concat([true_values_df, non_conflicting_df], ignore_index=True)
result_df['is.driver'] = result_df['is.driver'].replace({False: 'F', True: 'T'})

result_df = result_df.rename(columns={'sample_id': 'sampleID', 'cellular_prevalence':'CCF','cluster_id': 'cluster'})
result_df['tool'] = pd.Series(["pyclonevi" for x in range(len(result_df.index))])
ctree_input = result_df[['patientID','variantID','is.driver','is.clonal','cluster','nMuts','sampleID','CCF','tool']]



final_jt_file = sys.argv[3]
ctree_input.to_csv(final_jt_file, sep=",",index=False, header=True)


