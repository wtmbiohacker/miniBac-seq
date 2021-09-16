## Used to calculate the activities of transcriptional modules, including sigma factors and transcription factors
## Usage: python assess_module_activity.py -i TPM_rm_ncRNA.csv -t network_tf_gene.txt -s network_sigma_gene.txt -n log2 -f 5 -r 3 -o ./
## Written by Tianmin Wang
## Jintao Liu Lab, Tsinghua University
## Last updated May 06, 2021

from __future__ import print_function
import os
import sys
import numpy as np
import pandas as pd
import argparse

# process the argument
parser = argparse.ArgumentParser(description="assess module activity")
parser.add_argument('-i','--tpm_matrix', required=True, help='matrix of tmp of each gene at each condition')
parser.add_argument('-t','--tf_network', required=True, help='regulon information for each TF')
parser.add_argument('-s','--sigma_network', required=True, help='regulon information for each sigma factor')
parser.add_argument('-n','--normalization', default='log2', help='method to normalize the tpm, accepted method includes log2, log10 and no')
parser.add_argument('-f','--filter', type=int, default=10, help='threshold to filter tpm matrix, genes with tpm < this threshold at all conditions will be deleted')
parser.add_argument('-r','--regulon_size', type=int, default=3, help='size of the regulon used to filter those with regulated genes less than this threshold')
parser.add_argument('-o','--output_dir', default='./', help='directory for output')

args=parser.parse_args()

tpm_matrix = pd.read_csv(args.tpm_matrix, sep=',', index_col=0)
tf_network = args.tf_network
sigma_network = args.sigma_network
normalization = args.normalization
tpm_filter = args.filter
regulon_size = args.regulon_size
output_dir = args.output_dir
if normalization not in ['log2', 'log10', 'no']:
    print('incorrect normalization method to process tpm matrix!')
    sys.exit(-1)

# filter and normalize the tpm_matrix
filtered_genes = tpm_matrix[tpm_matrix >= tpm_filter].dropna(axis=0, how='all').index.tolist()
filtered_tpm_matrix = tpm_matrix[tpm_matrix.index.isin(filtered_genes)]

if normalization == 'log10':
    myData = filtered_tpm_matrix.apply(lambda x: np.log10(x + 1))
elif normalization == 'log2':
    myData = filtered_tpm_matrix.apply(lambda x: np.log2(x + 1))
else:
    myData = filtered_tpm_matrix

processed_tpm_matrix = output_dir + 'filtered_normalized_tpm_matrix.csv'
myData.to_csv(processed_tpm_matrix, sep=',')

# perform module activity calculation and visualization
os.system('python TF_activity.py -i %s -o %s -t %s -c %d'%(processed_tpm_matrix, output_dir, tf_network, regulon_size))
os.system('python sigma_activity.py -i %s -s %s -o %s'%(processed_tpm_matrix, sigma_network, output_dir))
os.system('Rscript plot_module_activity.R --sigma_activity sigma_activity_split.xlsx --TF_activity TF_activity.xlsx --output ./')

# relocate the output
os.system('mkdir module_activity/')
os.system('mv sigma*.pdf TF*.pdf sigma*xlsx TF*.xlsx filtered_normalized_tpm_matrix.csv module_activity/')
os.system('rm Rplots.pdf')
