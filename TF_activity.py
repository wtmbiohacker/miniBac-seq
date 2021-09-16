'''

Chai Ruochen
2021-02-21

This script is to calculate activity of transcription factor from an expression and transcription factor info matrix.

To execute:
python TF_activity.py -i <gene expression file.csv> -o <out_directory> -t <TF_info.txt(network_tf_gene.txt)> -c <No. of targets to cutoff>

'''

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="calculate TF activity")
parser.add_argument("-i",type=str,required=True,dest="gene_matrix",default=None,
                    help="gene expression matrix, csv file.")
parser.add_argument("-t",type=str,required=True,dest="tf_info",default=None,
                    help="TF info, txt file, downloaded from RegulonDB.")
parser.add_argument("-c",type=int,required=False,dest="target_cutoff",default=5,
                    help="cutoff of the number of targets a TF has")
parser.add_argument("-o",type=str,required=True,dest="out_dir",default='',
                    help="output directory.")

args=parser.parse_args()

if args.out_dir == '':
    OUT_DIR = os.getcwd()
else:
    OUT_DIR = args.out_dir

TF_info_file=args.tf_info
Gene_expression_matrix=args.gene_matrix
targets_cutoff=args.target_cutoff

#load file
df_expression = pd.read_csv(Gene_expression_matrix,index_col=0)
df_tf=pd.read_csv(TF_info_file,sep = "\t",header = None)
df_tf.drop(columns=[7],inplace=True)
df_tf.columns = ["TF_ID","TF_name","gene_ID","gene_name","effect","evidence","evidence_type"]
df_tf["effect"] = df_tf["effect"].replace("-","neg").replace("+","pos")
df_tf["TF_name"]=['_'.join(i) for i in df_tf[["TF_name","effect"]].values]

#calculate activity
df_merge = pd.merge(df_expression,df_tf,right_on=["gene_name"],left_index=True)
tf_n_targets = df_merge.groupby("TF_name").size()
tf_list = tf_n_targets.index[tf_n_targets >= targets_cutoff]
df_merge = df_merge[df_merge["TF_name"].isin(tf_list)]

df_tf_activity_median = df_merge.groupby("TF_name").median()
df_tf_activity_mean = df_merge.groupby("TF_name").mean()

#calculate rank
gene_number = df_expression.shape[0]
def get_rank(row):
    df_rank = df_expression.append(row).rank(method="min") / gene_number
    return df_rank.iloc[-1,]

df_tf_rank_mean = df_tf_activity_mean.apply(get_rank, axis=1)
df_tf_rank_median = df_tf_activity_median.apply(get_rank, axis=1)

# output
out_file = os.path.join(OUT_DIR,'TF_activity.xlsx')
with pd.ExcelWriter(out_file) as writer:
    df_tf_activity_mean.to_excel(excel_writer=writer, sheet_name="tf_activity_mean", index=True)
    df_tf_rank_mean.to_excel(excel_writer=writer, sheet_name="quantile_mean", index=True)
    df_tf_activity_median.to_excel(excel_writer=writer, sheet_name="tf_activity_median", index=True)
    df_tf_rank_median.to_excel(excel_writer=writer, sheet_name="quantile_median", index=True)
    df_merge[["gene_name","TF_name"]].groupby("TF_name").size().rename("n_target").to_excel(
        excel_writer=writer, sheet_name="n_target", index=True)
