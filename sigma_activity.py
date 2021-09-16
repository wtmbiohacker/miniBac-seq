'''

Chai Ruochen
2021-02-04

edit:2021-02-20

calculate sigma factor activity from an expression and a sigma infomation matrix(download from regulon DB)

To execute:
python sigma_activity.py -i <gene_expression_matrix.csv> -s <sigma_network.txt(network_sigma_gene.txt)> -o <out directory>

'''
import pandas as pd
import argparse
import os


parser = argparse.ArgumentParser(description= "calculate sigma factor activity")
#parser.add_argument("-i",type=str,default = None,required = True,dest = "gene_matrix",
#                    help="gene raw counts matrix csv file with genes on rows and samples on columns")
parser.add_argument("-i",type=str,default = None,required = True,dest = "gene_matrix",
                    help = "gene expression matrix")
parser.add_argument("-s",type=str,default = None,required = True,dest = "sigma_info",
                    help = "sigma regulation network file downloaded from regulonDB")
parser.add_argument("-o",type=str,default = "",dest = "out_dir",
                    help="Path to output file default as current directory")
args = parser.parse_args()

if args.out_dir == '':
    OUT_DIR = os.getcwd()
else:
    OUT_DIR = args.out_dir
Sigma_info_file = args.sigma_info
Gene_expression_matrix = args.gene_matrix


#load file
df_expression = pd.read_csv(Gene_expression_matrix,index_col=0)
df_sig_net = pd.read_csv(Sigma_info_file,sep = "\t",header = None)
df_sig_net.columns = ["sigma_factor","target","effect","evidence","evidence_type"]

for i in ["all","split"]:
    if i == "split": # split genes regulated by multiple sigma factors equally to each sigma factor
        df_sig_net = \
            df_sig_net.drop(columns = ["sigma_factor"]).\
                join(df_sig_net["sigma_factor"].str.split(", ",expand=True).stack().reset_index(level = 1,drop = True).
                rename("sigma_factor")).reset_index(drop = True).drop_duplicates(subset = ["target","sigma_factor"])
    #calculate activity by mean and median
    gene_number = df_expression.shape[0]
    df_merge = pd.merge(df_expression,df_sig_net,left_index=True,right_on="target")
    df_sigma_activity_mean = df_merge.groupby("sigma_factor").mean()
    df_sigma_activity_median = df_merge.groupby("sigma_factor").median()

    #calculate rank
    def get_rank(row):
        df_rank = df_expression.append(row).rank(method="min")/gene_number
        return df_rank.iloc[-1,]

    df_sigma_rank_mean = df_sigma_activity_mean.apply(get_rank,axis=1)
    df_sigma_rank_median = df_sigma_activity_median.apply(get_rank,axis=1)
    out_file = os.path.join(OUT_DIR,'sigma_activity_{0}.xlsx'.format(i))
    with pd.ExcelWriter(out_file) as writer:
        df_sigma_activity_mean.to_excel(excel_writer=writer, sheet_name="sigma_activity_mean", index=True)
        df_sigma_rank_mean.to_excel(excel_writer=writer, sheet_name="quantile_mean", index=True)
        df_sigma_activity_median.to_excel(excel_writer=writer, sheet_name="sigma_activity_median", index=True)
        df_sigma_rank_median.to_excel(excel_writer=writer, sheet_name="quantile_median", index=True)
        df_merge[["target","sigma_factor"]].groupby("sigma_factor").size().rename("n_target").to_excel(
            excel_writer=writer, sheet_name="n_target", index=True)
