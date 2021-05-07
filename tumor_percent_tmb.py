#! /usr/bin/env python
#-*-coding:utf-8-*-
import os,sys

import numpy as np
import pandas as pd

def cal_bins_tmb(slide_file,maf_bed,panel_length):
    """
    slide文件中有TOP和BOTTOM两个纯度值，要先求个均值，然后去匹配
    maf_bed文件，然后去对纯度排序分桶计算后得到每个桶区间的结果
    """
    slide_df = pd.read_csv(slide_file,sep="\t")
    slide_ser = slide_df.groupby('sample_submitter_id').apply(np.mean)
    slide_df_mean = slide_ser.reset_index()
    maf_df = pd.read_csv(maf_bed,header=None,names=["chr","start","end","sample_submitter_id"],sep="\t",low_memory=False)
    merge_df = pd.merge(maf_df,slide_df_mean,on='sample_submitter_id',how='left')
    all_rows = merge_df.shape[0]
    dropna_merge_df = merge_df.dropna()
    without_na_rows = dropna_merge_df.shape[0]
    na_rows = all_rows - without_na_rows
    sort_dropna_merge_df = dropna_merge_df.sort_values(by='percent_tumor_cells')
    #tumor_percent_uniq_values = sort_dropna_merge_df.percent_tumor_cells.unique()
    percent_tumor_count_df = sort_dropna_merge_df.groupby('percent_tumor_cells')[['sample_submitter_id']].count()
    percent_tumor_count_df.columns = ['count']
    with open("log.txt","a") as outf:
        outf.write(
            "没有纯度值而舍弃的位点个数：" + str(na_rows) + "\n"
        ) 
        outf.write(
            "纯度值分布：" + "\n"
        )
        percent_tumor_count_df.to_csv(outf,sep="\t")
    merge_df.to_csv("all_pos_table.xls",sep="\t",index=False,na_rep="NA")
    sort_dropna_merge_df.to_csv("sort_dropna_merge_table.xls",sep="\t",index=False,na_rep="NA")
    factor = pd.cut(sort_dropna_merge_df.percent_tumor_cells,bins=[0,13,17,23,27,98,100],labels=['a','10','b','20','c','100'])
    cut_bins_count_df = sort_dropna_merge_df.groupby(factor)[["percent_tumor_cells"]].count()
    cut_bins_count_df["tmb_value"] = cut_bins_count_df.percent_tumor_cells/(float(panel_length)/1000000)
    cut_bins_count_df.to_csv("bins_count_tmb.xls",sep="\t")

if __name__ == "__main__":
    slide_file = sys.argv[1]
    maf_bed = sys.argv[2]
    panel_length = sys.argv[3]
    cal_bins_tmb(slide_file,maf_bed,panel_length)
