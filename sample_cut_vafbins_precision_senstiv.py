#! /usr/bin/env python
#-*-coding:utf-8-*-

import sys,os
from bisect import bisect

import numpy as np
import pandas as pd

def get_skiprows(vcf_file):
    skip_rows = 0
    with open(vcf_file,"r") as inf:
        for each in inf:
            if each.strip().startswith("#"):
                skip_rows += 1
    return skip_rows -1
def create_bins():
    "此处不传参了，固定以0.001为step"
    m = np.arange(0,1,0.001)
    n = np.arange(0.001,1.001,0.001)
    bins = [(round(x,3),round(y,3)) for (x,y) in zip(m,n)]
    interval_index = pd.IntervalIndex.from_tuples(bins)
    return interval_index
def mycal(answer_df,vcf_df):
    TP = float((vcf_df.colcat.isin(answer_df.colcat)).sum())
    FN = float((~answer_df.colcat.isin(vcf_df.colcat)).sum())
    FP = float((~vcf_df.colcat.isin(answer_df.colcat)).sum())
    return [TP,FN,FP]

def cal_preci(x):
    try:
        PRECI = x.CUMSUM_TP/(x.CUMSUM_TP + x.CUMSUM_FP)
    except ZeroDivisionError:
        PRECI = 0
    return PRECI

def cal_recall(x):
    try:
        SENSITIVE = x.CUMSUM_TP/(x.CUMSUM_TP + x.CUMSUM_FN)
    except ZeroDivisionError:
        SENSITIVE = 0
    return SENSITIVE

def cal_fscore(x):
    try:
        FSCORE =  2 * x.PRECI * x.SENSITIVE/(x.PRECI+x.SENSITIVE)
    except ZeroDivisionError:
        FSCORE = 0
    return FSCORE

def sort_func(x):
    x = str(x).replace('(',"")
    x = x.split(",")[0]
    return float(x)

def get_median(x):
    x = str(x).replace('(',"")
    x = str(x).replace(']',"")
    a,b = x.split(",")
    avg = np.mean([float(a),float(b)])
    return round(avg,3)
    
def main(answer_file,vcf_file,outfile,sample_name,panel_name,purity):
    "按vaf0.01的梯度去计算准确度与敏感度"
    skip_rows = get_skiprows(vcf_file)
    answer_df = pd.read_csv(answer_file,sep="\t")
    vcf_df = pd.read_csv(vcf_file,sep="\t",skiprows=skip_rows,low_memory=False)
    answer_cols_num = answer_df.shape[1]
    vcf_cols_num = vcf_df.shape[1]
    answer_df = answer_df.iloc[:,[0,1,answer_cols_num-1]]
    vcf_df = vcf_df.iloc[:,[0,1,vcf_cols_num-1]]
    answer_df.columns = ['chrom','pos','vaf']
    vcf_df.columns = ['chrom','pos','Format']
    vcf_df['vaf'] = vcf_df.Format.map(lambda x:x.split(':')[4])
    vcf_df = vcf_df[['chrom','pos','vaf']]
    answer_df['colcat'] = answer_df.apply(lambda x: tuple([x.chrom,x.pos]),axis=1)
    vcf_df['colcat'] = vcf_df.apply(lambda x: tuple([x.chrom,x.pos]),axis=1)
    interval_index = create_bins()
    vcf_factor = pd.cut(vcf_df.vaf,bins=interval_index)
    vcf_dict = dict(list(vcf_df.groupby(vcf_factor)))
    result_df_index = []
    result_array = []
    vcf_dict_keys = vcf_dict.keys()
    sorted_dict_keys = sorted(vcf_dict_keys,key=sort_func)
    throwed_keys = []
    for each_key in sorted_dict_keys:
        vcf_bins_df = vcf_dict[each_key]
        if not vcf_bins_df.empty:
            result_df_index.append(each_key)
            cal_result = mycal(answer_df,vcf_bins_df)
            vcf_pos_count = vcf_bins_df.shape[0]
            cal_result.append(vcf_pos_count)
            result_array.append(cal_result)
        else:
            throwed_keys.append(each_key)
    result_df = pd.DataFrame(np.array(result_array),index=result_df_index,columns=["TP","FN","FP","COUNT"])
    result_df["CUMSUM_TP"] = result_df.TP.cumsum()
    result_df["CUMSUM_FN"] = answer_df.shape[0] - result_df["CUMSUM_TP"]
    result_df["CUMSUM_FP"] = result_df.COUNT - result_df.TP
    result_df["PRECI"] = result_df.apply(cal_preci,axis=1)
    result_df["SENSITIVE"] = result_df.apply(cal_recall,axis=1)
    result_df["FSCORE"] = result_df.apply(cal_fscore,axis=1)   
    result_df.sort_index(ascending=False,inplace=True)
    result_df.reset_index(inplace=True)
    result_df.columns = ["BINS","TP","FN","FP","COUNT","CUMSUM_TP","CUMSUM_FN","CUMSUM_FP","PRECI","SENSITIVE","FSCORE"]
    first_line = result_df.iloc[0,:]
    first_noempty_bin = first_line.iloc[0]
    idx = bisect(throwed_keys,first_noempty_bin) + 1
    add_keys = throwed_keys[idx:]
    first_line = first_line.to_frame()
    add_array = []
    for each in add_keys:
        add_array.append(first_line[1:])
    add_df = pd.concat(add_array,axis=1).T
    add_df['BINS'] = [str(x) for x in add_keys]
    add_df.set_index(['BINS'],drop=True,inplace=True)
    add_df.sort_index(ascending=False,inplace=True)
    add_df.reset_index(inplace=True)
    add_df = add_df[["BINS","TP","FN","FP","COUNT","CUMSUM_TP","CUMSUM_FN","CUMSUM_FP","PRECI","SENSITIVE","FSCORE"]]
    result_df = add_df.append(result_df)
    result_df['VAF'] = result_df.BINS.map(get_median)
    result_df['SAMPLE'] = sample_name
    result_df['COMPANY'] = panel_name
    result_df['PURITY'] = purity
    result_df = result_df[['VAF','BINS',"CUMSUM_TP","CUMSUM_FN","CUMSUM_FP","PRECI","SENSITIVE","FSCORE","SAMPLE","COMPANY","PURITY"]]
    result_df.to_csv(outfile,sep="\t",index=False)
            
if __name__ == "__main__":
    answer_file = sys.argv[1]
    vcf_file = sys.argv[2]
    outfile = sys.argv[3]
    sample_name = sys.argv[4]
    panel_name = sys.argv[5]
    purity = sys.argv[6]
    main(answer_file,vcf_file,outfile,sample_name,panel_name,purity)
