#! /usr/bin/env python
#-*-coding:utf-8-*-

import sys,os
import re

import numpy as np
import pandas as pd

def judge_company_value(x):
    if x.TMB_company >= x.threshold:
        return "TMB-H"
    else:
        return "TMB-L"

def judge_answer_value(x):
    if x.answer >= x.threshold:
        return "TMB-H"
    else:
        return "TMB-L"

def stat_compare(x):
    TP = ((x.answer_stat=="TMB-H")&(x.company_stat==x.answer_stat)).sum()
    TN=((x.answer_stat=="TMB-L")&(x.company_stat==x.answer_stat)).sum()
    FN=((x.answer_stat=="TMB-H")&(x.company_stat!=x.answer_stat)).sum()
    FP=((x.answer_stat=="TMB-L")&(x.company_stat!=x.answer_stat)).sum()
    try:
        PRECI = float(TP)/float((TP+FP))
    except ZeroDivisionError:
        PRECI = 0
    try:
        RECALL = float(TP)/float((TP+FN))
    except ZeroDivisionError:
        RECALL = 0
    try:
        FSCORE = 2*(float(PRECI)*float(RECALL))/(float(PRECI)+float(RECALL))
    except ZeroDivisionError:
        FSCORE = 0
    return TP,TN,FN,FP,PRECI,RECALL,FSCORE

def tmb_threshold(status_file,thresh_file,out_file,agg_outfile):
    """
    通过输入的阈值文件对状态文件进行样本的TMB值状态归类，经
    分类处理后输出文件
    """
    st_df = pd.read_csv(status_file,sep="\t")
    thresh_df = pd.read_csv(thresh_file,sep="\t",header=None,names=["groups","threshold"])
    st_df = st_df.iloc[:,[1,3,6,7]]
    st_df['groups'] = st_df.sample_name.map(lambda x:re.sub(r'\d+',"",x))
    merge_df = pd.merge(st_df,thresh_df,on="groups")
    merge_df['company_stat'] = merge_df.apply(judge_company_value,axis=1)
    merge_df['answer_stat'] = merge_df.apply(judge_answer_value,axis=1)
    merge_df['sample_code'] = merge_df.sample_name.map(lambda x:re.sub(r'[A-Z]{1}','',x))
    rst_df = merge_df.groupby(['companyNu','sample_code']).apply(stat_compare)
    agg_df = merge_df.groupby(['companyNu']).apply(stat_compare)
    rst_ser = rst_df.map(lambda x:re.sub(r'\(|\)','',str(x)))
    agg_ser = agg_df.map(lambda x:re.sub(r'\(|\)','',str(x)))
    rst_df = rst_ser.str.split(',',expand=True)
    agg_df = agg_ser.str.split(',',expand=True)
    rst_df.columns = ['TP','TN','FN','FP','PRECI','RECALL','FSCORE']
    agg_df.columns = ['TP','TN','FN','FP','PRECI','RECALL','FSCORE']
    rst_df.reset_index(inplace=True)
    agg_df.reset_index(inplace=True)
    rst_df.to_csv(out_file,sep="\t",index=False)
    agg_df.to_csv(agg_outfile,sep="\t",index=False)
    
if __name__ == "__main__":
    status_file = sys.argv[1]
    thresh_file = sys.argv[2]
    out_file = sys.argv[3]
    agg_outfile = sys.argv[4]
    tmb_threshold(status_file,thresh_file,out_file,agg_outfile)
