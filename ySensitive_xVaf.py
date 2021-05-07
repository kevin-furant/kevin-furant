#! /usr/bin/env python
#-*-coding:utf-8-*-

import sys,os

import numpy as np
import pandas as pd

def get_skiprows(vcf_file):
    skip_rows = 0
    with open(vcf_file,"r") as inf:
        if each.strip().startswith("#"):
            skip_rows += 1
    return skip_rows

def get_median(x):
    x = str(x).replace('(',"")
    x = str(x).replace(']',"")
    a,b = x.split(",")
    avg = np.mean([float(a),float(b)])
    return round(avg,3)

class sensitive_caculate(object):
    """
    vcf_file_list两列信息，第一列为8个样本的vcf路径，第二列为对应的公司_平台_样本名
    比如  nanda_illumina_E40, purity的值 比如 P40
    """
    def __init__(self,answer_file,vcf_file_list,purity):
        self.__answer = answer_file
        self.__list = vcf_file_list
        self.__bins = self.create_bins()
        self.__purity = purity
        self.__answer_df = self.get_answer_df
        self.__answer_len = self.__answer_df.shape[0]
        self.__vcf_dfs = self.get_vcf_dfs
    
    @staticmethod
    def create_bins():
        m = np.arange(0,1,0.001)
        n = np.arange(0.001,1.001,0.001)
        bins = [(round(x,3),round(y,3)) for (x,y) in zip(m,n)]
        interval_index = pd.IntervalIndex.from_tuples(bins)
        return interval_index

    @property
    def get_vcf_dfs(self):
        with open(self.__list,"r") as inf:
            for each in inf:
                vcf_file,prefix = each.strip().split("\t")
                skip_rows = get_skiprows(vcf_file)
                skip_rows = 1 if skip_rows == 0 else skip_rows
                vcf_df = pd.read_csv(vcf_file,sep="\t",skiprows=skip_rows,header=None,low_memory=False)
                vcf_cols_num = vcf_df.shape[1]
                vcf_df = vcf_df.iloc[:,[0,1,vcf_cols_num-1]]
                vcf_df.columns = ['chrom','pos','vaf']
                vcf_df.vaf = vcf_df.vaf.map(lambda x:x.split(':')[4])
                vcf_df.chrom = vcf_df.chrom.astype(np.string_)
                vcf_df['colcat'] = vcf_df.apply(lambda x: tuple([x.chrom,x.pos]),axis=1)
                vcf_factor = pd.cut(vcf_df.vaf,bins=self.__bins)
                vcf_df["tp"] = vcf_df.groupby(vcf_factor)["colcat"].apply(self.get_tp)
                vcf_df.tp.replace({0:np.nan},inplace=True)
                vcf_df = vcf_df[["tp"]]
                vcf_df.rename(columns={"tp":prefix + "_sensitive"},inplace=True)
                self.__vcf_dfs.append(vcf_df)
                
    @property
    def get_answer_df(self):
        answer_df = pd.read_csv(self.__answer,sep="\t")
        answer_cols_num = answer_df.shape[1]
        answer_df = answer_df.iloc[:,[0,1,answer_cols_num-1]]
        answer_df.columns = ['chrom','pos','vaf']
        answer_df.chrom = answer_df.chrom.astype(np.string_)
        answer_df['colcat'] = answer_df.apply(lambda x: tuple([x.chrom,x.pos]),axis=1)
        return answer_df
    
    @property
    def get_tp(x):
        return folat(x.colcat.isin(self.__answer_df.colcat)).sum()

    @staticmethod
    def get_mean(x):
        tp_min = np.min()
        x.fillna(tp_min)
        y=1-(x/self.__answer_len)
        return y

    def concat_dfs(self):
        merge_df = pd.concat(self.__vcf_dfs,axis=1)
        merge_df = merge_df.apply(self.get_mean,axis=1)
        col_names = merge_df.columns
        for each in col_names:
            each_df = merge_df.loc[:,each]
            each_df.reset_index(inplace=True)
            each_df.rename(columns={each_df.columns[0]:"bins"},inplace=True)
            each_df["vaf"] = each_df.bins.map(get_median)
            each_df["purity"] = self.__purity
            each_df = each_df["vaf","bins",each,"purity"]
            each_df.to_csv(each + ".rst",sep="\t",index=False)

if __name__ == "__main__":
    answer_file = sys.argv[1]
    vcf_list = sys.argv[2]
    purity = sys.argv[3]
    obj = sensitive_caculate(answer_file,vcf_list,purity)