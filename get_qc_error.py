#! /usr/bin/env python
#-*-coding:utf-8-*-

import os,sys
import argparse

import numpy as np
import pandas as pd

class samplesQmConcat(object):
    "对质控结果中的QM文件进行汇总，并将所有样本的数据合并成一个文件"
    def __init__(self,sample_list,indir,outdir,platform,prefix):
        self.__list = sample_list
        self.__indir = indir
        self.__outdir = outdir
        self.__platform = platform
        self.__prefix = prefix
        self.__dfs = []
        self.__record_dict = {}

    def search_qm(self):
        "按照目录格式生成QM文件对应路径的字典"
        with open(self.__list,"r") as inf:
            for each in inf:
                each = each.strip()
                qm_file = os.path.join(self.__outdir,each,"raw_" + each + ".QM")
                self.__record_dict[each] = qm_file

    def txt2df(self):
        self.search_qm()
        for each in self.__record_dict:
            df = pd.read_csv(self.__record_dict[each],sep="\t",header=None,names=["cycle","quality","error"])
            df = df.iloc[:,[0,2]]
            self.__dfs.append(df)

    def out2csv(self):
        self.txt2df()
        merge_df = pd.concat(self.__dfs,ignore_index=True)
        merge_df.sort_values(by=['cycle'],inplace=True)
        out_file = os.path.join(self.__outdir,self.__prefix + "_" + self.__platform  + "_merge_error.xls")
        merge_df.to_csv(out_file,sep="\t",index=False)

def init_args():
    parser = argparse.ArgumentParser(  
        description='收集QC中QM结果汇总成大表')  
    parser.add_argument(
        '--list','-l',type=str,required=True,help='给定样本列表')  
    parser.add_argument(  
        '--indir','-i',type=str,required=True,help='给定存放gvc结果的外层路径')  
    parser.add_argument(
        '--outdir','-o',type=str,default=os.getcwd(),help='给定输出文件生成的路径'
    )
    parser.add_argument(
        '--platform','-p',type=str,required=True,help='给定平台信息，illumina或者bgi'
    )
    parser.add_argument(
        '--prefix','-r',type=str,required=True,help='给定输出公司名称前缀，boke|nanda|twist|ajtk'
    )
    args = parser.parse_args() 
    return args

if __name__ == "__main__":
    args = init_args()
    sample_list = args.list
    indir = args.indir
    outdir = args.outdir
    platform = args.platform
    prefix = args.prefix
    obj =samplesQmConcat(sample_list,indir,outdir,platform,prefix)
    obj.out2csv()