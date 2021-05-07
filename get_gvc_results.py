#! /usr/bin/env python
#-*-coding:utf-8-*-
import os,sys
import re
import argparse

import numpy as np
import pandas as pd

class gvcResultCollector(object):
    "对gvc的结果进行汇总成表格"
    def __init__(self,sample_list,indir,outdir,platform,prefix):
        self.__list = sample_list
        self.__indir = indir
        self.__outdir = outdir
        self.__platform = platform
        self.__prefix = prefix
        self.__dfs = []
        self.__record_dict = {}

    def find_qctxt(self):
        "按照gvc的路径格式去指定到对应的样本路径，N只记录一次"
        with open(self.__list,"r") as inf:
            for each in inf:
                each = each.strip()
                normal_result = os.path.join(self.__indir,each,each + ".QC","N","QC.txt")
                tumor_result = os.path.join(self.__indir,each,each + ".QC","T","QC.txt")
                normal_sample_name = re.sub(r'\d+',"0",each)
                tumor_sample_name = each
                if not self.__record_dict.get(normal_sample_name):
                    self.__record_dict[normal_sample_name] = normal_result
                self.__record_dict[tumor_sample_name] = tumor_result
        
    def txt2df(self):
        self.find_qctxt()
        for each in self.__record_dict:
            sample = "_".join([self.__prefix,self.__platform,each])
            df = pd.read_csv(self.__record_dict[each],sep="\t",header=None,names=["idx",sample])
            df.set_index(df.idx,drop=False,inplace=True)  #测试的时候发现drop=True不起作用，所以索性设置为False
            df.drop(columns=["idx"],inplace=True)
            self.__dfs.append(df)

    def out2csv(self):
        self.txt2df()
        merge_df = pd.concat(self.__dfs,axis=1)
        out_file = os.path.join(self.__outdir,self.__prefix + "_" + self.__platform  + "_merge_gvc_result.xls")
        merge_df.to_csv(out_file,sep="\t")

def init_args():
    parser = argparse.ArgumentParser(  
        description='收集gvc结果汇总成大表')  
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
    obj = gvcResultCollector(sample_list,indir,outdir,platform,prefix)
    obj.out2csv()
