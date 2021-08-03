#!/usr/bin/evn python3
import sys,os
import subprocess
import shutil
import re
import pandas as pd
from multiprocessing import Pool
from configparser import ConfigParser
x = sys.argv[1] #Male_marker.txt
y=  sys.argv[2] #all_type
out=sys.argv[3] #sample.out
cfg = ConfigParser()
cfg.read("parameters.ini")
thread=int(cfg.get("parameter","THREADS"))
def callSSR(i):
    file=open(i,"r")
    lines = file.readlines()
    for line in lines:
        tmp=line.split("\t")
        gene=tmp[0]
        my_samtoools(gene)
    file.close
def my_samtoools(gene):
    cmd=' grep  '+gene+' '+y+'>>'+out
    run_command(cmd)
def main():
    df=pd.read_csv(x,sep='\t',header=None)
    len_df=len(df)
    numble=int(len_df/thread) #thread=20
    gtrd=pd.read_csv(x,chunksize=numble,sep='\t',header=None)
    a=0
    order=[]
    for i in gtrd :
        a=a+1
        order.append(str(a)+'.tmptxt')
        i.to_csv(f"{str(a)}.tmptxt",sep='\t',header=None,index=False)
    p=Pool(thread)
    p.map(callSSR,order)
    p.close()
    p.join()
    df1 = pd.read_csv(x,sep='\t',header=None)
    df2 = pd.read_csv(out,sep='\t',header=None)
    df = pd.merge(df1,df2,on=0,how="left")
    df=df.drop_duplicates()
    df.to_csv(x,sep='\t',header=False,index=False)
    for a in order:
        os.remove(a)
def run_command(cmd):
        print(cmd)
        return_code = subprocess.call(cmd, shell=True)

if __name__ == '__main__' :
    main()

