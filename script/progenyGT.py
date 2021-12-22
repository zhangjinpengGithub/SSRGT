#!/usr/bin/evn python3
import sys,os
import subprocess
import shutil
import re
import pandas as pd
from multiprocessing import Pool
from configparser import ConfigParser
x = sys.argv[1] 
y=  sys.argv[2] 
out=sys.argv[3] 
cfg = ConfigParser()
cfg.read("parameters.ini")
thread=cfg.getint("parameter","THREADS")
bcftools=cfg.get("folders","BCFTOOLS_FOLD")+'/bcftools'
DP=cfg.get("parameter","ALLELE_DEPTH")
DP=int(DP)
GQ=cfg.get("parameter","GQ")
def find_maxlen(s,a):
        c=re.finditer(r'(%s)+'%(a),s)
        find_maxlist=[]
        for i in c:
                find_maxlist.append(i.group())
        str2=max(find_maxlist,key=len, default='')
        return str2

def run_command(cmd):
    return_code = subprocess.call(cmd, shell=True)
def callSSR(i):
    file=open(i,"r")
    lines = file.readlines()
    for line in lines:
        tmp=line.strip().split("\t")
        a=tmp[3]
        chrname=tmp[0]+'-'+str(int(tmp[0].split(":")[1])+max(len(a.split('/')[0]),len(a.split('/')[1]))-1)
        myout=chrname+'.txt'
        cmd=bcftools+' view -r '+chrname+' '+y+' >'+myout
        run_command(cmd)
        ssr=tmp[1].split("(")[0]
        motif=a.split("/")[0]
        GMSSR(myout,tmp[0],ssr,motif)
        os.remove(myout)
    file.close
def GMSSR(geneout,name,ssr,motif):
    flag=0
    flag1=0
    stra=""
    a=int(geneout.split('-')[1].replace(".txt",""))-int(geneout.split('-')[0].split(':')[1])+1
    a1=name.split(":")[1]
    with open(geneout,"r") as f1:
        for line in f1:
            line=line.strip()
            if line.startswith("#"):
                continue
            else:
                tmp=line.split("\t")
                if(tmp[7].split(';')[0]=='INDEL' and str(int(tmp[1])+1)==a1):
                    if( int(tmp[7].split(';')[4].split('=')[1].split(',')[-1])>=DP and tmp[3].find(ssr)):
                        if(tmp[9].split(':')[0]=='0/0' ):
                            mydtct0=find_maxlen(tmp[3],ssr)
                            if(mydtct0.count(ssr)):
                                myout=name+'\t'+mydtct0+'/'+mydtct0
                                with open(out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='0/1' ):
                            mydtct0=find_maxlen(tmp[3],ssr)
                            mydtct2=find_maxlen(tmp[4],ssr)
                            if(mydtct0.count(ssr) and mydtct2.count(ssr)):
                                myout=name+'\t'+mydtct0+'/'+mydtct2
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='1/1' ):
                            mydtct2=find_maxlen(tmp[4],ssr)
                            if(mydtct2.count(ssr)):
                                myout=name+'\t'+mydtct2+'/'+mydtct2
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='1/2' ):
                            mydtct0=find_maxlen(tmp[4].split(",")[0],ssr)
                            mydtct2=find_maxlen(tmp[4].split(",")[1],ssr)
                            if(mydtct0.count(ssr) and mydtct2.count(ssr)):
                                myout=name+'\t'+mydtct0+'/'+mydtct2
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                else:
                    if(tmp[9].split(':')[0]=='0/0' and int(tmp[9].split(':')[1])>=DP):
                        flag=flag+1
                        if(flag==a):
                            myout=name+'\t'+motif+'/'+motif
                            with open (out,"a+") as f2:
                                f2.write(myout+"\n")
                            return (1)

    
def main():
    df=pd.read_csv(x,sep='\t',header=None)
    len_df=len(df)
    numble=int(len_df/thread) 
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
    df1 = pd.read_table(x,header=None)
    df2 = pd.read_table(out,header=None)    
    df = pd.merge(df1,df2,on=0,how="left")
    df=df.drop_duplicates()
    df.to_csv(x,sep='\t',header=False,index=False)
    for a in order:
        os.remove(a)
if __name__ == '__main__' :
    main()
