#!/usr/bin/evn python3
import sys,os,argparse
import subprocess
import shutil
import re
import pandas as pd
from multiprocessing import Pool
from configparser import ConfigParser
cfg = ConfigParser()
x = sys.argv[1] #female.01.txt
y=  sys.argv[2] #male.vcf.gz
out=sys.argv[3] #male.out
cfg.read("parameters.ini")
thread=cfg.getint("parameter","THREADS")
bcftools=cfg.get("folders","BCFTOOLS_FOLD")+'/bcftools'
DP=cfg.get("parameter","ALLELE_DEPTH")
DP=int(DP)
GQ=cfg.get("parameter","GQ")
def run_command(cmd):
       # print(cmd)
        return_code = subprocess.call(cmd, shell=True)
def callSSR(i):
    with open(i,"r") as f:
        lines = f.readlines()
        for line in lines:
            tmp=line.strip().split("\t")
            a=tmp[3]
            chrname=tmp[0]+'-'+str(int(tmp[0].split(":")[1])+max(len(a.split('/')[0]),len(a.split('/')[1]))-1)
            myout=chrname+'.txt'
            cmd=bcftools+' view -r '+chrname+' '+y+' >'+myout
            run_command(cmd)
            ssr=tmp[1].split("(")[0]
            motif=a.split("/")[0]
            motif2=a.split("/")[1]
            GMSSR(myout,tmp[0],ssr,motif,motif2)
            os.remove(myout)
def find_maxlen(s,a):
	c=re.finditer(r'(%s)+'%(a),s)
	find_maxlist=[]
	for i in c:
		find_maxlist.append(i.group())
	str2=max(find_maxlist,key=len, default='')
	return str2
def GMSSR(geneout,name,ssr,motif,motif2):
    flag=0
    flag1=0
    a=int(geneout.split('-')[1].replace(".txt",""))-int(geneout.split('-')[0].split(':')[1])+1
    stra=""
    with open(geneout,"r") as f1:
        for line in f1:
            line=line.strip()
            if line.startswith("#"):
                continue
            else:
                tmp=line.split("\t")
                if(tmp[7].split(';')[0]=='INDEL'):
                    if( int(tmp[7].split(';')[4].split('=')[1].split(',')[-1])>=DP and tmp[3].find(ssr)):
                        if(tmp[9].split(':')[0]=='0/0' ):
                            mydtct0=find_maxlen(tmp[3],ssr)
                            if(mydtct0==motif):
                                myout=name+'\t'+mydtct0+'/'+mydtct0+'\t'+'aa'
                                with open(out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='0/1' ):
                            mydtct0=find_maxlen(tmp[3],ssr)
                            mydtct2=find_maxlen(tmp[4],ssr)
                            if(motif==mydtct0 and motif2==mydtct2):
                                myout=name+'\t'+mydtct0+'/'+mydtct2+'\t'+'ab'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                            elif(motif==mydtct0 and motif2!=mydtct2 and mydtct2.find(ssr)):
                                myout=name+'\t'+mydtct0+'/'+mydtct2+'\t'+'ac'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='1/1' ): 
                            mydtct2=find_maxlen(tmp[4],ssr)
                            if(motif2==mydtct2):
                                myout=name+'\t'+mydtct2+'/'+mydtct2+'\t'+'aa'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                            elif(motif2!=mydtct2 and mydtct2.count(ssr)):
                                myout=name+'\t'+mydtct2+'/'+mydtct2+'\t'+'cc'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='1/2' ):
                            mydtct0=find_maxlen(tmp[4].split(",")[0],ssr)
                            mydtct2=find_maxlen(tmp[4].split(",")[1],ssr)
                            if(mydtct0==motif2 or mydtct2==motif2):
                                myout=name+'\t'+mydtct0+'/'+mydtct2+'\t'+'ac'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                            elif((mydtct0!=motif2 or mydtct2!=motif2) and mydtct0.count(ssr) and mydtct2.count(ssr)):
                                myout=name+'\t'+mydtct0+'/'+mydtct2+'\t'+'cd'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)

                else:
                    if(tmp[9].split(':')[0]=='0/0' and int(tmp[9].split(':')[1])>=DP):
                        flag=flag+1
                        if(flag==a):
                            myout=name+'\t'+motif+'/'+motif+'\t'+'aa'
                            with open (out,"a+") as f2:
                                f2.write(myout+"\n")
                                return (1)
                    if(tmp[7].split(';')[0]!='INDEL' and int(tmp[9].split(':')[-1])>=DP ):
                        if(tmp[9].count(':')==3):
                            if(int(tmp[9].split(":")[3])<int(GQ)):
                                continue
                        if(tmp[9].split(":")[0]=='1/1'):
                            flag1=flag1+1                 
                        if(tmp[4]=="."):
                            stra=stra+tmp[3]
                        else:
                            stra=stra+tmp[4]
    str1=motif
    mynum=int(geneout.split('-')[1].replace(".txt",""))-int(geneout.split('-')[0].split(':')[1])+1
    if(len(stra)<mynum):
        return (0)
    else:
        str2=find_maxlen(stra,ssr)
        if(len(str2)>7 and str1!=str2):
            if(motif2==str2 and flag1!=0):
                name=geneout.split('-')[0]
                myout=name+'\t'+str2+'/'+str2+'\t'+'aa'
                with open (out,"a+") as f2:
                    f2.write(myout+"\n")
     #           print(myout)
            elif(motif2==str2 and flag1==0):
                name=geneout.split('-')[0]
                myout=name+'\t'+str1+'/'+str2+'\t'+'ab'
     #           print(myout)
                with open (out,"a+") as f2:
                    f2.write(myout+"\n")

            

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
    for a in order:
        os.remove(a)

if __name__ == '__main__' :
    main()

