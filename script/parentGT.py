#!/usr/bin/evn python3
import sys,os,argparse
import subprocess
import shutil
import re
import pandas as pd
from multiprocessing import Pool
from configparser import ConfigParser
cfg = ConfigParser()
x = sys.argv[1] #SRR.txt
y=  sys.argv[2] #male.vcf.gz
myfile=sys.argv[3] #male.01.txt
cfg.read("parameters.ini")
thread=cfg.getint("parameter","THREADS")
bcftools=cfg.get("folders","BCFTOOLS_FOLD")+'/bcftools'
DP=cfg.get("parameter","ALLELE_DEPTH")
DP=int(DP)
GQ=cfg.get("parameter","GQ")
def ssrdict(i):
    chr_length = {}
    file = open(i, "r")
    lines = file.readlines()
    for line in lines:
        if line.startswith("number"):
            continue
        else:
            tmp=line.split("\t")
            a=tmp[1]+":"+tmp[6]
            b=tmp[2]+"_"+tmp[4]
            chr_length[a]=b
    return (chr_length)
mydict=ssrdict(x)
def run_command(cmd):
       # print(cmd)
        return_code = subprocess.call(cmd, shell=True)
def callSSR(i):
    with open(i,"r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("number"):
                continue
            else:
                tmp=line.strip().split("\t")
                chrname=tmp[1]+':'+tmp[6]+'-'+str(int(tmp[7])+6)
                myout=chrname+'.txt'
                cmd=bcftools+' view -r '+chrname+' '+y+' >'+myout
                run_command(cmd)
                ssr=tmp[2]
                motif=tmp[4]
                GMSSR(myout,ssr,motif)
                os.remove(myout)
def find_maxlen(s,a):
	c=re.finditer(r'(%s)+'%(a),s)
	find_maxlist=[]
	for i in c:
		find_maxlist.append(i.group())
	str2=max(find_maxlist,key=len, default='')
	return str2

def GMSSR(geneout,ssr,motif):
    a=""
    with open(geneout,"r") as f1:
        for line in f1:
            line=line.strip()
            if line.startswith("#"):
                continue
            else:
                tmp=line.split("\t")
                if(tmp[7].split(';')[0]=='INDEL' and tmp[1]==str(int(geneout.split("-")[0].split(":")[1])-1) and (tmp[9].split(":")[0]=='0/1' or tmp[9].split(":")[0]=='1/2')):
                    if (tmp[9].split(":")[0]=='0/1' ):
                        if(len(tmp[3])>7 and int(tmp[9].split(":")[3])>=int(GQ)):
                            name=tmp[0]+":"+str(int(tmp[1])+1)
                            if (name in mydict.keys()):
                                ssr=mydict[name].split("_")[0]
                                motif=mydict[name].split("_")[1]
                                str1=find_maxlen(tmp[3],ssr) #ssr=AT motif=5
                                str2=find_maxlen(tmp[4],ssr)
                                if (str1!=str2 and len(str1)==len(ssr*int(motif)) and str2!=""):
                                    myout=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(str1))+'\t'+str1+'/'+str2+'\t'+'ab'
                                    with open(myfile,"a+") as f2:
                                        f2.write(myout+"\n")
                                        return (1)
                    elif(tmp[9].split(":")[0]=='1/2' ):
                        if(len(tmp[3])>7 and int(tmp[9].split(":")[3])>=int(GQ)):
                            name=tmp[0]+":"+str(int(tmp[1])+1)
                            if (name in mydict.keys()):
                                ssr=mydict[name].split("_")[0]
                                motif=mydict[name].split("_")[1]
                                str1=find_maxlen(tmp[4].split(",")[0],ssr)
                                str2=find_maxlen(tmp[4].split(",")[1],ssr)
                                if(str1.count(ssr) and str2.count(ssr)):
                                    myout=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(ssr*int(motif)))+'\t'+str1+'/'+str2+'\t'+'ab'
                                    with open(myfile,"a+") as f2:
                                        f2.write(myout+"\n")
                                        return (1)
                else:
                    if( tmp[7].split(';')[0]!='INDEL' and int(tmp[9].split(':')[-1])>=DP and tmp[9].split(":")[0]!='1/1'): #tmp[9].split(":")[0]=='0/0' and
                        if(tmp[9].count(':')==3):
                            if(int(tmp[9].split(":")[3])<int(GQ)):
                                continue
                        if(tmp[4]=="."):
                            a=a+tmp[3]
                        else:
                            a=a+tmp[4]
    str1=ssr*int(motif)
    mynum=int(geneout.split('-')[1].replace(".txt",""))-int(geneout.split('-')[0].split(':')[1])+1
    if(len(a)<mynum):
        return (0)
    else:
        str2=find_maxlen(a,ssr)
        if(len(str2)>7 and str1!=str2):
            name=geneout.split('-')[0]
            myout=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(ssr*int(motif)))+'\t'+str1+'/'+str2+'\t'+'ab'
        #    print(myout)
            with open(myfile,"a+") as f2:
                f2.write(myout+"\n")
                    
    

def main():
    file = open(myfile, "w")
    file.close()
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

