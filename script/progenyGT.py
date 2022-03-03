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
ref=sys.argv[4]
cfg = ConfigParser()
cfg.read("parameters.ini")
thread=cfg.getint("parameters","THREADS")
bcftools=cfg.get("folders","BCFTOOLS_FOLD")+'/bcftools'
DP=cfg.get("parameters","ALLELE_DEPTH")
DP=int(DP)
GQ=cfg.get("parameters","GQ")
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
        newmyout=chrname+'.txt'
        cmd=bcftools+' view -r '+chrname+' '+y+' >'+newmyout
        run_command(cmd)
        ssr=tmp[1].split("(")[0]
        motif=a.split("/")[0]
        GMSSR(newmyout,tmp[0],ssr,motif)
        if(os.path.isfile(newmyout)):
            os.remove(newmyout)
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
                    if( tmp[3].find(ssr)):
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
                            if(mydtct2.count(ssr) and len(mydtct2)>5):
                                myout=name+'\t'+mydtct2+'/'+mydtct2
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='1/2' ):
                            mydtct0=find_maxlen(tmp[4].split(",")[0],ssr)
                            mydtct2=find_maxlen(tmp[4].split(",")[1],ssr)
                            if(mydtct0.count(ssr) and mydtct2.count(ssr) and len(mydtct2)>5):
                                myout=name+'\t'+mydtct0+'/'+mydtct2
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                    else:
                        return (0)   
                elif(tmp[7].split(';')[0]=='INDEL' and int(tmp[1])+1<int(a1) ):
                    myflag=0
                    str2=""
                    tmpstr=""
                    myout=geneout.replace('.txt','')
                    cmd=bcftools+' norm -f '+ref+' -m +both -r '+myout+' '+y+' 1 >'+myout +" 2>log"
                    run_command(cmd)
                               
                    with open(myout,"r") as f3:
                        tmplist=[]
                        f3=f3.readlines()
                        for line3 in f3:
                            line3=line3.strip()
                            if line3.startswith("#"):
                                continue
                            else:
                                tmp2=line3.split("\t")
                                tmplist.append(tmp2[1])
                                if(tmp2[7].split(';')[0]=='INDEL' and (tmp2[9].split(':')[0]=='0/1' ) and int(tmp2[1])+1<int(a1)):
                                    if(len(tmp2[3])>len(tmp2[4])):
                                        tmpstr=tmp2[3]
                    num=int(tmplist[1])-int(tmplist[0])
                    if(num>0 and len(tmpstr)>num):
                            tmpstr=tmpstr[num:]
                            myflag=len(tmpstr)
                    else:
                        if(os.path.isfile(myout)):
                            os.remove(myout)                       
                        return (0)
                    with open(myout,"r") as f3:
                        
                        f3=f3.readlines()
                        for line3 in f3:
                            line3=line3.strip()
                            if line3.startswith("#"):
                                continue
                            else:
                                tmp2=line3.split("\t")
                                if(tmp2[7].split(';')[0]=='INDEL' and (tmp2[9].split(':')[0]=='0/1')):
                                    str2=tmp2[4]
                                elif(tmp2[9].split(':')[0]=='0/0'):
                                    if(myflag>0):
                                        myflag=myflag-1
                                    elif(myflag==0):
                                        str2=str2+tmp2[3]
                                elif(tmp2[9].split(':')[0]=='1/1'):
                                    str2=str2+tmp2[4]
                    str2=find_maxlen(str2,ssr)
                    if(len(str2)!=len(motif)):
                        myout2=name+'\t'+motif+'/'+str2
                        if(str2==""):
                            if(os.path.isfile(myout)):
                                os.remove(myout)
                                return (0)
                        with open (out,"a+") as f2:
                            f2.write(myout2+"\n")
                            if(os.path.isfile(myout)):
                                os.remove(myout)
                            return (1)
                    if(os.path.isfile(myout)):
                        os.remove(myout)
                                                                                                                 
                    return (0)
                elif(tmp[7].split(';')[0]!='INDEL' and tmp[9].split(':')[0]=='0/0' ):
                    str2=""
                    myflag=0
                    flag11=0
                    indelflag=0
                    myout=geneout.replace('.txt','')
                    cmd=bcftools+' norm -f '+ref+' -m +both -r '+myout+' '+y+' 1 >'+myout +" 2>log"
                    run_command(cmd)
                    with open(myout,"r") as f3:
                        f3=f3.readlines()
                        for line3 in f3:
                            line3=line3.strip()
                            if line3.startswith("#"):
                                continue
                            else:
                                tmp2=line3.split("\t")
                                if(tmp2[9].split(':')[0]=='0/0'):
                                    flag=flag+1
                                    if(flag==a):
                                        myout2=name+'\t'+motif+'/'+motif
                                        with open (out,"a+") as f2:
                                            f2.write(myout2+"\n")
                                            if(os.path.isfile(myout)):
                                                os.remove(myout)
                                            return (1)
                                    str2=str2+tmp2[3]
                                elif(tmp2[7].split(';')[0]=='INDEL' and (tmp2[9].split(':')[0]=='0/1') or (tmp2[9].split(':')[0]=='1/1' and (len(tmp2[4])!=1 or len(tmp2[3])!=1 ))):
                                    indelflag=1
                                    if(tmp2[9].split(':')[0]=='1/1'):
                                        flag11=1
                                    if(len(tmp2[3])<len(tmp2[4])):
                                        str2=str2+tmp[4]
                                    elif(len(tmp2[3])>len(tmp2[4])):
                                        myflag=len(tmp2[3])-len(tmp2[4])
                                        str2=str2+tmp[4]
                                elif(myflag>0):
                                    myflag=myflag-1
                        if( flag11!=1 and indelflag==1):
                            str2=find_maxlen(str2,ssr)
                            myout2=name+'\t'+motif+'/'+str2
                            with open (out,"a+") as f2:
                                f2.write(myout2+"\n")
                                if(os.path.isfile(myout)):
                                    os.remove(myout)
                                return (1)
                        elif(flag11==1):
                            if(len(str2)==a or len(str2)==len(motif) or str2.count('.')):
                                os.remove(myout)
                                return (0)
                            str2=find_maxlen(str2,ssr)
                            myout2=name+'\t'+str2+'/'+str2
                            with open (out,"a+") as f2:
                                f2.write(myout2+"\n")
                                if(os.path.isfile(myout)):
                                    os.remove(myout)
                                return (1)
                    if(os.path.isfile(myout)):
                        os.remove(myout)
                    return(0)


    
def main():
    df=pd.read_csv(x,sep='\t',header=None)
    len_df=len(df)
    thread=cfg.getint("parameters","THREADS")
    if(len_df<thread):
        thread=len_df
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
    df1 = pd.read_csv(x,header=None,sep='\t')
    df2 = pd.read_csv(out,header=None,sep='\t')    
    df = pd.merge(df1,df2,on=0,how="left")
    df=df.drop_duplicates()
    df.to_csv(x,sep='\t',header=False,index=False)
    for a in order:
        os.remove(a)
if __name__ == '__main__' :
    main()
