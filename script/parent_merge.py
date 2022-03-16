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
ref=sys.argv[4]
cfg.read("parameters.ini")
thread=cfg.getint("parameters","THREADS")
bcftools=cfg.get("folders","BCFTOOLS_FOLD")+'/bcftools'
DP=cfg.get("parameters","ALLELE_DEPTH")
DP=int(DP)
GQ=cfg.get("parameters","GQ")
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
            ssr=tmp[1].split("(")[0] #CA
            motif=a.split("/")[0] #CACACAC
            motif2=a.split("/")[1] #CACA
            GMSSR(myout,tmp[0],ssr,motif,motif2) 
            os.remove(myout)
def find_maxlen(s,a):
        s=s.upper()
        c=re.finditer(r'(%s)+'%(a),s)
        find_maxlist=[]
        for i in c:
            find_maxlist.append(i.group())
        str2=max(find_maxlist,key=len, default='')
        if(len(str2)<6):
            return ('0')
        return str2
def GMSSR(geneout,name,ssr,motif,motif2):
    hangnum=abs(int(geneout.split("-")[0].split(":")[1])-int(geneout.split("-")[1].split(".")[0]))
    flag=0
    flag1=0
    a=int(geneout.split('-')[1].replace(".txt",""))-int(geneout.split('-')[0].split(':')[1])+1
    a1=name.split(":")[1] #chr1:position
    stra=""
    with open(geneout,"r") as f1:
        f1=f1.readlines()
        for line in f1:
            line=line.strip()
            if line.startswith("#"):
                continue
            else:
                tmp=line.split("\t")
                if(tmp[7].split(';')[0]=='INDEL' and str(int(tmp[1])+1)==a1 ):
                    if( int(tmp[7].split(';')[3].split('=')[1])>=DP and tmp[3].find(ssr)):
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
                            elif(motif==mydtct0 and motif2!=mydtct2!='0' and mydtct2.find(ssr)): #
                                myout=name+'\t'+mydtct0+'/'+mydtct2+'\t'+'ac'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='1/1' ): 
                            mydtct2=find_maxlen(tmp[4],ssr)
                            if(motif2==mydtct2 and mydtct2!='0'):
                                myout=name+'\t'+mydtct2+'/'+mydtct2+'\t'+'aa'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                            elif(motif2!=mydtct2 and mydtct2.count(ssr) and len(mydtct2)>5):  ###
                                myout=name+'\t'+mydtct2+'/'+mydtct2+'\t'+'cc'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                        elif(tmp[9].split(':')[0]=='1/2' ):
                            mydtct0=find_maxlen(tmp[4].split(",")[0],ssr)
                            mydtct2=find_maxlen(tmp[4].split(",")[1],ssr)
                            if((mydtct0==motif2 or mydtct2==motif2 ) and  mydtct2!='0'):
                                myout=name+'\t'+mydtct0+'/'+mydtct2+'\t'+'ac'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                            elif((mydtct0!=motif2 or mydtct2!=motif2) and mydtct0.count(ssr) and mydtct2.count(ssr) and len(mydtct2)>5):
                                myout=name+'\t'+mydtct0+'/'+mydtct2+'\t'+'cd'
                                with open (out,"a+") as f2:
                                    f2.write(myout+"\n")
                                return (1)
                    else:
                        return (0)
                elif(tmp[7].split(';')[0]!='INDEL' and tmp[9].split(':')[0]=='0/0' ):
                    str2=""
                    myflag=0
                    flag11=0
                    hang=0
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
                                hang=hang+1
                                if(tmp2[9].split(':')[0]=='0/0'):
                                    flag=flag+1
                                    if(flag==a):
                                        myout2=name+'\t'+motif+'/'+motif+'\t'+'aa'
                                        with open (out,"a+") as f2:
                                            f2.write(myout2+"\n")
                                            os.remove(myout)
                                            return (1)
                                    str2=str2+tmp2[3]
                                elif(tmp2[7].split(';')[0]=='INDEL' and (tmp2[9].split(':')[0]=='0/1') or (tmp2[9].split(':')[0]=='1/1')):
                                    if(tmp2[9].split(':')[0]=='1/1'):
                                        flag11=1
                                    if(len(tmp2[3])<len(tmp2[4])):
                                        str2=str2+tmp[4]
                                    elif(len(tmp2[3])>len(tmp2[4])):
                                        myflag=len(tmp2[3])-len(tmp2[4])
                                        str2=str2+tmp[4]
                                elif(myflag>0):
                                    myflag=myflag-1
                        if(hang<hangnum):
                            os.remove(myout)
                            return (0)
                        if(hang>hangnum):
                            if(abs(len(str2)-len(motif))==1):
                                os.remove(myout)
                                return (0)
                        str2=find_maxlen(str2,ssr)
                        if(str2=="0"):
                            os.remove(myout)
                            return (0)
                        if(str2==motif2 and flag11!=1):
                            myout2=name+'\t'+motif+'/'+motif2+'\t'+'ab'
                            with open (out,"a+") as f2:
                                f2.write(myout2+"\n")
                                os.remove(myout)
                                return (1)
                        elif(str2!=motif2 and flag11!=1):
                            myout2=name+'\t'+motif+'/'+str2+'\t'+'ac'
                            with open (out,"a+") as f2:
                                f2.write(myout2+"\n")
                                os.remove(myout)
                                return (1)
                        elif(flag11==1):
                            if(str2==motif2):
                                myout2=name+'\t'+motif2+'/'+motif2+'\t'+'aa'
                                with open (out,"a+") as f2:
                                    f2.write(myout2+"\n")
                                    os.remove(myout)
                                    return (1)
                            elif(str2!=motif2):
                                myout2=name+'\t'+str2+'/'+str2+'\t'+'cc'
                                with open (out,"a+") as f2:
                                    f2.write(myout2+"\n")
                                    os.remove(myout)
                                    return (1)
                    if(os.path.isfile(myout)):
                        os.remove(myout)
                    return(0)
                elif(tmp[7].split(';')[0]!='INDEL' and tmp[9].split(':')[0]=='1/1' ):
                    return(0)

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

