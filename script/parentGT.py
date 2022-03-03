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
ref=sys.argv[4]
cfg.read("parameters.ini")
thread=cfg.getint("parameters","THREADS")
bcftools=cfg.get("folders","BCFTOOLS_FOLD")+'/bcftools'
DP=cfg.get("parameters","ALLELE_DEPTH")
DP=int(DP)
GQ=cfg.get("parameters","GQ")
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
                chrname=tmp[1]+':'+tmp[6]+'-'+str(int(tmp[7]))
                myout=chrname+'.txt'
                cmd=bcftools+' norm -f '+ref+' -m +both -r '+chrname+' '+y+' 1 >'+myout +" 2>log"
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
        if(len(str2)<6):
            return ('0')
        return str2

def GMSSR(geneout,ssr,motif):
    str2=""
    a=""
    indelflag=0
    hangnum=abs(int(geneout.split("-")[0].split(":")[1])-int(geneout.split("-")[1].split(".")[0]))
    with open(geneout,"r") as f1:
        f1=f1.readlines()
        for line in f1:
            line=line.strip()
            if line.startswith("#"):
                continue
            else:
                tmp=line.split("\t")
                if(tmp[7].split(';')[0]=='INDEL' and tmp[1]==str(int(geneout.split("-")[0].split(":")[1])-1) and (tmp[9].split(":")[0]=='0/1' or tmp[9].split(":")[0]=='1/2')):
                    if (tmp[9].split(":")[0]=='0/1' ):
                        if(int(tmp[9].split(":")[-1])>=int(GQ)):
                            name=tmp[0]+":"+str(int(tmp[1])+1)
                            if (name in mydict.keys()):
                                ssr=mydict[name].split("_")[0] #AT
                                motif=mydict[name].split("_")[1]  #5
                                str1=ssr*int(motif)
                                if(len(tmp[3])<len(tmp[4])): #gain
                                    for line2 in f1:
                                        line2=line2.strip()
                                        if line2.startswith("#"):
                                            continue
                                        else:
                                            tmp2=line2.split("\t")
                                            if(tmp2[4]=="."):
                                                str2=str2+tmp2[3]
                                            else:
                                                str2=str2+tmp2[4]
                                    str2=find_maxlen(str2,ssr)
                                    if (str1!=str2 and str2!="" and str2!="0"):
                                        myout=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(str1))+'\t'+str1+'/'+str2+'\t'+'ab'
                                        with open(myfile,"a+") as f2:
                                            f2.write(myout+"\n")
                                            return (1)
                                elif(len(tmp[3])>len(tmp[4])): #loss
                                    num=len(tmp[3])-len(tmp[4])+1
                                    for line2 in f1:
                                        line2=line2.strip()
                                        if line2.startswith("#"):
                                            continue
                                        else:
                                            tmp2=line2.split("\t")
                                            if(tmp2[4]=="." and num==0):
                                                str2=str2+tmp2[3]
                                            elif(tmp2[4]!="." and num==0):
                                                str2=str2+tmp2[4]
                                            if(num!=0):
                                                num=num-1
                                    str2=find_maxlen(str2,ssr)
                                    if (str1!=str2 and str2!="" and str2!="0"):
                                        myout=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(str1))+'\t'+str1+'/'+str2+'\t'+'ab'
                                        with open(myfile,"a+") as f2:
                                            f2.write(myout+"\n")
                                            return (1)
                    elif(tmp[9].split(":")[0]=='1/2'):
                        if(len(tmp[4])<7):
                            return (0)
                        myout=geneout.replace('.txt','')
                        cmd=bcftools+' view -r '+myout+' '+y+' >'+myout
                        run_command(cmd)
                        with open(myout,"r") as f3:
                            f3=f3.readlines()
                            for line3 in f3:
                                line3=line3.strip()
                                if line3.startswith("#"):
                                    continue
                                else:
                                    tmp3=line3.split("\t")
                                    if(int(tmp3[9].split(":")[-1])>=int(GQ)):
                                        name=tmp3[0]+":"+str(int(tmp3[1])+1)
                                        if (name in mydict.keys() and tmp3[9].split(":")[0]=='1/2'): ###tmp3
                                            #####
                                            ssr=mydict[name].split("_")[0]
                                            motif=mydict[name].split("_")[1]
                                            str1=find_maxlen(tmp3[4].split(",")[0],ssr)
                                            str2=find_maxlen(tmp3[4].split(",")[1],ssr)
                                            if(str1.count(ssr) and str2.count(ssr) and  str2!="0"):
                                                myout2=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(ssr*int(motif)))+'\t'+str1+'/'+str2+'\t'+'ab'
                                                with open(myfile,"a+") as f2:
                                                    f2.write(myout2+"\n")
                                                    os.remove(myout)
                                                    return (1)
                                        elif(int(geneout.split("-")[0].split(":")[1])-int(tmp3[1])>1):  ###
                                            pass
                                        else:
                                            os.remove(myout)
                                            return (0)
                        os.remove(myout)
                        return (0)



                elif((tmp[7].split(';')[0]=='INDEL' and tmp[1]!=str(int(geneout.split("-")[0].split(":")[1])-1) and (tmp[9].split(":")[0]=='0/1' ))):
                    num=len(tmp[3])-len(tmp[4])
                    num2=int(geneout.split("-")[0].split(":")[1])-int(tmp[1])
                    name=geneout.split("-")[0]
                    str1=ssr*int(motif)
                    if(num>0):
                        for line2 in f1:
                            line2=line2.strip()
                            if line2.startswith("#"):
                                continue
                            else:
                                tmp2=line2.split("\t")
                                if(tmp2[4]=="."):
                                    str2=str2+tmp2[3]
                                else:
                                    str2=str2+tmp2[4]
                                    str2=str2[num2:]
                                    str2=tmp[4]+str2
                                    str2=find_maxlen(str2,ssr)
                                    if (str1!=str2 and str2!=""  and str2!="0"):
                                        myout=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(str1))+'\t'+str1+'/'+str2+'\t'+'ab'
                                        with open(myfile,"a+") as f2:
                                            f2.write(myout+"\n")
                                            return (1)
                else:
                    str1=ssr*int(motif)
                    if( tmp[7].split(';')[0]!='INDEL' and int(tmp[9].split(':')[-1])>=DP and tmp[9].split(":")[0]=='0/0'):
                        flag=0
                        hang=0
                        indelflag=0
                        for line2 in f1:
                            line2=line2.strip()
                            if line2.startswith("#"):
                                continue
                            else:
                                tmp2=line2.split("\t")
                                hang=hang+1
                                if(tmp2[7].split(';')[0]=='INDEL' and  tmp[9].split(":")[0]!='0/0'):  ###
                                    indelflag=1
                                if(tmp2[4]=="." and tmp2[9].split(":")[0]=='0/0' and flag==0):
                                    str2=str2+tmp2[3]
                                elif(tmp2[7].split(';')[0]=='INDEL' and tmp2[9].split(":")[0]=='0/1'):
                                    if(len(tmp2[3])<len(tmp2[4])):
                                        str2=str2+tmp[4]
                                    elif(len(tmp2[3])>len(tmp2[4])):
                                        flag=len(tmp2[3])-len(tmp2[4])
                                        str2=str2+tmp[4] 
                                elif(flag>0):
                                    flag=flag-1
                        if(indelflag==0):
                            return (0)
                        if(hang<hangnum):
                            return (0)
                        if(str2=="0"):
                            return (0)
                        str2=find_maxlen(str2,ssr)
                        if (str1!=str2 and str2!="" and str2!="0"):
                            name=geneout.split("-")[0]
                            myout=name+'\t'+ssr+'('+str(motif)+')'+'\t'+str(len(str1))+'\t'+str1+'/'+str2+'\t'+'ab'
                            with open(myfile,"a+") as f2:
                                f2.write(myout+"\n")
                                return (1)      
                        return (0)  
                    elif(tmp[7].split(';')[0]=='INDEL' and tmp[9].split(":")[0]=='1/1'):
                        return (0)

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

