#!/usr/bin/evn python3
import sys,os
import subprocess
import shutil
import re
import pandas as pd
from multiprocessing import Pool
from configparser import ConfigParser
x = sys.argv[1] #Male_marker.txt
y=  sys.argv[2] #progan1_filter.sorted.bam
out=sys.argv[3] #
cfg = ConfigParser()
cfg.read("parameters.ini")
samtools=cfg.get("folders","SAMTOOLS")
thread=cfg.getint("parameter","THREADS")
Depth_Of_Coverage=cfg.getint("parameter","Depth_Of_Coverage")
homozygotes=cfg.get("parameter","Frequency_Of_Homozygotes")
heterozygotes=cfg.get("parameter","Frequency_Of_Heterozygotes")
homozygotes=float(homozygotes)
heterozygotes=float(heterozygotes)
major_allele=heterozygotes
def run_command(cmd):
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
def find_maxlen(s,a):
        c=re.finditer(r'(%s)+'%(a),s)
        find_maxlist=[]
        for i in c:
                find_maxlist.append(i.group())
        if (len(find_maxlist)==1 ):
            str2=max(find_maxlist,key=len, default='')
            return (str2)
        elif (len(find_maxlist)==2):
            if (len(find_maxlist[0])>3 and len(find_maxlist[1])>3):
                start=int(s.find(a))
                end=int(s.rfind(a))+len(a)
                return(s[start:end])
            else:
                str2=max(find_maxlist,key=len, default='')
                return (str2)
        else:
            return (0)

def ssr(line,ref_end,motif_str):
    a=line.upper()
    ref_start=0
    ref_end=ref_end
    a=a[ref_start:ref_end]
    len_a=ref_end-ref_start
    if (" " not in a and len(re.findall(r"[A,T,C,G,*]",a))==len_a and a.count(motif_str)!=0 and a.count(',')==0):
#    if " " not in a:
        a=a.replace("*","")
        str2=find_maxlen(a,motif_str)#example: a=ATATCTATATAT motif_str=AT  return ATATAT
        return (str2)
    else:
        return(0)
def GMSSR(gene,motif,motif2,refmotif):
    motif_str=''
    motif_repeat=0
    order=[]
    re_search=re.findall('\w+',refmotif)
    motif_str=re_search[0]
    motif_repeat=int(re_search[1])
    refmotif=motif_str*motif_repeat
    len1=len(refmotif)
    with open(gene,"r") as f:
        next(f)
        lines = f.readlines()
        for line in lines:
            refstr=line
            break
    ref_start=refstr.find("N",0)
    ref_end=ref_start+len1+refstr.count("*",ref_start)+3
    with open(gene,"r") as f:
        lines = f.readlines()[3:]
        flen=len(lines)
        if (flen < Depth_Of_Coverage):
            print("wrong")
        else:
            for line in lines:
                if (line.startswith(" ")):
                    continue
                else:
                    str_numble=ssr(line,ref_end,motif_str)
                    if(str_numble!=0 and str_numble!='None'):
                        order.append(str(str_numble))
            order_len=len(order)
            if (order_len==0):
                print("wrong")
            else:
                myset = set(order)
                dict={}
                for item in myset:
                    dict.update({item:order.count(item)})
                dict={k:v for k, v in dict.items() if v>=3}
                dict=sorted(dict.items(),key=lambda x:x[1],reverse=True)#order by value , reverse=True is from  largest to smallest
                my_dtct=str(dict).strip('{}[]').replace("\'", "").replace("(", "").replace(")", "")
                print(my_dtct)
                if (len(dict)==0):
                    return (0)
                else:
                    b=re.findall(r'\d+',my_dtct)
                    order_len=0
                    for i in range(len(b)):
                        order_len=int(b[i])+ order_len
                print(order_len)
                order_len=float(order_len)
                mydtct0=my_dtct.split(",")[0]
                gene_rep=gene.replace(".txt", "")
                first=int(my_dtct.split(",")[1])
                if (my_dtct.count(",")==1  and first>4 and len(my_dtct.split(",")[0])>7):
                    myout=gene_rep+'\t'+mydtct0+'/'+mydtct0
                    with open(out,"a+") as f2:
                        f2.write(myout+"\n")
                elif(my_dtct.count(",")>1 and (len(mydtct0)>7 or len(my_dtct.split(",")[2])>7)):
                    first=float(my_dtct.split(",")[1])#27.0
                    second=float(my_dtct.split(",")[3])#20.0
                    mydtct2=my_dtct.split(",")[2].lstrip()
                    if(first/order_len>homozygotes  and int(first)>4 and len(my_dtct.split(",")[0])>7):
                        myout=gene_rep+'\t'+mydtct0+'/'+mydtct0
                        with open(out,"a+") as f2:
                            f2.write(myout+"\n")
                    elif(first/order_len>=major_allele and second/order_len>=heterozygotes and int(first)>2 and int(second)>2 ):
                        myout=gene_rep+'\t'+mydtct0+'/'+mydtct2
                        with open(out,"a+") as f2:
                            f2.write(myout+"\n")                    
def my_samtoools(gene,motif2,motif3,motif):
    global samtools
    geneout=gene+'.txt'
    cmd=samtools+' tview -p '+gene+' -d T '+y+'>'+geneout
    run_command(cmd)
    GMSSR(geneout,motif2,motif3,motif)
    os.remove(geneout)
def callSSR(i):
    file=open(i,"r")
    lines = file.readlines()
    for line in lines:
        tmp=line.split("\t")
        gene=tmp[0]
        motif=tmp[1]
        motif2=tmp[3].split("/")[0]	 
        motif3=tmp[3].split("/")[1]
        my_samtoools(gene,motif2,motif3,motif)
    file.close
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
    df1 = pd.read_table(x,header=None)
    df2 = pd.read_table(out,header=None)    
    df = pd.merge(df1,df2,on=0,how="left")
    df=df.drop_duplicates()
    df.to_csv(x,sep='\t',header=False,index=False)
    for a in order:
        os.remove(a)
if __name__ == '__main__' :
    main()
