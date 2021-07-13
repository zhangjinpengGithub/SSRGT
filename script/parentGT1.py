#!/usr/bin/evn python3
import sys,os
import subprocess
import shutil
import re
import pandas as pd
from multiprocessing import Pool
from configparser import ConfigParser
cfg = ConfigParser()
cfg.read("parameters.ini")

x = sys.argv[1] #ref.SSRs
if (sys.argv[1].startswith('--help') or sys.argv[1].startswith('-h')):
        print ('''
###########################################################################################################
Description: This script is used to generate all SSR genotypes of a species by SSR loci information of the reference sequences.
Contact: Chunfa Tong  <tongchf@njfu.edu.cn>
Version: 1.0
Usage:python parentGT1.py ref.SSRs species.sorted.bam  heterozygotes.out all_SSR_types.txt
Parameter:
'ref.SSRs' file is SSR locus information derived from the reference sequence of the SSRMMD software.
'species.sorted.bam' file is generated by SAMtools software.
'heterozygotes.out' and 'all_SSR_types.txt'file are output files.
###########################################################################################################''')

        sys.exit()
y=  sys.argv[2] #bam
out=sys.argv[3]#out 01
outall=sys.argv[4]#all_type

samtools=cfg.get("folders","SAMTOOLS")
thread=cfg.getint("parameter","THREADS")
Depth_Of_Coverage=cfg.getint("parameter","Depth_Of_Coverage")
homozygotes=cfg.get("parameter","Frequency_of_Homozygotes")
heterozygotes=cfg.get("parameter","Frequency_of_Heterozygotes")
homozygotes=float(homozygotes)
heterozygotes=float(heterozygotes)
major_allele=heterozygotes

with open(x,'r') as f:
    next(f)
    for line in f:
        tmp=line.split("\t")
        genename=tmp[1]+':'+str(int(tmp[6])-3)
        motif=tmp[2]+'_'+tmp[4]
        myout=genename+'\t'+motif
        with open('my_id',"a+") as f2:
            f2.write(myout+"\n")

def run_command(cmd):
        print(cmd)
        return_code = subprocess.call(cmd, shell=True)

def find_maxlen(s,a):
        c=re.finditer(r'(%s)+'%(a),s)
        find_maxlist=[]
        for i in c:
                find_maxlist.append(i.group())
 #       print(find_maxlist) ###
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
    ref_end=ref_end ## #+3
###    b=a
    a=a[ref_start:ref_end]
    len_a=ref_end-ref_start
    if (" " not in a and len(re.findall(r"[A,T,C,G,*]",a))==len_a and a.count(motif_str)!=0 and a.count(',')==0):
        a=a.replace("*","") ###b.rep->a.rep
        print (a)
        str2=find_maxlen(a,motif_str)
#example: a=ATATCTATATAT motif_str=AT  return ATATCTATATAT

        return (str2)
    else:
        return(0)
        

def GMSSR(gene,c):
    motif_str=''
    motif_repeat=0
    motif=''
    motif_str=c.split("_")[0]
    motif_repeat=int(c.split("_")[1])
    motif=motif_str*motif_repeat
    SSRprint(motif_str,gene,motif_repeat,motif)
def SSRprint(motif_str,gene,motif_repeat,motif):
    order=[]
    len1=len(motif) #8
    with open(gene,"r") as f:
        lines = f.readlines()[1:]
        for line in lines:
            refstr=line.strip()
            break
    ref_start=refstr.find("N",0)
    ref_end=ref_start+len1+refstr.count("*",ref_start)+3  #+4 
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
                dict=sorted(dict.items(),key=lambda x:x[1],reverse=True)
                #order by value , reverse=True is from  largest to smallest
                my_dtct=str(dict).strip('{}[]').replace("\'", "").replace("(", "").replace(")", "")
                print(my_dtct) # 4, 27, 5, 20, 3, 1  #AT(4) 27 ,AT(5) 20
           #     print(re.findall(r"[A,T,C,G]",my_dtct))
#TCTCTCTCTCTCTCTCTCTCTC, 8, TCTCTCTCTCTCTCGCTCTCTCTC, 6, TCTCTCGCTCTCTCTC, 2, TCTCTCTCTCTCTCTCTCTC, 1
                gene_rep=gene.replace(".txt", "")
                if (len(dict)!=0):
                    first=int(my_dtct.split(",")[1])
                    mydtct0=my_dtct.split(",")[0]
                    b=re.findall(r'\d+',my_dtct)
                    order_len=0
                    for i in range(len(b)):
                        order_len=int(b[i])+ order_len
                    print(order_len)
                    order_len=float(order_len)
                    if (my_dtct.count(",")==1  and first>4 and len(my_dtct.split(",")[0])>7):
                        myout=gene_rep+'\t'+motif_str+'('+str(motif_repeat)+')'+'\t'+str(len1)+'\t'+mydtct0+'/'+mydtct0+'\t'+'aa'
                        with open(outall,"a+") as f3:
                            f3.write(myout+"\n")
                    elif(my_dtct.count(",")>1 and (len(my_dtct.split(",")[0])>7 or len(my_dtct.split(",")[2])>7)):
                        first=float(my_dtct.split(",")[1])
                        second=float(my_dtct.split(",")[3])
                        mydtct2=my_dtct.split(",")[2].lstrip()
                        if(first/order_len>=major_allele and second/order_len>=heterozygotes and int(first)>2 and int(second)>2 and first>second ):
                            myout=gene_rep+'\t'+motif_str+'('+str(motif_repeat)+')'+'\t'+str(len1)+'\t'+mydtct0+'/'+mydtct2+'\t'+'ab'
                            with open(out,"a+") as f2:
                                f2.write(myout+"\n")
                            with open(outall,"a+") as f3:
                                f3.write(myout+"\n")
                        elif(first/order_len>homozygotes and int(first)>4 and len(my_dtct.split(",")[0])>7):    
                            myout=gene_rep+'\t'+motif_str+'('+str(motif_repeat)+')'+'\t'+str(len1)+'\t'+mydtct0+'/'+mydtct0+'\t'+'aa'
                            with open(outall,"a+") as f3:
                                f3.write(myout+"\n")
                    

def my_samtoools(gene,c):
    global samtools
    geneout=gene+'.txt'
    cmd=samtools+' tview -p '+gene+' -d T '+y+'>'+geneout
    run_command(cmd)
    GMSSR(geneout,c)
    os.remove(geneout)
def callSSR(i):
    file=open(i,"r")
    lines = file.readlines()
    for line in lines:
        line=line.strip()
        tmp=line.split("\t")
        gene=tmp[0]
        my_samtoools(gene,tmp[1])
    file.close
def main():
    
    x='my_id'
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
    os.remove('my_id')	
    
if __name__ == '__main__' :
    main()
