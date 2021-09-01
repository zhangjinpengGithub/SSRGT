#!/usr/bin/evn python3
import sys,os,argparse
import subprocess
import shutil
import re
import pandas as pd
from multiprocessing import Pool
from configparser import ConfigParser
cfg = ConfigParser()
cfg.read("parameters.ini")
parser = argparse.ArgumentParser(
     description='''
Example: python parent_merge.py -parent1 male01.txt -parent2 femaleallSSR_type.txt -o female.out
''')
parser.add_argument('-parent1','--Parent1',required=True)
parser.add_argument('-parent2','--Parent2',required=True)
parser.add_argument('-o','--out',required=True)

args = parser.parse_args()

x = args.Parent1 #male01.txt
y=  args.Parent2 #femaleallSSR_type.txt
out=args.out #out 
thread=cfg.getint("parameter","THREADS")
def run_command(cmd):
        print(cmd)
        return_code = subprocess.call(cmd, shell=True)

def GMSSR(gene,motif,motif2,refmotif,mydtct0,mydtct2):
    gene_rep=gene.replace(".txt", "")
    if((mydtct0==mydtct2) and (mydtct0==motif or mydtct0==motif2)):
        myout=gene_rep+'\t'+mydtct0+'/'+mydtct0+'\t'+'aa'
        with open(out,"a+") as f2:
            f2.write(myout+"\n")
    elif(mydtct0!=motif and mydtct0!=motif2 and (mydtct0==mydtct2)):
        myout=gene_rep+'\t'+mydtct0+'/'+mydtct0+'\t'+'cc'
        with open(out,"a+") as f2:
            f2.write(myout+"\n")
    elif(((mydtct0==motif and mydtct2==motif2 )or (mydtct2==motif and mydtct0==motif2))and (mydtct0!=mydtct2)):
        myout=gene_rep+'\t'+mydtct0+'/'+mydtct2+'\t'+'ab'
        with open(out,"a+") as f2:
            f2.write(myout+"\n")
    elif((mydtct0!=mydtct2) and ((mydtct0==motif and mydtct2!=motif2)or(mydtct2==motif and mydtct0!=motif2)or (mydtct0==motif2 and mydtct2!=motif) or (mydtct2==motif2 and mydtct0!=motif))):
        myout=gene_rep+'\t'+mydtct0+'/'+mydtct2+'\t'+'ac'
        with open(out,"a+") as f2:
            f2.write(myout+"\n")
    elif(mydtct0!=motif!= mydtct2!=motif2):
        myout=gene_rep+'\t'+mydtct0+'/'+mydtct2+'\t'+'cd'
        with open(out,"a+") as f2:
            f2.write(myout+"\n")
def my_samtoools(gene,motif,motif2,refmotif):
#    global samtools
    geneout=gene+'.txt'
    cmd=' grep  '+gene+' '+y+'>'+geneout
    run_command(cmd)
    with open(geneout,"r") as f2:
        lines=f2.readlines()
        flen=len(lines)
        if (flen < 1):
            print("wrong")
        elif(flen >= 1):
            for line in lines:
                tmp=line.strip().split('\t')
                a=tmp[3].split("/")[0] 
                b=tmp[3].split("/")[1]
                GMSSR(geneout,motif,motif2,refmotif,a,b)
    os.remove(geneout)
def callSSR(i):
    file=open(i,"r")
    lines = file.readlines()
    for line in lines:
        tmp=line.split("\t")
        gene=tmp[0]
        refmotif=tmp[1]
        motif=tmp[3].split("/")[0]
        motif2=tmp[3].split("/")[1]	 #AT(10)
        my_samtoools(gene,motif,motif2,refmotif)
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
    for a in order:
        os.remove(a)

if __name__ == '__main__' :
    main()

