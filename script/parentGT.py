#!/usr/bin/evn python3
import sys,os
import subprocess
import shutil
import re
import pandas as pd
from configparser import ConfigParser
cfg = ConfigParser()
cfg.read("parameters.ini")
x = sys.argv[1]#vcf
out=sys.argv[2]#out 01
outall=sys.argv[3]#all_type
thread=cfg.get("parameter","THREADS")
DP=cfg.get("parameter","Depth Of Coverage")
Alleles_Quality=cfg.get("parameter","Alleles quality score")
DP=int(DP)
Alleles_Quality=float(Alleles_Quality)
def call(x):
    with open(x,'r') as f:
        for line in f:
            if (line.startswith("#")):
                continue
            else:
                tmp=line.strip().split("\t")
                if(tmp[9]=='.'):
                    continue
                elif(tmp[9].find('/')):
                    GT=tmp[9].split(":")[0]
                    myDP=int(tmp[9].split(":")[1])
                    Q=float(tmp[9].split(":")[2])
                    if (myDP>=DP and Q>=Alleles_Quality):
                        mydtct0=""
                        mydtct2=""
                        gene_rep=tmp[0]+":"+tmp[1]
                        motif_str=tmp[7].split(";")[1].split("=")[1]
                        motif_repeat=tmp[7].split(";")[3].split("=")[1]
                        len1=len(tmp[3])
                        if(GT=='0/0'):
                            mydtct0=tmp[3]
                            myout=gene_rep+'\t'+motif_str+'('+str(motif_repeat)+')'+'\t'+str(len1)+'\t'+mydtct0+'/'+mydtct0+'\t'+'aa'
                            with open(outall,"a+") as f3:
                                f3.write(myout+"\n")
                        elif(GT=='1/1'):
                            mydtct0=tmp[4]
                            if(len(mydtct0)>5):
                                myout=gene_rep+'\t'+motif_str+'('+str(motif_repeat)+')'+'\t'+str(len1)+'\t'+mydtct0+'/'+mydtct0+'\t'+'aa'
                                with open(outall,"a+") as f3:
                                    f3.write(myout+"\n")
                        elif(GT=='0/1' or GT=='1/0'):
                            mydtct0=tmp[3]
                            mydtct2=tmp[4]
                            if(len(mydtct2)>5):
                                myout=gene_rep+'\t'+motif_str+'('+str(motif_repeat)+')'+'\t'+str(len1)+'\t'+mydtct0+'/'+mydtct2+'\t'+'ab'
                                with open(out,"a+") as f2:
                                    f2.write(myout+"\n")
                                with open(outall,"a+") as f3:
                                    f3.write(myout+"\n")
                        elif(GT=='1/2'):
                            mydtct0=tmp[4].split(",")[0]
                            mydtct2=tmp[4].split(",")[1]
                            if(len(mydtct2)>5 and len(mydtct0)>5):
                                myout=gene_rep+'\t'+motif_str+'('+str(motif_repeat)+')'+'\t'+str(len1)+'\t'+mydtct0+'/'+mydtct2+'\t'+'ab'
                                with open(out,"a+") as f2:
                                    f2.write(myout+"\n")
                                with open(outall,"a+") as f3:
                                    f3.write(myout+"\n")
                    else:
                        continue
call(x)
