#!/usr/bin/evn python3
import sys,os
import subprocess
import argparse
import shutil
import re
import pandas as pd
from scipy import stats
from multiprocessing import Pool
import time
import os.path
from configparser import ConfigParser
cfg = ConfigParser()
cfg.read("parameters.ini")
bwa=cfg.get("folders","BWA_FOLD")
samtools=cfg.get("folders","SAMTOOLS_FOLD")+'/samtools'
bcftools=cfg.get("folders","BCFTOOLS_FOLD")+'/bcftools'
thread=cfg.getint("parameters","THREADS")
GQ=cfg.get("parameters","GQ")
pvalue=cfg.get("parameters","PVALUE")
MISPCT=cfg.get("parameters","MISS_GENOTYPES")
data_fold=cfg.get("folders","RADDATA_FOLD")
GENOME=data_fold+'/'+cfg.get("data files","REFERENCE_FILE")
progenyfold=data_fold
parentfold=data_fold
script=cfg.get("folders","SSRGT_FOLD")+'/script'
SSRMMD=script+'/SSRMMD.pl'
DP=cfg.get("parameters","ALLELE_DEPTH")
DP=int(DP)
###########################################################################################################
Home=os.environ['HOME']
if (bwa.startswith("~")):
	bwa=bwa.replace("~",Home)
if (samtools.startswith("~")):
        samtools=samtools.replace("~",Home)
if (data_fold.startswith("~")):
        data_fold=data_fold.replace("~",Home)
if (script.startswith("~")):
        script=script.replace("~",Home)
files =os.listdir(bwa)
for file in files:
	if (file.endswith('bwa-mem2')):
		bwa=bwa+'/bwa-mem2'
	elif(file.endswith('bwa')):
		bwa=bwa+'/bwa'
def run_command(cmd):
#	print(cmd)
	return_code = subprocess.call(cmd, shell=True)
if(os.path.isfile(GENOME)):
        cmd1 = 'cp '+GENOME+ ' ./'
        run_command(cmd1)
        GENOME='./'+cfg.get("data files","REFERENCE_FILE")
def get_SSR(reader1):
	cmd1 = 'perl '  +SSRMMD+' '+' -f1 '+reader1+' -l 100 -ss 1  -mo '+ motif +' -t '+str(thread)+' -o ./'
	run_command(cmd1)
	cmd2 = samtools + ' faidx ' + reader1  
	run_command(cmd2)
	files = os.listdir()
	for file in files:
		if (file.find('SSRs')!=-1):
        		SSR=file
	minmotif=motif.split(',')[0].split('=')[1]
	file = open(SSR, "r")
	file1 = open('SSR.txt', "w")
	file1.close()
	lines = file.readlines()
	for line in lines:
		line=line.strip()
		if line.startswith("number"):
			continue
		else:
			tmp=line.split("\t")
			if (tmp[3]=='1' and int(tmp[4])<int(minmotif)):
				continue
			else:
				with open('SSR.txt',"a+") as f2:
					f2.write(line+"\n")
	file.close()
def find_parent(reader1):
        	
        cmd = bwa + ' index '	+reader1
        run_command(cmd)
        cmd = bwa +  ' mem  -M  '+reader1+' '+male1+' '+male2+' > male.sam'
        run_command(cmd)
        cmd=samtools+' sort  '+' male.sam > male.bam'
        run_command(cmd)
        cmd=samtools+' index  male.bam'
        run_command(cmd)
        cmd=bcftools +" mpileup  -Obuzv -a AD,INFO/AD -f "+reader1+" male.bam  | "+bcftools +" call -m -f gq -Oz -o male.vcf.gz"
        run_command(cmd)
        cmd=bcftools+" index male.vcf.gz" 
        run_command(cmd)
        cmd = bwa +  ' mem  -M  '+reader1+' '+female1+' '+female2+' > female.sam'
        run_command(cmd)
        cmd=samtools+' sort   '+' female.sam >female.bam'
        run_command(cmd)
        cmd=samtools+' index  female.bam'
        run_command(cmd)
        cmd=bcftools +" mpileup  -Obuzv -a AD,INFO/AD -f "+reader1+" female.bam  | "+bcftools +" call -m -f gq -Oz -o female.vcf.gz"
        run_command(cmd)
        cmd=bcftools+" index female.vcf.gz"
        run_command(cmd)
def get_proganyID(list1,proganyID):
        with open(proganyID) as f:
                for line in f:
                        a=line.strip().split("\t")
                        list1.append(a[0])
def callparent2():
        cmd='python   '+script+'/parentGT.py '+' SSR.txt male.vcf.gz male.01.txt'
        run_command(cmd)
        cmd='python   '+script+'/parentGT.py '+' SSR.txt female.vcf.gz female.01.txt'
        run_command(cmd)
        cmd='python   '+script+'/parent_merge.py '+' female.01.txt male.vcf.gz male.out'
        run_command(cmd)
        cmd='python   '+script+'/parent_merge.py '+' male.01.txt female.vcf.gz female.out'
        run_command(cmd)
        cmd='less -S female.out |sort|uniq > '+'female.sort.out'
        run_command(cmd)
        cmd='less -S male.out |sort|uniq > '+'male.sort.out'
        run_command(cmd)
        df1 = pd.read_csv('female.01.txt',sep='\t',header=None)
        df2 = pd.read_csv('male.sort.out',sep='\t',header=None)
        data = pd.merge(df1,df2,on=0,how="left")
        df=data.dropna(axis=0,how='any')
        df.to_csv('Female_marker.txt',sep='\t',header=False,index=False)
        df1 = pd.read_csv('male.01.txt',sep='\t',header=None)
        df2 = pd.read_csv('female.sort.out',sep='\t',header=None)
        data = pd.merge(df1,df2,on=0,how="left")
        df=data.dropna(axis=0,how='any')
        df.to_csv('Male_marker.txt',sep='\t',header=False,index=False)

def Mapping_Pool():
    df=pd.read_csv(progenyID,sep='\t',header=None)
    len_df=len(df)
    if(len_df<thread):
        mythread=len_df
    else:
        mythread=thread
    numble=int(len_df/mythread)
    gtrd=pd.read_csv(progenyID,chunksize=numble,sep='\t',header=None)
    a=0
    order=[]
    for i in gtrd :
        a=a+1
        order.append(str(a)+'.tmptxt')
        i.to_csv(f"{str(a)}.tmptxt",sep='\t',header=None,index=False)
    p=Pool(mythread)
    p.map(Mapping_progeny,order)
    p.close()
    p.join()
    for a in order:
        os.remove(a)
def Mapping_progeny(i):
        with open(i) as f:
                for line in f:
                        tmp = line.strip().split('\t')
                        cmd = bwa +  ' mem  -M  '+GENOME+' '+progenyfold+'/'+tmp[1]+' '+progenyfold+'/'+tmp[2]+' > '+tmp[0]+'.sam'
                        run_command(cmd)
                        cmd=samtools+' sort   '+tmp[0]+'.sam'+' >'+tmp[0]+'.sort.bam'
                        run_command(cmd)
                        cmd=samtools+' index  ' +tmp[0]+'.sort.bam'
                        run_command(cmd)
                        cmd=bcftools +" mpileup  -Obuzv -a AD,INFO/AD -f "+GENOME+" "+tmp[0]+".sort.bam"+" | "+bcftools +" call -m -f gq -Oz -o "+tmp[0]+".vcf.gz"
                        run_command(cmd)
                        cmd=bcftools+" index "+tmp[0]+".vcf.gz"
                        run_command(cmd)
                        tmpsam=tmp[0]+'.sam'
                        os.remove(tmpsam)

def SSRGM(reader1):
	df=pd.read_csv('Female_marker.txt','\t',header=None)
	if (df.shape[1]>7):
		df=df.iloc[:,0:7]
		df.to_csv('Female_marker.txt', header=None, sep='\t', index=False)
	df=pd.read_csv('Male_marker.txt','\t',header=None)
	if (df.shape[1]>7):
		df=df.iloc[:,0:7]
		df.to_csv('Male_marker.txt', header=None, sep='\t', index=False)
	with open(progenyID) as f:
		for line in f:
			tmp = line.strip().split('\t')
			cmd='python '+script+'/progenyGT.py '+'Female_marker.txt '+tmp[0]+'.vcf.gz '+tmp[0]+'.out'
			run_command(cmd)
			tmout=tmp[0]+'.out'
			os.remove(tmout)
			cmd='python '+script+'/progenyGT.py '+'Male_marker.txt '+tmp[0]+'.vcf.gz '+tmp[0]+'.out2'
			run_command(cmd)
			tmout=tmp[0]+'.out2'
			os.remove(tmout)
	
def genetype(x,list1,progenyID,out):
        df1 = pd.read_csv(x,header=None,sep='\t')
        df1=df1.values.tolist()
        df1=pd.DataFrame(df1)
        df2 = df1.fillna('--')
        namelist=['SSR_site','SSR_type','Number','Male','Male_SSR_type','Female','Female_SSR_type']
        list1=namelist+list1
        df2.columns=list1
        df2['progeny_All_types_of_SSR']='nan'
        list2=df2.values.tolist()
 #       print(list2)
        myc=''
        progany= pd.read_csv(progenyID,sep='\t')
        progany=progany.values.tolist()
        for i in range(len(list2)):
                mylist=[]
                for j in range(len(progany)+1):
                        if (list2[i][4]=='ab' and list2[i][6]=='ab' ):
                            myab=list2[i][3]
                            a=list2[i][7+j]
                            if (a=='--'):
                                genety='--'
                                mylist.append(genety)
                            elif (a.split("/")[0]!=a.split("/")[1] and (a.split("/")[0]==myab.split("/")[0] or a.split("/")[0]==myab.split("/")[1]) and (a.split("/")[1]==myab.split("/")[0] or a.split("/")[1]==myab.split("/")[1])):
                                genety='ab'
                                mylist.append(genety)
                            elif(a.split("/")[0]==a.split("/")[1] == myab.split("/")[0]):
                                genety='aa'
                                mylist.append(genety)
                            elif(a.split("/")[0]==a.split("/")[1] == myab.split("/")[1]):
                                genety='bb'
                                mylist.append(genety)
                            else:
                                genety='--'
                                mylist.append(genety)
                        elif(list2[i][4]=='ab' and list2[i][6]=='aa'):
                            a=list2[i][7+j]
                            myab=list2[i][3]
                            myab1=list2[i][5]
                            if (a=='--'):
                                genety='--'
                                mylist.append(genety)
                            elif (a.split("/")[0]!=a.split("/")[1] and (a.split("/")[0]==myab.split("/")[0] or a.split("/")[0]==myab.split("/")[1]) and (a.split("/")[1]==myab.split("/")[0] or a.split("/")[1]==myab.split("/")[1])):
                                genety='ab'
                                mylist.append(genety)
                            elif(a.split("/")[0]==a.split("/")[1] == myab1.split("/")[0]):
                                genety='aa'
                                mylist.append(genety)
                            else:
                                genety='--'
                                mylist.append(genety)
                        elif (list2[i][4]=='aa' and list2[i][6]=='bc' ):
                            a=list2[i][7+j]
                            myab=list2[i][3]#aa
                            myab1=list2[i][5]#bc
                            if (a=='--'):
                                genety='--'
                                mylist.append(genety)
                            elif (a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[0] ) or (a.split("/")[0]==myab1.split("/")[0] and a.split("/")[1]==myab.split("/")[0] ))):
                                genety='ab'
                                mylist.append(genety)
                            elif (a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[1] ) or (a.split("/")[0]==myab1.split("/")[1] and a.split("/")[1]==myab.split("/")[0] ))):
                                genety='ac'
                                mylist.append(genety)
                            else:
                                genety='--'
                                mylist.append(genety)
                        elif(list2[i][4]=='ab' and list2[i][6]=='cc'):
                            a=list2[i][7+j]
                            myab1=list2[i][3]#ab
                            myab=list2[i][5]#cc
                            if (a=='--'):
                                genety='--'
                                mylist.append(genety)
                            elif (a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[0] ) or (a.split("/")[0]==myab1.split("/")[0] and a.split("/")[1]==myab.split("/")[0] ))):
                                genety='ac'
                                mylist.append(genety)
                            elif (a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[1] ) or (a.split("/")[0]==myab1.split("/")[1] and a.split("/")[1]==myab.split("/")[0] ))):
                                genety='bc'
                                mylist.append(genety)
                            else:
                                genety='--'
                                mylist.append(genety)
                        elif(list2[i][4]=='ab' and list2[i][6]=='ac'):
                            a=list2[i][7+j]
                            mya=list2[i][3].split("/")[0]
                            myb=list2[i][3].split("/")[1]
                            if(list2[i][5].split("/")[0]==mya ):
                                myc=list2[i][5].split("/")[1]
                            elif(list2[i][5].split("/")[1]==mya ):
                                myc=list2[i][5].split("/")[0]
                            if (a=='--'):
                                genety='--'
                                mylist.append(genety)
                            ###aa ac ab bc
                            elif(a==mya+"/"+mya):
                                genety='aa'
                                mylist.append(genety)
                            elif(a==mya+"/"+myc or a==myc+"/"+mya):
                                genety='ac'
                                mylist.append(genety)
                            elif(a==mya+"/"+myb or a==myb+"/"+mya):
                                genety='ab'
                                mylist.append(genety)
                            elif(a==myc+"/"+myb or a==myb+"/"+myc):
                                genety='bc'
                                mylist.append(genety)
                            else:
                                genety='--'
                                mylist.append(genety)
                        elif(list2[i][4]=='ab' and list2[i][6]=='cd'):
                            a=list2[i][7+j]
                            mya=list2[i][3].split("/")[0]
                            myb=list2[i][3].split("/")[1]
                            myc=list2[i][5].split("/")[0]
                            myd=list2[i][5].split("/")[1]
                            if (a=='--'):
                                genety='--'
                                mylist.append(genety)
                            ###ac ad bc bd
                            elif(a==mya+"/"+myc or a==myc+"/"+mya):
                                genety='ac'
                                mylist.append(genety)
                            elif(a==mya+"/"+myd or a==myd+"/"+mya):
                                genety='ad'
                                mylist.append(genety)
                            elif(a==myb+"/"+myc or a==myc+"/"+myb):
                                genety='bc'
                                mylist.append(genety)
                            elif(a==myb+"/"+myd or a==myd+"/"+myb):
                                genety='bd'
                                mylist.append(genety)
                            else:
                                genety='--'
                                mylist.append(genety)

                myset=set(mylist)
                dict={}
                for item in myset:
                        dict.update({item:mylist.count(item)})
                list2[i][7+j+1]=str(dict).strip('{}')
        list1.append('progany_All_types_of_SSR')
        df1=pd.DataFrame(list2)
        df1.columns=list1
        df1=df1[df1['progany_All_types_of_SSR']!='']
        df1.to_csv(out,sep='\t',index=False)

def abxac(x,out2,y):
    with open (x,'r') as f:
        lines = f.readlines()
        first_line = lines[0]
        line_num=first_line.count('\t')
    progany_numble=line_num-7
    file=open(x,'r')
    open2=open(out2,'w')
    lines=file.readlines()
    for line in lines:
        tmp= line.strip().split('\t')
        if (tmp[4]=='ab' and tmp[6]=='ac' and ('aa' in tmp[line_num]) and ('ac' in tmp[line_num]) and ('ab' in tmp[line_num]) and ('bc' in tmp[line_num])):
             a=re.findall(r'aa\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype00=int(b)
             a=re.findall(r'ac\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype01=int(b)
             a=re.findall(r'ab\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype02=int(b)
             a=re.findall(r'bc\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype03=int(b)

             if (float(progany_numble-genetype00-genetype01-genetype02-genetype03)<float(float(MISPCT)*progany_numble)):
                 my_out=tmp[0]+"\t"+str(genetype00)+"\t"+str(genetype01)+"\t"+str(genetype02)+"\t"+str(genetype03)+"\n"
                 open2.write(my_out)
        elif (tmp[4]=='ab' and tmp[6]=='cd' and ('ac' in tmp[line_num]) and ('ad' in tmp[line_num]) and ('bc' in tmp[line_num]) and ('bd' in tmp[line_num])):
             a=re.findall(r'ac\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype00=int(b)
             a=re.findall(r'ad\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype01=int(b)
             a=re.findall(r'bc\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype02=int(b)
             a=re.findall(r'bd\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype03=int(b)

             if (float(progany_numble-genetype00-genetype01-genetype02-genetype03)<float(float(MISPCT)*progany_numble)):
                 my_out=tmp[0]+"\t"+str(genetype00)+"\t"+str(genetype01)+"\t"+str(genetype02)+"\t"+str(genetype03)+"\n"
                 open2.write(my_out)
    open2.close()
    file.close()
    myout2=out2+'.out'
    with open(out2,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 1):
            cmd='python ' +script+'/abxcd.py '+out2+' '+myout2
            run_command(cmd)
    if(os.path.isfile(myout2)==True):
        file = open(myout2, "r")
        z=y+'.joinmap'
        lines = file.readlines()
        flen=len(lines)
        if (flen < 2):
            print(myout2)
            return 0
        for line in lines:
            line=line.strip()
            tmp=line.split("\t")
            if (float(tmp[6])>float(pvalue)):
                get(tmp[0],y,z)
        file.close()
        if (y=='abxac.txt'):
            abxac_map(z)
        elif(y=='abxcd.txt'):
            abxcd_map(z)
def aaxbc(x,out2,y):
    with open (x,'r') as f:
        lines = f.readlines()
        first_line = lines[0]
        line_num=first_line.count('\t')
    progany_numble=line_num-7
    file=open(x,'r')
    open2=open(out2,'w')
    lines=file.readlines()
    for line in lines:
        tmp= line.strip().split('\t')
        if (tmp[4]=='aa' and tmp[6]=='bc' and ('ab' in tmp[line_num]) and ('ac' in tmp[line_num])):
             a=re.findall(r'ab\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype00=int(b)
             a=re.findall(r'ac\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype01=int(b)
             if (float(progany_numble-genetype00-genetype01)<float(float(MISPCT)*progany_numble)):
                 my_out=tmp[0]+"\t"+str(genetype00)+"\t"+str(genetype01)+"\n"
                 open2.write(my_out)
        elif(tmp[4]=='ab' and tmp[6]=='cc' and ('ac' in tmp[line_num]) and ('bc' in tmp[line_num])):
             a=re.findall(r'ac\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype00=int(b)
             a=re.findall(r'bc\': \d+',tmp[line_num])
             b="".join(a)
             b="".join(re.findall(r' \d+',b))
             genetype01=int(b)
             if (float(progany_numble-genetype00-genetype01)<float(float(MISPCT)*progany_numble)):
                 my_out=tmp[0]+"\t"+str(genetype00)+"\t"+str(genetype01)+"\n"
                 open2.write(my_out)

    open2.close()
    file.close()
    myout2=out2+'.out'
    with open(out2,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 1):
            cmd='python ' +script+'/aaxbc.py '+out2+' '+myout2
            run_command(cmd)
    if(os.path.isfile(myout2)==True):
        file = open(myout2, "r")
        z=y+'.joinmap'
        lines = file.readlines()
        flen=len(lines)
        if (flen < 2):
            return 0
        for line in lines:
            line=line.strip()
            tmp=line.split("\t")
            if (float(tmp[4])>float(pvalue)):
                get(tmp[0],y,z)
        file.close()
        if (y=='aaxbc.txt'):
            aaxbc_map(z)
        elif(y=='abxcc.txt'):
            abxcc_map(z)

def abxac_map(z):
    if(os.path.isfile(z)!=True):
        return 0
    with open(z,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 0):
            df1 = pd.read_csv(z,header=None,sep='\t')
            list1=df1.values.tolist()
            for i in range(len(list1)):
                 list1[i][1]='abxac'
                 list1[i][0]=str(list1[i][0])
                 list1[i][2]=''
            df1=pd.DataFrame(list1)
            df1=df1.sort_index(axis=1)
            df1 = df1.fillna('--')
            list2=df1.values.tolist()
            progany= pd.read_csv(progenyID,sep='\t')
            progany=progany.values.tolist()
            for i in range(len(list2)):
                for j in range(len(progany)+1):
                    a=list2[i][7+j]
                    mya=list2[i][3].split("/")[0]
                    myb=list2[i][3].split("/")[1]
                    myc=''
                    if(list2[i][5].split("/")[0]==mya ):
                        myc=list2[i][5].split("/")[1]
                    elif(list2[i][5].split("/")[1]==mya ):
                        myc=list2[i][5].split("/")[0]
                    if (a!='--' and (a==mya+"/"+mya )):
                        list2[i][7+j]='aa'
                    elif (a!='--' and (a==mya+"/"+myc or a==myc+"/"+mya) ):
                        list2[i][7+j]='ac'
                    elif (a!='--' and (a==myc+"/"+myb or a==myb+"/"+myc) ):
                        list2[i][7+j]='bc'
                    elif (a!='--' and (a==mya+"/"+myb or a==myb+"/"+mya) ):
                        list2[i][7+j]='ba'
                    else:
                        list2[i][7+j]='--'

                list2[i][3]='(ee,eg,fg,fe)'
            df1=pd.DataFrame(list2)
            list3 = []
            get_proganyID(list3,progenyID)
            namelist=['SSR_Locus','Segregation pattern','Phase','classification','4','5','6']
            list3=namelist+list3
            df1.columns=list3
            df1=df1.drop(['4','5','6'],axis=1)
            df1.to_csv('abxac.joinmap.txt',sep='\t',header=True,index=False)
            os.remove(z)
def abxcd_map(z):
    if(os.path.isfile(z)!=True):
        return 0

    with open(z,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 0):
            df1 = pd.read_csv(z,header=None,sep='\t')
            list1=df1.values.tolist()
            for i in range(len(list1)):
                 list1[i][1]='abxcd'
                 list1[i][0]=str(list1[i][0])
                 list1[i][2]=''
            df1=pd.DataFrame(list1)
            df1=df1.sort_index(axis=1)
            df1 = df1.fillna('--')
            list2=df1.values.tolist()
            progany= pd.read_csv(progenyID,sep='\t')
            progany=progany.values.tolist()
            for i in range(len(list2)):
                for j in range(len(progany)+1):
                    a=list2[i][7+j]
                    mya=list2[i][3].split("/")[0]
                    myb=list2[i][3].split("/")[1]
                    myc=list2[i][5].split("/")[0]
                    myd=list2[i][5].split("/")[1]
                    if (a!='--' and (a==mya+"/"+myc or a==myc+"/"+mya)):
                        list2[i][7+j]='ac'
                    elif (a!='--' and (a==mya+"/"+myd or a==myd+"/"+mya) ):
                        list2[i][7+j]='ad'
                    elif (a!='--' and (a==myb+"/"+myc or a==myc+"/"+myb) ):
                        list2[i][7+j]='bc'
                    elif (a!='--' and (a==myb+"/"+myd or a==myd+"/"+myb) ):
                        list2[i][7+j]='bd'
                    else:
                        list2[i][7+j]='--'

                list2[i][3]='(ac,ad,bc,bd)'
            df1=pd.DataFrame(list2)
            list3 = []
            get_proganyID(list3,progenyID)
            namelist=['SSR_Locus','Segregation pattern','Phase','classification','4','5','6']
            list3=namelist+list3
            df1.columns=list3
            df1=df1.drop(['4','5','6'],axis=1)
            df1.to_csv('abxcd.joinmap.txt',sep='\t',header=True,index=False)
            os.remove(z)
def aaxbc_map(z):
    if(os.path.isfile(z)!=True):
        return 0
    with open(z,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 0):
            df1 = pd.read_csv(z,header=None,sep='\t')
            list1=df1.values.tolist()
            for i in range(len(list1)):
                 list1[i][1]='aaxbc'
                 list1[i][0]=str(list1[i][0])
                 list1[i][2]=''
            df1=pd.DataFrame(list1)
            df1=df1.sort_index(axis=1)
            df1 = df1.fillna('--')
            list2=df1.values.tolist()
            progany= pd.read_csv(progenyID,sep='\t')
            progany=progany.values.tolist()
            for i in range(len(list2)):
                for j in range(len(progany)+1):
                    a=list2[i][7+j]
                    myab=list2[i][3]#aa
                    myab1=list2[i][5]#bc
                    if (a!='--' and a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[0] ) or (a.split("/")[0]==myab1.split("/")[0] and a.split("/")[1]==myab.split("/")[0] ))):
                        list2[i][7+j]='ab'
                    elif (a!='--' and a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[1] ) or (a.split("/")[0]==myab1.split("/")[1] and a.split("/")[1]==myab.split("/")[0]))):
                        list2[i][7+j]='ac'
                    else:
                        list2[i][7+j]='--'
                list2[i][3]='(nn,np)'
            df1=pd.DataFrame(list2)
            list3 = []
            get_proganyID(list3,progenyID)
            namelist=['SSR_Locus','Segregation pattern','Phase','classification','4','5','6']
            list3=namelist+list3
            df1.columns=list3
            df1=df1.drop(['4','5','6'],axis=1)
            df1.to_csv('aaxbc.joinmap.txt',sep='\t',header=True,index=False)
            os.remove(z)
def abxcc_map(z):
    if(os.path.isfile(z)!=True):
        return 0
    with open(z,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 0):
            df1 = pd.read_csv(z,header=None,sep='\t')
            list1=df1.values.tolist()
            for i in range(len(list1)):
                 list1[i][1]='abxcc'
                 list1[i][0]=str(list1[i][0])
                 list1[i][2]=''
            df1=pd.DataFrame(list1)
            df1=df1.sort_index(axis=1)
            df1 = df1.fillna('--')
            list2=df1.values.tolist()
            progany= pd.read_csv(progenyID,sep='\t')
            progany=progany.values.tolist()
            for i in range(len(list2)):
                for j in range(len(progany)+1):
                    a=list2[i][7+j]
                    myab1=list2[i][3]#ab
                    myab=list2[i][5]#cc
                    if (a!='--' and a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[0] ) or (a.split("/")[0]==myab1.split("/")[0] and a.split("/")[1]==myab.split("/")[0] ))):
                        list2[i][7+j]='bc'
                    elif(a!='--' and a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[1] ) or (a.split("/")[0]==myab1.split("/")[1] and a.split("/")[1]==myab.split("/")[0] ))):
                        list2[i][7+j]='ac'
                    else:
                        list2[i][7+j]='--'
                list2[i][3]='(lm,ll)'
            df1=pd.DataFrame(list2)
            list3 = []
            get_proganyID(list3,progenyID)
            namelist=['SSR_Locus','Segregation pattern','Phase','classification','4','5','6']
            list3=namelist+list3
            df1.columns=list3
            df1=df1.drop(['4','5','6'],axis=1)
            df1.to_csv('abxcc.joinmap.txt',sep='\t',header=True,index=False)
            os.remove(z)
def print_out(line_num,tmp,progany_numble,out2):
        if ('ab' in tmp[line_num] and 'aa' in tmp[line_num]):
                a=re.findall(r'aa\': \d+',tmp[line_num])
                b="".join(a)
                b="".join(re.findall(r' \d+',b))
                genetype00=int(b)
                a=re.findall(r'ab\': \d+',tmp[line_num])
                b="".join(a)
                b="".join(re.findall(r' \d+',b))
                genetype01=int(b)
                if (float(progany_numble-genetype00-genetype01)<float(float(MISPCT)*progany_numble)):
                        my_out=tmp[0]+"\t"+str(genetype00)+"\t"+str(genetype01)+"\n"
                        out2.write(my_out)

def print_out3(line_num,tmp,open3,progany_numble):
        a=re.findall(r'aa\': \d+',tmp[line_num])
        b="".join(a)
        b="".join(re.findall(r' \d+',b))
        genetype00=int(b)
        a=re.findall(r'ab\': \d+',tmp[line_num])
        b="".join(a)
        b="".join(re.findall(r' \d+',b))
        genetype01=int(b)
        a=re.findall(r'bb\': \d+',tmp[line_num])
        b="".join(a)
        b="".join(re.findall(r' \d+',b))
        genetype11=int(b)
        if (float(progany_numble-genetype00-genetype11-genetype01)<float(float(MISPCT)*progany_numble)):
            my_out=tmp[0]+"\t"+str(genetype00)+"\t"+str(genetype01)+"\t"+str(genetype11)+"\n"
            open3.write(my_out)


def pchisq(x,out2,out3):
	list1 = []
	get_proganyID(list1,progenyID)
	if (x=='Male_marker.out'):
		genetype('Male_marker.txt',list1,progenyID,'Male_marker.out')
	elif(x=='Female_marker.out'):
		genetype('Female_marker.txt',list1,progenyID,'Female_marker.out')
	myout2=out2+'.out'
	myout3=out3+'.out'
	with open (x,'r') as f:
		lines = f.readlines()
		first_line = lines[0]
		line_num=first_line.count('\t')
	progany_numble=line_num-7
	file=open(x,'r')
	open2=open(out2,'w')
	open3=open(out3,'w')
	lines=file.readlines()
	for line in lines:
		tmp= line.strip().split('\t')
		if (tmp[4]=='ab' and tmp[6]=='aa'):
			print_out(line_num,tmp,progany_numble,open2)               
		elif(tmp[4]=='ab' and tmp[6]=='ab'):
			if ('aa' in tmp[line_num] and 'ab' in tmp[line_num] and 'bb' in tmp[line_num] ):
				print_out3(line_num,tmp,open3,progany_numble)
	file.close()
	open2.close()
	open3.close()
	cmd='python ' +script+'/pchisq.py '+out2+' '+out3+' '+myout2+' '+myout3
	run_command(cmd)
def ALL_type(y):
	df1 = pd.read_table(y,header=None)
	list1=df1.values.tolist()
	df1=pd.DataFrame(list1)
	df1 = df1.fillna('--')
	list2=df1.values.tolist()
	progany= pd.read_table(progenyID)
	progany=progany.values.tolist()
	for i in range(len(list2)):
		ssr=''.join(re.findall(r'[ATCG]', list2[i][1]))
		if (change(list2[i][3].split('/')[0],ssr) and change(list2[i][3].split('/')[1],ssr) and change(list2[i][5].split('/')[0],ssr) and change(list2[i][5].split('/')[1],ssr)):
			list2[i][3]=change(list2[i][3].split('/')[0],ssr)+'/'+change(list2[i][3].split('/')[1],ssr)
			list2[i][5]=change(list2[i][5].split('/')[0],ssr)+'/'+change(list2[i][5].split('/')[1],ssr)
		for j in range(len(progany)+1):            
			a=list2[i][7+j]
			if (a!='--' and change(a.split('/')[0],ssr) and change(a.split('/')[1],ssr)): #10/15
				list2[i][7+j]=change(a.split('/')[0],ssr)+'/'+change(a.split('/')[1],ssr)   
	df1=pd.DataFrame(list2)
	list3 = []
	get_proganyID(list3,progenyID)
	if (y=='Male_marker.txt'):
		namelist=['SSR_Locus','SSR_type','SSR_Number','Male_alleles','Male_SSR_type','Female_alleles','Female_SSR_type']
	elif(y=='Female_marker.txt'):
		namelist=['SSR_Locus','SSR_type','SSR_Number','Female_alleles','Female_SSR_type','Male_alleles','Male_SSR_type']
	else:
		print("Male_marker.txt or Female_marker.txt is wrong")
	list3=namelist+list3
	df1.columns=list3
	if (y=='Male_marker.txt'):
		df1.to_csv('Male_marker.txt.1',sep='\t',header=False,index=False)
	elif(y=='Female_marker.txt'):
		df1.to_csv('Female_marker.txt.1',sep='\t',header=False,index=False)

def change(s,a):
    c=re.finditer(r'(%s)+'%(a),s)
    find_list=[]
    for i in c:
        find_list.append(i.group())
    if (len(find_list)==1):
        gt=a+'('+str(s.count(a))+')'
        return (gt)
    elif (len(find_list)==2):
        start=len(find_list[0])
        end=len(s)-len(find_list[1])
        gt=a+'('+str(find_list[0].count(a))+')'+'('+s[start:end]+')'+a+'('+str(find_list[1].count(a))+')'
        return (gt)
    elif (len(find_list)>2):
        return 0
def get(x,y,z):
    f = open(y,'r')
    lines = f.readlines()
    for line in lines:
        if x in  line:
            with open(z, 'a') as file:
                file.write(line)

def write2map(x,y,z,hybrid_pchis,hybrid_map,abxaa_Fslinkmap,abxab_Fslinkmap):
    file = open(x, "r")
    lines = file.readlines()
    for line in lines:
        line=line.strip()
        tmp=line.split("\t")
        if (float(tmp[4])>float(pvalue)):
            get(tmp[0],y,z)
    file.close()
    if(os.path.isfile(z)==True):
        df1 = pd.read_csv(z,header=None,sep='\t')
        list1=df1.values.tolist()
        progany= pd.read_csv(progenyID,sep='\t')
        progany=progany.values.tolist()
        for i in range(len(list1)):
            list1[i][1]='abxaa'
            list1[i][0]=str(list1[i][0])
            list1[i][2]=''
        df1=pd.DataFrame(list1)
        df1=df1.sort_index(axis=1)
        df1 = df1.fillna('--')
        list2=df1.values.tolist()
        for i in range(len(list2)):
            myab=list2[i][3]
            for j in range(len(progany)+1):
                a=list2[i][7+j]
                if (a!='--' and a.split("/")[0]!=a.split("/")[1] and (a.split("/")[0]==myab.split("/")[0] or a.split("/")[0]==myab.split("/")[1]) and (a.split("/")[1]==myab.split("/")[0] or a.split("/")[1]==myab.split("/")[1])):
                    list2[i][7+j]='ab'
                elif(a!='--' and a.split("/")[0]==a.split("/")[1] and (a.split("/")[0] == list2[i][5].split("/")[0])):
                    list2[i][7+j]='aa'
                else:
                    list2[i][7+j]='--'
            list2[i][3]='(aa,ab)'
        df1=pd.DataFrame(list2)
        list3 = []
        get_proganyID(list3,progenyID)
        namelist=['SSR_Locus','Segregation_pattern','Phase','classification','4','5','6']
        list3=namelist+list3
        df1.columns=list3
        df1=df1.drop(['4','5','6'],axis=1)
        df1.to_csv(z,sep='\t',header=True,index=False)
        Fslinkmap = pd.read_csv(z,header=None,sep='\t')
        list1=Fslinkmap.values.tolist()
        for i in range(1,len(list1)):
            list1[i][0]='*'+list1[i][0]
            list1[i][1]='abxaa'
        Fslinkmap=pd.DataFrame(list1)
        df1=Fslinkmap.drop([2,3],axis=1)
        df1.to_csv(abxaa_Fslinkmap,sep='\t',header=None,index=False)
        if (y=='Male_marker.txt'):
            cmd="sed 's/\t/  /g' Male_abxaa_Fslinkmap.txt >aaxab_Fslinkmap.txt"
            run_command(cmd)
            cmd="sed -i 's/abxaa/aaxab/g' aaxab_Fslinkmap.txt"
            run_command(cmd)
            cmd="sed -i 's/abxaa/aaxab/g' aaxab.joinmap.txt"
            run_command(cmd)
        elif(y=='Female_marker.txt'):
            cmd="sed  's/\t/  /g' Female_abxaa_Fslinkmap.txt >abxaa_Fslinkmap.txt"
            run_command(cmd)       
    if(os.path.isfile(hybrid_pchis)):
        # print("There were no data consistent with abxab segregation pattern")
#    else:
        file = open(hybrid_pchis, "r")
        lines = file.readlines()
        for line in lines:
            line=line.strip()
            tmp=line.split("\t")
            if (float(tmp[5])>float(pvalue)):
                get(tmp[0],y,hybrid_map)
        file.close()
        if(os.path.isfile(hybrid_map)):
#            print("There were no data consistent with abxab segregation pattern")
#        else:
            abxab(hybrid_map,abxab_Fslinkmap)
    
def abxab(hybrid_map,abxab_Fslinkmap):
    with open(hybrid_map,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 0):
            df1 = pd.read_csv(hybrid_map,header=None,sep='\t')
            list1=df1.values.tolist()
            for i in range(len(list1)):
                 list1[i][1]='abxab'
                 list1[i][0]=str(list1[i][0])
                 list1[i][2]=''
            df1=pd.DataFrame(list1)
            df1=df1.sort_index(axis=1)
            df1 = df1.fillna('--')
            list2=df1.values.tolist()
            progany= pd.read_csv(progenyID,sep='\t')
            progany=progany.values.tolist()
            for i in range(len(list2)):
                myab=list2[i][3]
                for j in range(len(progany)+1):
                    a=list2[i][7+j]
                    if (a!='--' and a.split("/")[0]!=a.split("/")[1] and (a.split("/")[0]==myab.split("/")[0] or a.split("/")[0]==myab.split("/")[1]) and (a.split("/")[1]==myab.split("/")[0] or a.split("/")[1]==myab.split("/")[1])):
                        list2[i][7+j]='ab'
                    elif(a!='--' and a.split("/")[0]==a.split("/")[1] and (a.split("/")[0] == myab.split("/")[0] or a.split("/")[0] == myab.split("/")[1])):
                        if(a.split("/")[0] == myab.split("/")[0]):
                            list2[i][7+j]='aa'
                        elif(a.split("/")[0] == myab.split("/")[1]):
                            list2[i][7+j]='bb'
                    else:
                        list2[i][7+j]='--'
                list2[i][3]='(aa,ab,bb)'
            df1=pd.DataFrame(list2)
            list3 = []
            get_proganyID(list3,progenyID)
            namelist=['SSR_Locus','Segregation pattern','Phase','classification','4','5','6']
            list3=namelist+list3
            df1.columns=list3
            df1=df1.drop(['4','5','6'],axis=1)
            df1.to_csv(hybrid_map,sep='\t',header=True,index=False)
            Fslinkmap = pd.read_csv(hybrid_map,header=None,sep='\t')
            list1=Fslinkmap.values.tolist()
            for i in range(1,len(list1)):
                list1[i][0]='*'+list1[i][0]
                list1[i][1]='abxab'
            Fslinkmap=pd.DataFrame(list1)
            df1=Fslinkmap.drop([2,3],axis=1)
            df1.to_csv(abxab_Fslinkmap,sep='\t',header=None,index=False)
            if(hybrid_map=='Female.abxab.joinmap.txt'):
                cmd="sed  's/\t/  /g' Female_abxab_Fslinkmap.txt >abxab_Fslinkmap.txt "
                run_command(cmd)
                cmd="less -S Female.abxab.joinmap.txt >abxab.joinmap.txt"
                run_command(cmd)
def get_segregation(x,y):
    f=open(x,'r')
    lines=f.readlines()
    for line in lines:
        tmp=line.strip().split('\t')
        if (tmp[4]=='ab' and tmp[6]=='aa'):
            with open('abxaa.txt', 'a+') as file:
                file.write(line)
        elif(tmp[4]=='ab' and tmp[6]=='ab'):
            with open('abxab.txt', 'a+') as file:
                file.write(line)
        elif(tmp[4]=='ab' and tmp[6]=='cd'):
            with open('abxcd.txt', 'a+') as file:
                file.write(line)
        elif(tmp[4]=='ab' and tmp[6]=='cc'):
            with open('abxcc.txt', 'a+') as file:
                file.write(line)
        elif(tmp[4]=='ab' and tmp[6]=='ac'):
            with open('abxac.txt', 'a+') as file:
                file.write(line)
    f.close()
    df1 = pd.read_csv(y,header=None,sep='\t')
    list1=df1.values.tolist()
    for i in range(len(list1)):
        tmp=list1[i][3]
        list1[i][3]=list1[i][5]
        list1[i][5]=tmp
        tmp=list1[i][4]
        list1[i][4]=list1[i][6]
        list1[i][6]=tmp
    df1=pd.DataFrame(list1)
    df1.to_csv('Male_marker.txt.2',sep='\t',header=False,index=False)
    f=open('Male_marker.txt.2','r')
    lines=f.readlines()
    for line in lines:
        tmp=line.strip().split('\t')
        if (tmp[4]=='aa' and tmp[6]=='ab'):
            with open('aaxab.txt', 'a+') as file:
                file.write(line)
        elif(tmp[4]=='cc' and tmp[6]=='ab'):
            with open('aaxbc.txt', 'a+') as file:
                file.write(line)

    f.close()
    if(os.path.isfile('aaxbc.txt')==True):
        df1 = pd.read_csv('aaxbc.txt',header=None,sep='\t')
        list1=df1.values.tolist()
        for i in range(len(list1)):
            list1[i][4]='aa'
            list1[i][6]='bc'
        df1=pd.DataFrame(list1)
        df1.to_csv('aaxbc.txt',sep='\t',header=False,index=False)
    list1 = []
    get_proganyID(list1,progenyID)
    if(os.path.isfile('abxaa.txt')==True):
        chang_head('abxaa.txt')
    if(os.path.isfile('aaxab.txt')==True):
        chang_head('aaxab.txt')
    if(os.path.isfile('abxab.txt')==True):
        chang_head('abxab.txt')
    if(os.path.isfile('abxcd.txt')==True):
        genetype('abxcd.txt',list1,progenyID,'abxcd.marker.out')
        abxac('abxcd.marker.out','abxcd.marker.out.pchis','abxcd.txt')
        chang_head('abxcd.txt')
    if(os.path.isfile('aaxbc.txt')==True):
        genetype('aaxbc.txt',list1,progenyID,'aaxbc.marker.out')
        aaxbc('aaxbc.marker.out','aaxbc.marker.out.pchis','aaxbc.txt')
        chang_head('aaxbc.txt')
    if(os.path.isfile('abxcc.txt')==True):
        genetype('abxcc.txt',list1,progenyID,'abxcc.marker.out')
        aaxbc('abxcc.marker.out','abxcc.marker.out.pchis','abxcc.txt')
        chang_head('abxcc.txt')
    if(os.path.isfile('abxac.txt')==True):
        genetype('abxac.txt',list1,progenyID,'abxac.marker.out')
        abxac('abxac.marker.out','abxac.marker.out.pchis','abxac.txt')
        chang_head('abxac.txt')

def chang_head(x):
    list1 = []
    get_proganyID(list1,progenyID)
    df1 = pd.read_csv(x,header=None,sep='\t')
    df1=df1.values.tolist()
    df1=pd.DataFrame(df1)
    namelist=['SSR_site','Ref_SSR_type','Number','Female','Female_segregation_type','Male','Male_segregation_type',]
    list1=namelist+list1
    df1.columns=list1
    df1.to_csv(x,sep='\t',header=True,index=False)
def getfastiq(x):
    global male1,male2,female1,female2
    Myproid=open('Myprogeny.id','w')
    Myproid.close()
    with open(x,'r') as f:
        i=0
        for line in f:
            i=i+1
            if (line.startswith("[data files]")):
                flag=i+1
                break
    with open(x,'r') as f:
        f =f.readlines()[flag+2:]
        for line in f:
            line=line.strip()
            myout=line.split(":")[0]+'\t'+line.split()[0].split(":")[1]+'\t'+line.split()[1]
            with open('Myprogeny.id',"a+") as f2:
                f2.write(myout+"\n")
    with open(x,'r') as f:
        f =f.readlines()[flag:flag+2]
        for line in f:
            line=line.strip()
            if (line.split(":")[0]=='MALE'):
                male1=parentfold+'/'+line.split()[0].split(":")[1]
                male2=parentfold+'/'+line.split()[1]
            if (line.split(":")[0]=='FEMALE'):
                female1=parentfold+'/'+line.split()[0].split(":")[1]
                female2=parentfold+'/'+line.split()[1]
    if(os.path.isfile('Myprogeny.id')):
        return ('./Myprogeny.id')

def changaaxab(dict,df1,z):
    df1.insert(1, 'E', 1)
    list1 = []
    get_proganyID(list1,progenyID)
    num=len(list1)
    flag=5+num
    list1 = []
    for i in range(10+num):
        list1.append(str(i))
    df1.columns=list1
    df1.drop([0],axis = 0, inplace = True)
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][1]="aa"
        list2[i][2]="ab"
        list2[i][3]=dict[list2[i][0]].split("/")[0]
        list2[i][4]=dict[list2[i][0]].split("/")[1]
        list2[i][flag]=int(list2[i][flag])
        list2[i][flag+1]=int(list2[i][flag+1])
        list2[i][flag+2]=num-list2[i][flag]-list2[i][flag+1]
    df1=pd.DataFrame(list2)
    order=[0,3,4,flag,flag+1,flag+2,flag+4,flag+3,1,2]
    list1 = []
    for i in range(5,flag):
         list1.append(i)
    order=order+list1
    df1=df1[order]
    list3 = []
    get_proganyID(list3,progenyID)
    namelist=['Position','Allele_a' ,'Allele_b','aa','ab','--','x2','p-value','female','male']
    list3=namelist+list3
    df1.columns=list3
    df1.to_csv(z,sep='\t',header=True,index=False)
def changabxaa(dict,df1,z):
    df1.insert(1, 'E', 1)
    list1 = []
    get_proganyID(list1,progenyID)
    num=len(list1)
    flag=5+num
    list1 = []
    for i in range(10+num):
        list1.append(str(i))
    df1.columns=list1
    df1.drop([0],axis = 0, inplace = True)
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][1]="ab"
        list2[i][2]="aa"
        list2[i][3]=dict[list2[i][0]].split("/")[0]
        list2[i][4]=dict[list2[i][0]].split("/")[1]
        list2[i][flag]=int(list2[i][flag])
        list2[i][flag+1]=int(list2[i][flag+1])
        list2[i][flag+2]=num-list2[i][flag]-list2[i][flag+1]
    df1=pd.DataFrame(list2)
    order=[0,3,4,flag,flag+1,flag+2,flag+4,flag+3,1,2]
    list1 = []
    for i in range(5,flag):
         list1.append(i)
    order=order+list1
    df1=df1[order]
    list3 = []
    get_proganyID(list3,progenyID)
    namelist=['Position','Allele_a' ,'Allele_b','ab','aa','--','x2','p-value','female','male']
    list3=namelist+list3
    df1.columns=list3
    df1.to_csv(z,sep='\t',header=True,index=False)
def changabxab(dict,df1,z):
    df1.insert(1, 'E', 1)
    list1 = []
    get_proganyID(list1,progenyID)
    num=len(list1)
    flag=5+num
    list1 = []
    for i in range(11+num):
        list1.append(str(i))
    df1.columns=list1
    df1.drop([0],axis = 0, inplace = True)
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][1]="ab"
        list2[i][2]="ab"
        list2[i][3]=dict[list2[i][0]].split("/")[0]
        list2[i][4]=dict[list2[i][0]].split("/")[1]
        list2[i][flag]=int(list2[i][flag])
        list2[i][flag+1]=int(list2[i][flag+1])
        list2[i][flag+2]=int(list2[i][flag+2])
        list2[i][flag+3]=num-list2[i][flag]-list2[i][flag+1]-list2[i][flag+2]
    df1=pd.DataFrame(list2)
    order=[0,3,4,flag,flag+1,flag+2,flag+3,flag+5,flag+4,1,2]
    list1 = []
    for i in range(5,flag):
         list1.append(i)
    order=order+list1
    df1=df1[order]
    list3 = []
    get_proganyID(list3,progenyID)
    namelist=['Position','Allele_a' ,'Allele_b','aa','ab','bb','--','x2','p-value','female','male']
    list3=namelist+list3
    df1.columns=list3
    df1.to_csv(z,sep='\t',header=True,index=False)
def changaaxbc(dict,df1,z):
    df1.insert(1, 'E', 1)
    df1.insert(2, 'f', 2)
    list1 = []
    get_proganyID(list1,progenyID)
    num=len(list1)
    flag=6+num
    list1 = []
    for i in range(11+num):
        list1.append(str(i))
    df1.columns=list1
    df1.drop([0],axis = 0, inplace = True)
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][1]="aa"
        list2[i][2]="bc"
        list2[i][3]=dict[list2[i][0]].split("/")[0]
        list2[i][4]=dict[list2[i][0]].split("/")[1]
        list2[i][5]=dict[list2[i][0]].split("/")[2]
        list2[i][flag]=int(list2[i][flag])
        list2[i][flag+1]=int(list2[i][flag+1])
        list2[i][flag+2]=num-list2[i][flag]-list2[i][flag+1]
    df1=pd.DataFrame(list2)
    order=[0,3,4,5,flag,flag+1,flag+2,flag+4,flag+3,1,2]
    list1 = []
    for i in range(6,flag):
         list1.append(i)
    order=order+list1
    df1=df1[order]
    list3 = []
    get_proganyID(list3,progenyID)
    namelist=['Position','Allele_a' ,'Allele_b','Allele_c','aa','bc','--','x2','p-value','female','male']
    list3=namelist+list3
    df1.columns=list3
    df1.to_csv(z,sep='\t',header=True,index=False)
def changabxac(dict,df1,z):
    df1.insert(1, 'E', 1)
    df1.insert(2, 'f', 2)
    list1 = []
    get_proganyID(list1,progenyID)
    num=len(list1)
    flag=6+num
    list1 = []
    for i in range(13+num):
        list1.append(str(i))
    df1.columns=list1
    df1.drop([0],axis = 0, inplace = True)
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][1]="ab"
        list2[i][2]="ac"
        list2[i][3]=dict[list2[i][0]].split("/")[0]
        list2[i][4]=dict[list2[i][0]].split("/")[1]
        list2[i][5]=dict[list2[i][0]].split("/")[2]
        list2[i][flag]=int(list2[i][flag])
        list2[i][flag+1]=int(list2[i][flag+1])
        list2[i][flag+2]=int(list2[i][flag+2])
        list2[i][flag+3]=int(list2[i][flag+3])
        list2[i][flag+4]=num-list2[i][flag]-list2[i][flag+1]-list2[i][flag+2]-list2[i][flag+3]
    df1=pd.DataFrame(list2)
    order=[0,3,4,5,flag,flag+1,flag+2,flag+3,flag+4,flag+6,flag+5,1,2]
    list1 = []
    for i in range(6,flag):
         list1.append(i)
    order=order+list1
#    print(order)
    df1=df1[order]
    list3 = []
    get_proganyID(list3,progenyID)
    namelist=['Position','Allele_a' ,'Allele_b','Allele_c','aa','ac','ab','bc','--','x2','p-value','female','male']
    list3=namelist+list3
    df1.columns=list3
    df1.to_csv(z,sep='\t',header=True,index=False)
def changabxcc(dict,df1,z):
    df1.insert(1, 'E', 1)
    df1.insert(2, 'f', 2)
    list1 = []
    get_proganyID(list1,progenyID)
    num=len(list1)
    flag=6+num
    list1 = []
    for i in range(11+num):
        list1.append(str(i))
    df1.columns=list1
    df1.drop([0],axis = 0, inplace = True)
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][1]="ab"
        list2[i][2]="cc"
        list2[i][3]=dict[list2[i][0]].split("/")[0]
        list2[i][4]=dict[list2[i][0]].split("/")[1]
        list2[i][5]=dict[list2[i][0]].split("/")[2]
        list2[i][flag]=int(list2[i][flag])
        list2[i][flag+1]=int(list2[i][flag+1])
        list2[i][flag+2]=num-list2[i][flag]-list2[i][flag+1]
    df1=pd.DataFrame(list2)
    order=[0,3,4,5,flag,flag+1,flag+2,flag+4,flag+3,1,2]
    list1 = []
    for i in range(6,flag):
         list1.append(i)
    order=order+list1
    df1=df1[order]
    list3 = []
    get_proganyID(list3,progenyID)
    namelist=['Position','Allele_a' ,'Allele_b','Allele_c','ac','bc','--','x2','p-value','female','male']
    list3=namelist+list3
    df1.columns=list3
    df1.to_csv(z,sep='\t',header=True,index=False)
def changabxcd(dict,df1,z):
    df1.insert(1, 'E', 1)
    df1.insert(2, 'f', 2)
    df1.insert(3, 'g', 3)
    list1 = []
    get_proganyID(list1,progenyID)
    num=len(list1)
    flag=7+num
    list1 = []
    for i in range(14+num):
        list1.append(str(i))
    df1.columns=list1
    df1.drop([0],axis = 0, inplace = True)
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][1]="ab"
        list2[i][2]="cd"
        list2[i][3]=dict[list2[i][0]].split("/")[0]
        list2[i][4]=dict[list2[i][0]].split("/")[1]
        list2[i][5]=dict[list2[i][0]].split("/")[2]
        list2[i][6]=dict[list2[i][0]].split("/")[3]
        list2[i][flag]=int(list2[i][flag])
        list2[i][flag+1]=int(list2[i][flag+1])
        list2[i][flag+2]=int(list2[i][flag+2])
        list2[i][flag+3]=int(list2[i][flag+3])
        list2[i][flag+4]=num-list2[i][flag]-list2[i][flag+1]-list2[i][flag+2]-list2[i][flag+3]
    df1=pd.DataFrame(list2)
    order=[0,3,4,5,6,flag,flag+1,flag+2,flag+3,flag+4,flag+6,flag+5,1,2]
    list1 = []
    for i in range(7,flag):
         list1.append(i)
    order=order+list1
#    print(order)
    df1=df1[order]
    list3 = []
    get_proganyID(list3,progenyID)
    namelist=['Position','Allele_a' ,'Allele_b','Allele_c','Allele_d','ac','ad','bc','bd','--','x2','p-value','female','male']
    list3=namelist+list3
    df1.columns=list3
    df1.to_csv(z,sep='\t',header=True,index=False)
def Merge(x,y,z):
    if(os.path.isfile(x) and os.path.isfile(y) and os.path.isfile(z)):
        dict={}
        with open (z) as f:
            lines=f.readlines()[1:]
            for line in lines:
                tmp=line.strip().split('\t')
                if(z=='abxaa.txt' or z=='abxab.txt'):
                    dict[tmp[0]]=tmp[3]
                elif(z=='aaxab.txt'):
                    dict[tmp[0]]=tmp[5]
                elif(z=='aaxbc.txt'):
                    myout=tmp[3].split('/')[0]+'/'+tmp[5]
                    dict[tmp[0]]=myout
                elif(z=='abxac.txt'):
                    if(tmp[5].split('/')[0]!=tmp[3].split('/')[0]!=tmp[3].split('/')[1] ):
                        dict[tmp[0]]=tmp[3]+'/'+tmp[5].split('/')[0]
                    elif(tmp[5].split('/')[1]!=tmp[3].split('/')[0]!=tmp[3].split('/')[1] ):
                        dict[tmp[0]]=tmp[3]+'/'+tmp[5].split('/')[1]
                elif(z=='abxcc.txt'):
                    myout=tmp[3]+'/'+tmp[5].split('/')[0]
                    dict[tmp[0]]=myout
                elif(z=='abxcd.txt'):
                    myout=tmp[3]+'/'+tmp[5]
                    dict[tmp[0]]=myout
        df1 = pd.read_csv(x,sep='\t',header=None)
        df2 = pd.read_csv(y,sep='\t',header=None)
        data = pd.merge(df1,df2,on=0,how="left")
        if(z=='abxaa.txt'):
            changabxaa(dict,data,z)
        elif(z=='aaxab.txt'):
            changaaxab(dict,data,z)
        elif(z=='aaxbc.txt'):
            changaaxbc(dict,data,z)
        elif(z=='abxab.txt'):
            changabxab(dict,data,z)
        elif(z=='abxac.txt'):
            changabxac(dict,data,z)
        elif(z=='abxcc.txt'):
            changabxcc(dict,data,z)
        elif(z=='abxcd.txt'):
            changabxcd(dict,data,z)
    else:
        if(os.path.isfile(z)):
            os.remove(z)
def populationtype(p):
    Merge('abxaa.joinmap.txt','Female_pure_pchis.out','abxaa.txt')
    Merge('aaxab.joinmap.txt', 'Male_pure_pchis.out' ,'aaxab.txt')
    Merge('aaxbc.joinmap.txt', 'aaxbc.marker.out.pchis.out' ,'aaxbc.txt')
    Merge('abxab.joinmap.txt','Female_hybrid_pchis.out','abxab.txt')
    Merge('abxac.joinmap.txt','abxac.marker.out.pchis.out','abxac.txt')
    Merge('abxcd.joinmap.txt','abxcd.marker.out.pchis.out','abxcd.txt')
    Merge('abxcc.joinmap.txt','abxcc.marker.out.pchis.out','abxcc.txt')
    os.mkdir('./WorkingDirectory')
    if (p == 'F2'):
        files =os.listdir()
        count=0
        for file in files:
            if (file.startswith('abxab')):
                cmd='mv a* ./WorkingDirectory'    
                run_command(cmd)
                cmd='mv ./WorkingDirectory/abxab* ./'
                run_command(cmd)
                break
            else:
                conut=count+1
                if (conut==len(files) or ('abxab.txt' not in files) ):
                    print("No SSR genotypes are found that fit Mendelian segregation ratios for abxab.")
                    break
    elif(p.split(':')[0] == 'BC'):
        if(p=='BC:female' or p=='BC:male'):
            if (p.split(':')[1] == 'male'):
                files =os.listdir()
                conut=0 
                for file in files:
                    if (file.startswith('abxaa')):
                        cmd='mv a* ./WorkingDirectory'
                        run_command(cmd)
                        cmd='mv ./WorkingDirectory/abxaa* ./'
                        run_command(cmd)
                        break
                    else:
                        conut=conut+1
                        if(conut==len(files) or ('abxaa.txt' not in files)):
                            print("No SSR genotypes are found that fit Mendelian segregation ratios for abxaa.")
                            break
            elif(p.split(':')[1] == 'female'):
                files =os.listdir()
                count=0
                for file in files:
                    if (file.startswith('aaxab')):
                        cmd='mv a* ./WorkingDirectory'
                        run_command(cmd)
                        cmd='mv ./WorkingDirectory/aaxab* ./'
                        run_command(cmd)
                        break
                    else:
                        count=count+1
                        if(count==len(files) or ('aaxab.txt' not in files)):
                            print("No SSR genotypes are found that fit Mendelian segregation ratios for aaxab.")
                            break
        else:
            sys.exit(' Error!Please input BC:female or BC:male')
    cmd='mv  *marker*  *map*  ./WorkingDirectory'
    run_command(cmd)
    cmd='mv F* M*  ./WorkingDirectory'
    run_command(cmd)
    if (os.path.exists('./WorkingDirectory/Female_marker.txt') and os.path.exists('./WorkingDirectory/Male_marker.txt')):
        shutil.move('./WorkingDirectory/Female_marker.txt','./')
        shutil.move('./WorkingDirectory/Male_marker.txt','./')
            
    files =os.listdir()
    for file in files:
        if (file.startswith('female') or file.startswith('male') ):
            cmd='mv female* male* ./WorkingDirectory'
            run_command(cmd)
            break
    files =os.listdir('./WorkingDirectory')
    file="".join(files)
    if (re.findall("fq",file)):
        cmd='mv ./WorkingDirectory/*fq  ./'
        run_command(cmd)
    cmd='rm -rf ./WorkingDirectory'
    run_command(cmd)
    if(p=='CP'):
        files=os.listdir()
        if('abxaa.txt' not in files):
            print("No SSR genotypes are found that fit Mendelian segregation ratios for abxaa.")
        if('aaxab.txt' not in files):
            print("No SSR genotypes are found that fit Mendelian segregation ratios for aaxab.")
        if('abxab.txt' not in files):
            print("No SSR genotypes are found that fit Mendelian segregation ratios for abxab.")


def filter(population):
        if (population == 'CP' or population == 'F2'or population == 'BC:female' or  population == 'BC:male'):
                if(os.path.isfile('SSR.txt')==False):
                        sys.exit(' Error!Please do not execute the --nosearch command,this step cannot be skipped now.')
                if(os.path.isfile('Female_marker.txt')==False or os.path.isfile('Male_marker.txt')==False):
                        sys.exit(' Error!Please do not execute the --nocatalog command,this step cannot be skipped now.')
                flag1=0
                with open(progenyID,'r') as f:
                        for line in f:
                                flag1=flag1+1
                                if(os.path.isfile(line.split('\t')[0]+'.sort.bam')==False ):
                                       sys.exit(' Error!Please do not execute the --nomap command,this step cannot be skipped now.')
                df1 = pd.read_table('Male_marker.txt',header=None)
                df1=df1.values.tolist()
                df1=pd.DataFrame(df1)
                if (df1.shape[1]<8):
                        sys.exit(' Error!Please do not execute the --nocall command,this step cannot be skipped now.')
                files =os.listdir()
                for file in files:
                        if (file.startswith('a') ):
                                cmd='rm a* '
                                run_command(cmd)
                                break
                pchisq('Male_marker.out','Male_pure_pchis','Male_hybrid_pchis')
                pchisq('Female_marker.out','Female_pure_pchis','Female_hybrid_pchis')
                ALL_type('Male_marker.txt')
                ALL_type('Female_marker.txt')
                get_segregation('Female_marker.txt.1','Male_marker.txt.1')
                write2map('Male_pure_pchis.out','Male_marker.txt','aaxab.joinmap.txt','Male_hybrid_pchis.out','Male.abxab.joinmap.txt','Male_abxaa_Fslinkmap.txt','Male_abxab_Fslinkmap.txt')
                write2map('Female_pure_pchis.out','Female_marker.txt','abxaa.joinmap.txt','Female_hybrid_pchis.out','Female.abxab.joinmap.txt','Female_abxaa_Fslinkmap.txt','Female_abxab_Fslinkmap.txt')
                populationtype(population)


def main(args):
	global progenyID,motif
	progenyID=getfastiq("parameters.ini")
	motif = args.motif
	population=args.population
	if (args.nosearch):
		get_SSR(GENOME)
	if(args.nocatalog):
		if(os.path.isfile('SSR.txt')==False):
			sys.exit(' Error!Please do not execute the --nosearch command,this step cannot be skipped now.')
		find_parent(GENOME)
		print("This step is calling SSR genotypes for the parents, please wait!")
		callparent2()
	if(args.nomap):
		Mapping_Pool()
	if(args.nocall):
		if(os.path.isfile('Female_marker.txt')==False or os.path.isfile('Male_marker.txt')==False):
			sys.exit(' Error!Please do not execute the --nocatalog command,this step cannot be skipped now.')
		with open(progenyID,'r') as f:
			for line in f:
				if(os.path.isfile(line.split('\t')[0]+'.sort.bam')==False ):
					sys.exit(' Error!Please do not execute the --nomap command,this step cannot be skipped now.')
		print("This step is calling SSR genotypes for the progeny, please wait!")
		SSRGM(GENOME)
	if(args.nofilter):
		if (population == 'CP' or population == 'F2'or population == 'BC:female' or  population == 'BC:male'):
			filter(population)
		else:
			sys.exit(' Error!Please input population type [CP,F2,BC:female,BC:male]')
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		prog="python SSRGT.py",
		formatter_class=argparse.RawTextHelpFormatter,
		description='''
-------------------------------------------------------------------------------------------------------
SSRGT provides a usable and accurate SSR genotyping platform for distant hybridization populations. 
It has the advantage that a large number of SSR loci segregated by Mendelian segregation ratio in a
population can be discovered with a simple command. 
It will facilitate the construction of SSR genetic linkage maps,
locating quantitative trait loci, marker-assisted selection breeding and various genetic studies.			

Contact: Chunfa Tong  <tongchf@njfu.edu.cn>
Version: 1.0
Usage: python SSRGT.py  [Options]

-------------------------------------------------------------------------------------------------------
''')
	parser.add_argument(
		'-mo', '--motif',
		default="1=11,2=5,3=4,4=4,5=4,6=4",
		help="set the threshold of motif, default: 1=11,2=5,3=4,4=4,5=4,6=4 \nleft  of equal : length of motif \nright of equal : the minimum number of repeat")
	parser.add_argument(
		'-p', '--population',
		default='CP',
		help="set the population type [CP,F2,BC:female,BC:male], default : CP \nIf your population type  is CP or F2, you should choose 'CP' or 'F2', if it's BC population type and the maternal parent is recurrent parent, choose 'BC:female', if it's BC population type and the paternal parent is recurrent parent, choose 'BC:male'")
	parser.add_argument("--nosearch", help = "skip the step for searching SSR in reference sequence.", default='True',action='store_false')
	parser.add_argument("--nocatalog", help = "skip the step for generating parental SSR catalogs.", default='True',action='store_false')
	parser.add_argument("--nomap", help = "skip the step for mapping the progeny reads to the reference sequence.", default='True',action='store_false')
	parser.add_argument("--nocall", help = "skip the step for calling  SSR genotypes for each individual.", default='True',action='store_false')
	parser.add_argument("--nofilter", help = "skip the step for filtering  SSR genotype data.", default='True',action='store_false')
	args = parser.parse_args()
	filename='./parameters.ini'
	if(os.path.isfile(filename)==False):
        	print("\n This program needs a parameter file,namely \"parameters.ini\".\n")
        	sys.exit(' Error!Please check the parameters.ini')
	if(os.path.isfile(GENOME)==False):
        	sys.exit(' Error!Please check  the line of REFERENCE_FILE in the parameters.ini')
	if(os.path.isfile(bwa)==False):
        	sys.exit(' Error!Please check  the line of BWA_FOLD in the parameters.ini')
	if(os.path.isfile(samtools)==False):
        	sys.exit(' Error!Please check  the line of SAMTOOLS_FOLD in the parameters.ini')
	if(os.path.exists(parentfold)==False ):
        	sys.exit(' Error!Please check  the line of RADDATA_FOLD in the parameters.ini')
	if(os.path.exists(progenyfold)):
		if(os.path.exists(script)):
			main(args)
			print("All tasks done")
		else:
			print(' Error!Please check  the line of SSRGT_FOLD in the parameters.ini. Please check if the script file directory exists.')
	else:
		print(' Error!Please check  the line of RADDATA_FOLD in the parameters.ini')
