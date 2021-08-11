#!/usr/bin/evn python3
import sys,os
import subprocess
import argparse
import shutil
import re
import pandas as pd
from multiprocessing import Pool
import time
import os.path
import xlwt
from configparser import ConfigParser
cfg = ConfigParser()
cfg.read("parameters.ini")
bwa=cfg.get("folders","BWA")
GangSTR=cfg.get("folders","GangSTR")
samtools=cfg.get("folders","SAMTOOLS")
thread=cfg.getint("parameter","THREADS")
Stutter=cfg.get("parameter","Stutter")
MAPQ=cfg.getint("parameter","MAPQ")
pvalue=cfg.get("parameter","PVALUE")
MISPCT=cfg.get("parameter","MISS_GENOTYPES")
GENOME=cfg.get("data files","REFERENCE_GENOME")
progenyfold=cfg.get("folders","PROGENY")
script=cfg.get("folders","SCRIPT")
SSRMMD=script+'/SSRMMD.pl'
PICARD=script+'/picard-tools-1.119'
DP=cfg.get("parameter","Depth_Of_Coverage")
Alleles_Quality=cfg.get("parameter","Alleles quality score")
DP=int(DP)
Alleles_Quality=float(Alleles_Quality)
male1=cfg.get("data files","MALEPARENT_1")
male2=cfg.get("data files","MALEPARENT_2")
female1=cfg.get("data files","FEMALEPARENT_1")
female2=cfg.get("data files","FEMALEPARENT_2")
###########################################################################################################

def run_command(cmd):
	print(cmd)
	return_code = subprocess.call(cmd, shell=True)

def get_SSR(reader1):
	cmd1 = 'perl '  +SSRMMD+' '+' -f1 '+reader1+' -l 100 -ss 1  -mo '+ motif +' -t '+str(thread)+' -o ./'
	run_command(cmd1)
	cmd2 = samtools + ' faidx ' + reader1  
	run_command(cmd2)
	SSR=reader1+'.SSRs'
	file = open(SSR, "r")
	lines = file.readlines()
	for line in lines:
		line=line.strip()
		if line.startswith("number"):
			continue
		else:
			tmp=line.split("\t")
			geneID=tmp[1]+'\t'+tmp[6]+'\t'+tmp[7]+'\t'+tmp[3]+'\t'+tmp[2]
			with open('myregions',"a+") as f2:
				f2.write(geneID+"\n")
	file.close()
def find_parent(reader1):
	cmd = bwa + ' index '	+reader1
	run_command(cmd)
	cmd='java -Xmx4g -jar -XX:ParallelGCThreads='+ str(thread)+' '+PICARD+'/CreateSequenceDictionary.jar R='+reader1+' O='+reader1+'.dict'
	run_command(cmd)
	cmd = bwa +  ' mem  -M -t '+ str(thread)+'  -R ' +"'@RG\\tID:EV\\tPL:ILLUMINA\\tLB:EV\\tSM:EV'"+' '+reader1+' '+male1+' '+male2+' > male.sam'
	run_command(cmd)
	cmd=samtools+" view -@ "+str(thread)+" -bS male.sam -q "+str(MAPQ)+" -o  male.bam"
	run_command(cmd)
	cmd='java -jar -XX:ParallelGCThreads='+ str(thread)+' '+PICARD+'/SortSam.jar INPUT=male.bam OUTPUT=male.sort.bam SORT_ORDER=coordinate'
	run_command(cmd)
	cmd='java -jar -XX:ParallelGCThreads='+str(thread)+' '+PICARD+'/MarkDuplicates.jar INPUT=male.sort.bam OUTPUT=male.sort.dedup.bam METRICS_FILE=male_dedup.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true REMOVE_DUPLICATES=true'
	run_command(cmd)
	cmd = bwa +  ' mem  -M -t '+ str(thread)+'  -R ' +"'@RG\\tID:EV\\tPL:ILLUMINA\\tLB:EV\\tSM:EV'"+' '+reader1+' '+female1+' '+female2+' > female.sam'
	run_command(cmd)
	cmd=samtools+" view -@ "+str(thread)+" -bS female.sam -q "+str(MAPQ)+" -o  female.bam"
	run_command(cmd)
	cmd='java -jar -XX:ParallelGCThreads='+str(thread)+' '+PICARD+'/SortSam.jar INPUT=female.bam OUTPUT=female.sort.bam SORT_ORDER=coordinate'
	run_command(cmd)
	cmd='java -jar -XX:ParallelGCThreads='+str(thread)+' '+PICARD+'/MarkDuplicates.jar INPUT=female.sort.bam OUTPUT=female.sort.dedup.bam METRICS_FILE=female_dedup.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true REMOVE_DUPLICATES=true'
	run_command(cmd)

	
def get_proganyID(list1,proganyID):
        with open(proganyID) as f:
                for line in f:
                        a=line.strip().split("\t")
                        list1.append(a[0])

def SSRGM(reader1):
			
	cmd='python   '+script+'/parent_merge.py '+'-parent1 male.vcf01.txt.out -parent2 femaleallSSR_type.txt -o female.out'
	run_command(cmd)
	cmd='python   '+script+'/parent_merge.py '+'-parent1 female.vcf01.txt.out -parent2 maleallSSR_type.txt -o  male.out'
	run_command(cmd)
	cmd='less -S female.out |sort|uniq > '+'female.sort.out'
	run_command(cmd)
	cmd='less -S male.out |sort|uniq > '+'male.sort.out'
	run_command(cmd)
	df1 = pd.read_csv('female.vcf01.txt.out',sep='\t',header=None)
	df2 = pd.read_csv('male.sort.out',sep='\t',header=None)
	data = pd.merge(df1,df2,on=0,how="left")
	df=data.dropna(axis=0,how='any')
	df.to_csv('Female_marker.txt',sep='\t',header=False,index=False)
	df1 = pd.read_csv('male.vcf01.txt.out',sep='\t',header=None)
	df2 = pd.read_csv('female.sort.out',sep='\t',header=None)
	data = pd.merge(df1,df2,on=0,how="left")
	df=data.dropna(axis=0,how='any')
	df.to_csv('Male_marker.txt',sep='\t',header=False,index=False)
	
	os.remove('female.sam')
	os.remove('male.sam')
	os.remove('female.bam')
	os.remove('male.bam')	
		
	with open(progenyID) as f:
		for line in f:
			tmp = line.strip().split('\t')
			cmd = bwa +  ' mem  -M -t '+ str(thread)+'  -R ' +"'@RG\\tID:EV\\tPL:ILLUMINA\\tLB:EV\\tSM:EV'"+' '+reader1+' '+progenyfold+'/'+tmp[1]+' '+progenyfold+'/'+tmp[2]+' > '+tmp[0]+'.sam'
			run_command(cmd)
			cmd=samtools+" view -@ "+str(thread)+" -bS "+tmp[0]+'.sam'+" -q "+str(MAPQ)+" -o "+ tmp[0]+'.bam'
			run_command(cmd)
			cmd='java -jar -XX:ParallelGCThreads='+str(thread)+' '+PICARD+'/SortSam.jar INPUT='+tmp[0]+'.bam'+ ' OUTPUT='+tmp[0]+'.sort.bam'+' SORT_ORDER=coordinate'
			run_command(cmd)
			cmd='java -jar -XX:ParallelGCThreads='+str(thread)+' '+PICARD+'/MarkDuplicates.jar INPUT='+tmp[0]+'.sort.bam'+' OUTPUT='+tmp[0]+'.sort.dedup.bam METRICS_FILE='+tmp[0]+'.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true REMOVE_DUPLICATES=true'
			run_command(cmd)
			tmpsam=tmp[0]+'.sam'
			tmpbam=tmp[0]+'.bam'
			tmpfix=tmp[0]+'.metrics'
			tmpfixsort=tmp[0]+'.sort.bam'
			os.remove(tmpsam)
			os.remove(tmpbam)
			os.remove(tmpfix)
			os.remove(tmpfixsort)
		
	
	cmd="awk '{print $1}' Female_marker.txt |sort >Female_mark.idtxt "
	run_command(cmd)
	cmd="awk '{print $1}' Male_marker.txt |sort >Male_mark.idtxt "
	run_command(cmd)
	cmd="sort Female_mark.idtxt Male_mark.idtxt  | uniq  > mergeid.txt "
	run_command(cmd)
	parentmyregion('mergeid.txt')
	df=pd.read_csv(progenyID,sep='\t',header=None)
	len_df=len(df)
	if (len_df<thread):
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
	p.map(callprogeny,order)
	p.close()
	p.join()
	for a in order:
		os.remove(a)
	
	with open(progenyID) as f:
		for line in f:
			tmp = line.strip().split('\t')
			cmd='python '+script+'/progenyGT.py '+'Female_marker.txt '+tmp[0]+'.all_type '+tmp[0]+'.out'
			run_command(cmd)
			tmout=tmp[0]+'.out'
			os.remove(tmout)
			cmd='python '+script+'/progenyGT.py '+'Male_marker.txt '+tmp[0]+'.all_type '+tmp[0]+'.out2'
			run_command(cmd)
			tmout=tmp[0]+'.out2'
			os.remove(tmout)
		
	list1 = []
	get_proganyID(list1,progenyID)
	genetype('Male_marker.txt',list1,progenyID,'Male_marker.out')
	genetype('Female_marker.txt',list1,progenyID,'Female_marker.out')
	
def parentmyregion(c):
    file=open("myregions.tmp","w")
    file.close
    file=open("regious.tmp","w")
    file.close
    cmd="sed 's/\t/:/g' myregions >myregions.tmp"
    run_command(cmd)
    with open(c,'r') as f:
        for line in f:
            tmp=line.strip().split("\t")
            cmd=' grep -w '+tmp[0]+' myregions.tmp' +' >>regious.tmp'
            run_command(cmd)
    cmd="less -S regious.tmp|sort |uniq|sed  's/:/\t/g' >parentmyregions"
    run_command(cmd)

def callprogeny(i):
    file=open(i,"r")
    lines = file.readlines()
    for line in lines:
        tmp=line.split("\t")
        gene=tmp[0]
        my_progeny(gene)
    file.close
def my_progeny(x):
    y=x+'.sort.dedup.bam'
    cmd=GangSTR +' --bam '+ y + ' --ref '+ GENOME +' --regions parentmyregions'+ ' --stutterup '+Stutter+' --stutterdown '+Stutter+' --out '+ x +' '+wgs
    run_command(cmd)
    vcf=x+'.vcf'
    outall=x+'.alltype'
    callprogeny_vcf(vcf,outall)
    cmd='less -S '+outall+' |sort|uniq > '+x+'.all_type'
    run_command(cmd)
    cmd=' rm '+outall
    run_command(cmd)
def callprogeny_vcf(x,outall):
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
                            myout=gene_rep+'\t'+mydtct0+'/'+mydtct0
                            with open(outall,"a+") as f3:
                                f3.write(myout+"\n")
                        elif(GT=='1/1'):
                            mydtct0=tmp[4]
                            if(len(mydtct0)>5):
                                myout=gene_rep+'\t'+mydtct0+'/'+mydtct0
                                with open(outall,"a+") as f3:
                                    f3.write(myout+"\n")
                        elif(GT=='0/1' or GT=='1/0'):
                            mydtct0=tmp[3]
                            mydtct2=tmp[4]
                            if(len(mydtct2)>5):
                                myout=gene_rep+'\t'+mydtct0+'/'+mydtct2
                                with open(outall,"a+") as f3:
                                    f3.write(myout+"\n")
                        elif(GT=='1/2'):
                            mydtct0=tmp[4].split(",")[0]
                            mydtct2=tmp[4].split(",")[1]
                            if(len(mydtct2)>5 and len(mydtct0)>5):
                                myout=gene_rep+'\t'+mydtct0+'/'+mydtct2
                                with open(outall,"a+") as f3:
                                    f3.write(myout+"\n")
                    else:
                        continue

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
    cmd='Rscript ' +script+'/abxcd.R '+out2+' '+myout2
    run_command(cmd)
    if(os.path.isfile(myout2)==True):
        file = open(myout2, "r")
        z=y+'.joinmap'
        lines = file.readlines()
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
    cmd='Rscript ' +script+'/aaxbc.R '+out2+' '+myout2
    run_command(cmd)
    if(os.path.isfile(myout2)==True):
        file = open(myout2, "r")
        z=y+'.joinmap'
        lines = file.readlines()
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
                 list1[i][1]='<efxeg>'
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
                        list2[i][7+j]='ee'
                    elif (a!='--' and (a==mya+"/"+myc or a==myc+"/"+mya) ):
                        list2[i][7+j]='eg'
                    elif (a!='--' and (a==myc+"/"+myb or a==myb+"/"+myc) ):
                        list2[i][7+j]='fg'
                    elif (a!='--' and (a==mya+"/"+myb or a==myb+"/"+mya) ):
                        list2[i][7+j]='fe'
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
                 list1[i][1]='<abxcd>'
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
                 list1[i][1]='<nnxnp>'
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
                        list2[i][7+j]='nn'
                    elif (a!='--' and a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[1] ) or (a.split("/")[0]==myab1.split("/")[1] and a.split("/")[1]==myab.split("/")[0]))):
                        list2[i][7+j]='np'
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
                 list1[i][1]='<lmxll>'
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
                        list2[i][7+j]='ll'
                    elif(a!='--' and a.split("/")[0]!=a.split("/")[1] and ((a.split("/")[0]==myab.split("/")[0] and a.split("/")[1]==myab1.split("/")[1] ) or (a.split("/")[0]==myab1.split("/")[1] and a.split("/")[1]==myab.split("/")[0] ))):
                        list2[i][7+j]='lm'
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
	cmd='Rscript ' +script+'/pchisq.R '+out2+' '+out3+' '+myout2+' '+myout3
	run_command(cmd)
def ALL_type(y):
	df1 = pd.read_csv(y,header=None,sep='\t')
	list1=df1.values.tolist()
	df1=pd.DataFrame(list1)
	df1 = df1.fillna('--')
	list2=df1.values.tolist()
	progany= pd.read_csv(progenyID,sep='\t')
	progany=progany.values.tolist()
	for i in range(len(list2)):
		ssr=''.join(re.findall(r'[atcg]', list2[i][1]))
		if (change(list2[i][3].split('/')[0],ssr) and change(list2[i][3].split('/')[1],ssr) and change(list2[i][5].split('/')[0],ssr) and change(list2[i][5].split('/')[1],ssr)):
			list2[i][3]=change(list2[i][3].split('/')[0],ssr)+'/'+change(list2[i][3].split('/')[1],ssr)
			list2[i][5]=change(list2[i][5].split('/')[0],ssr)+'/'+change(list2[i][5].split('/')[1],ssr)
		for j in range(len(progany)+1):            
			a=list2[i][7+j]
			if (a!='--'):
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
		df1.to_excel("./Male_alltype.xls", index=False)
	elif(y=='Female_marker.txt'):
		df1.to_csv('Female_marker.txt.1',sep='\t',header=False,index=False)
		df1.to_excel("./Female_alltype.xls", index=False)

def change(s,a):
    c=re.finditer(r'(%s)+'%(a),s)
    find_list=[]
    for i in c:
        find_list.append(i.group())
    if (len(find_list)==1):
        gt=a+'('+str(s.count(a))+')'
        return (gt)
    else:
        return (0)

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
    df1 = pd.read_csv(z,header=None,sep='\t')
    list1=df1.values.tolist()
    progany= pd.read_csv(progenyID,sep='\t')
    progany=progany.values.tolist()
    
    for i in range(len(list1)):
        list1[i][1]='<abxaa>'
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
    
    if(os.path.isfile(hybrid_pchis)==False):
        print("There were no data consistent with abxab segregation pattern")
    else:
        file = open(hybrid_pchis, "r")
        lines = file.readlines()
        for line in lines:
            line=line.strip()
            tmp=line.split("\t")
            if (float(tmp[5])>float(pvalue)):
                get(tmp[0],y,hybrid_map)
        file.close()
        if(os.path.isfile(hybrid_map)==False):
            print("There were no data consistent with abxab segregation pattern")
        else:
            abxab(hybrid_map,abxab_Fslinkmap)
    
def abxab(hybrid_map,abxab_Fslinkmap):
    with open(hybrid_map,"r") as f2:
        lines = f2.readlines()
        flen=len(lines)
        if (flen > 0):
            df1 = pd.read_csv(hybrid_map,header=None,sep='\t')
            list1=df1.values.tolist()
            for i in range(len(list1)):
                 list1[i][1]='<abxab>'
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
    cmd='cat aaxab.txt aaxbc.txt abxaa.txt abxab.txt abxac.txt abxcc.txt abxcd.txt >>allgenotype.txt'
    run_command(cmd)
    list1 = []
    get_proganyID(list1,progenyID)
    df1 = pd.read_csv('allgenotype.txt',header=None,sep='\t')
    df1=df1.values.tolist()
    df1=pd.DataFrame(df1)
    namelist=['SSR_site','Ref_SSR_type','Number','Female','Female_segregation_type','Male','Male_segregation_type',]
    list1=namelist+list1
    df1.columns=list1
    df1.to_excel('allgenotype.xls', index=False,header=True)
    chang_head('allgenotype.txt')
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
def call(x):
	y=x+'.sort.dedup.bam'
	cmd=GangSTR +' --bam '+ y + ' --ref '+ GENOME +' --regions myregions'+ ' --stutterup '+Stutter+' --stutterdown '+Stutter+' --out '+ x +' '+wgs 
	run_command(cmd)
	vcf=x+'.vcf'
	out=x+'.vcf'+'01.txt'
	out2=x+'allSSR_type.txt'
	cmd='python '+script+'/parentGT.py'+ ' ' +vcf + '  '+out +' '+out2
	run_command(cmd)
	myout=out+'.out'
	cmd='less -S '+out+ ' |sort|uniq > '+myout
	run_command(cmd)
def callparent(order):
	p=Pool(2)
	p.map(call,order)
	p.close()
	p.join()
def getfastiq(x):
    
    Myproid=open('Myprogeny.id','w')
    Myproid.close()
    with open(x,'r') as f:
        i=0
        for line in f:
            i=i+1
            if (line.startswith("[progeny fastq files]")):
                flag=i
                break
    with open(x,'r') as f:
        f =f.readlines()[flag:]
        for line in f:
            line=line.strip()
            myout=line.split(":")[0]+'\t'+line.split()[0].split(":")[1]+'\t'+line.split()[1]
            with open('Myprogeny.id',"a+") as f2:
                f2.write(myout+"\n")
    
    if(os.path.isfile('Myprogeny.id')):
        return ('./Myprogeny.id')
def main(args):
	global progenyID,motif,wgs 
	progenyID=getfastiq("parameters.ini")
	motif = args.motif
	wgs = args.WGS
	if (wgs=='True'):
		wgs=''
	elif(wgs=='False'):
		wgs=' --targeted --coverage  --nonuniform '
	get_SSR(GENOME)
	find_parent(GENOME)
	tasks=['male','female']
	callparent(tasks)
	SSRGM(GENOME)
	pchisq('Male_marker.out','Male_pure_pchis','Male_hybrid_pchis')
	pchisq('Female_marker.out','Female_pure_pchis','Female_hybrid_pchis')
	ALL_type('Male_marker.txt')
	ALL_type('Female_marker.txt')
	get_segregation('Female_marker.txt.1','Male_marker.txt.1')
	write2map('Male_pure_pchis.out','Male_marker.txt','aaxab.joinmap.txt','Male_hybrid_pchis.out','Male.abxab.joinmap.txt','Male_abxaa_Fslinkmap.txt','Male_abxab_Fslinkmap.txt')
	write2map('Female_pure_pchis.out','Female_marker.txt','abxaa.joinmap.txt','Female_hybrid_pchis.out','Female.abxab.joinmap.txt','Female_abxaa_Fslinkmap.txt','Female_abxab_Fslinkmap.txt')
	os.mkdir('./WorkingDirectory')	
	cmd='mv *.bam *.bai ./WorkingDirectory'
	run_command(cmd)
	cmd='mv female* male*  Male* Female*   ./WorkingDirectory'
	run_command(cmd)
	cmd='mv  *.vcf *marker* ./WorkingDirectory'
	run_command(cmd)	
	cmd='rm *.all_type  *.tab'
	run_command(cmd)
if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		prog="python SSRGT.py",
		formatter_class=argparse.RawTextHelpFormatter,
		description='''
-------------------------------------------------------------------------------------------------------
SSRGT is a software used to genotype and SSR calling across a hybrid population with genome resequencing data.			
Contact: Chunfa Tong  <tongchf@njfu.edu.cn>
Version: 1.0
Usage: python SSRGT.py  [Options]

-------------------------------------------------------------------------------------------------------
''')
	parser.add_argument(
		'-mo', '--motif',
		default="1=10,2=5,3=4,4=4,5=4,6=4",
		help="set the threshold of motif, default: 1=10,2=5,3=4,4=4,5=4,6=4 \nleft  of equal : length of motif \nright of equal : the minimum number of repeat")
	parser.add_argument(
		'-wgs', '--WGS', 
		default='False',
		help="set the sequencing data type, default : False \nIf your sequencing data is RAD-seq or WES data, you should choose '-wgs False', if it's whole genome sequencing, choose '-wgs True'")
	args = parser.parse_args()
	filename='./parameters.ini'
	if(os.path.isfile(filename)==False):
        	print("\n This program needs a parameter file,namely \"parameters.ini\".\n")
        	sys.exit(' Error!Please check the parameters.ini')
	if(os.path.isfile(GENOME)==False):
        	sys.exit(' Error!Please check  the line of REFERENCE_GENOME in the parameters.ini')
	if(os.path.isfile(bwa)==False):
        	sys.exit(' Error!Please check  the line of BWA in the parameters.ini')
	if(os.path.isfile(samtools)==False):
        	sys.exit(' Error!Please check  the line of SAMTOOLS in the parameters.ini')
	if(os.path.isfile(male1)==False ):
        	sys.exit(' Error!Please check  the line of MALEPARENT_1 in the parameters.ini')
	if(os.path.isfile(male2)==False ):
        	sys.exit(' Error!Please check  the line of MALEPARENT_2 in the parameters.ini')
	if(os.path.isfile(female1)==False ):
        	sys.exit(' Error!Please check  the line of FEMALEPARENT_1 in the parameters.ini')
	if(os.path.isfile(female2)==False ):
        	sys.exit(' Error!Please check  the line of FEMALEPARENT_2 in the parameters.ini')
	if(os.path.exists(progenyfold)):
		if(os.path.exists(script)):
			main(args)
			print("All tasks done")
		else:
			print(' Error!Please check  the line of SCRIPT in the parameters.ini')
	else:
		print(' Error!Please check  the line of PROGENY in the parameters.ini')
