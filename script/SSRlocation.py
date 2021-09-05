#!/usr/bin/python
import sys, os,argparse,re
import subprocess
import pandas as pd
from multiprocessing import Pool
parser = argparse.ArgumentParser(
     description='''
Example: python SSRlocation.py -g ../p.gff3 -s ../femaleallSSR_type.txt -t 2 -o female.location
''')
parser.add_argument('-g','--gff',help='Please input GFF file',required=True)
parser.add_argument('-s','--SSR_types',help='Please input SSR_types.txt file',required=True)
parser.add_argument('-t','--threads',help='Please input threads',required=True)
parser.add_argument('-o','--out',help='Please input out_put file_name',required=True)

args = parser.parse_args()
def run_command(cmd):
        print(cmd)
        return_code = subprocess.call(cmd, shell=True)
def location(Chr,genestart):
	with open(args.gff,'r') as f:
		i=0
		for line in f:
			if line.startswith("#"):
				continue
			else:
				tmp=line.split("\t")
				gffnum=int(''.join(re.findall(r'\d+',tmp[0])))
				gffstr=re.findall(r'\D+',tmp[0])
				chrstr=re.findall(r'\D+',Chr)
				chrnum=int(''.join(re.findall(r'\d+',Chr)))
				if (tmp[0]==Chr and genestart>int(tmp[3]) and genestart<int(tmp[4]) ):
					myout=tmp[0]+'\t'+str(genestart)+'\t'+tmp[2]+'\t'+tmp[8].split(";")[0].split("=")[1]+'\t'+tmp[3]+'\t'+tmp[4]
					with open(args.out,'a') as f3:
						f3.write(myout+"\n")
						i=i+1
						if(i==4):
							break
					
def callSSR(i):
    file=open(i,"r")
    lines = file.readlines()
    for line in lines:
        strings = line.split("\t")[0]
        Chr=strings.split(":")[0]
        genestart=int(strings.split(":")[1])
        location(Chr,genestart)
def main():
    df=pd.read_csv(args.SSR_types,sep='\t',header=None)
    len_df=len(df)
    numble=int(len_df/int(args.threads)) 
    gtrd=pd.read_csv(args.SSR_types,chunksize=numble,sep='\t',header=None)
    a=0
    order=[]
    for i in gtrd :
        a=a+1
        order.append(str(a)+'.tmptxt')
        i.to_csv(f"{str(a)}.tmptxt",sep='\t',header=None,index=False)
    p=Pool(int(args.threads))
    p.map(callSSR,order)
    p.close()
    p.join()
    for a in order:
        os.remove(a)	
    cmd='less -S '+args.out+' |sort|uniq > '+args.out+'.sort.txt'
    run_command(cmd)
if __name__ == '__main__' :
    main()
