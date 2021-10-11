import sys,os
import pandas as pd
from scipy import stats
x=sys.argv[1]
y=sys.argv[2]
df1 = pd.read_csv(x,header=None,sep='\t')
df1.insert(5,'5',[l for l in range(len(df1))])
df1.insert(6,'6',[h for h in range(len(df1))])
list2=df1.values.tolist()
for i in range(len(list2)):
    list2[i][5]=(int(list2[i][1])+int(list2[i][2])+int(list2[i][3])+int(list2[i][4]))/4
    e1=list2[i][5]
    a=int(list2[i][1])
    b=int(list2[i][2])
    c=int(list2[i][3])
    d=int(list2[i][4])
    list2[i][6]=stats.chisquare([a,b,c,d], f_exp = [e1,e1,e1,e1]).pvalue
df1=pd.DataFrame(list2)
df1.to_csv(y,sep='\t',index=False,header=False)

