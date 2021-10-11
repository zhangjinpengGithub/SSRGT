import sys,os
import pandas as pd
from scipy import stats
x=sys.argv[1]
y=sys.argv[2]
df1 = pd.read_csv(x,header=None,sep='\t')
df1.insert(3,'3',[l for l in range(len(df1))])
df1.insert(4,'4',[k for k in range(len(df1))])
list2=df1.values.tolist()
for i in range(len(list2)):
    list2[i][3]=(int(list2[i][1])+int(list2[i][2]))/2
    e1=list2[i][3]
    a=int(list2[i][1])
    b=int(list2[i][2])
    list2[i][4]=stats.chisquare([a,b], f_exp = [e1,e1]).pvalue
df1=pd.DataFrame(list2)
df1.to_csv(y,sep='\t',index=False,header=False)
