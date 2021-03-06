import sys,os
import pandas as pd
from scipy import stats
pure_pchis=sys.argv[1]
hybrid_pchis=sys.argv[2]
y=sys.argv[3]
z=sys.argv[4]
if(os.path.isfile(pure_pchis)):
    df1 = pd.read_csv(pure_pchis,header=None,sep='\t')
    df1.insert(3,'3',[l for l in range(len(df1))])
    df1.insert(4,'4',[k for k in range(len(df1))])
    df1.insert(5,'5',[k for k in range(len(df1))])
    list2=df1.values.tolist()
    for i in range(len(list2)):
        list2[i][3]=(int(list2[i][1])+int(list2[i][2]))/2
        e1=list2[i][3]
        a=int(list2[i][1])
        b=int(list2[i][2])
        list2[i][4]=stats.chisquare([a,b], f_exp = [e1,e1]).pvalue
        list2[i][5]=stats.chisquare([a,b], f_exp = [e1,e1]).statistic
    df1=pd.DataFrame(list2)
    df1.to_csv(y,sep='\t',index=False,header=False)
else:
    print("No SSR genotypes are found that fit Mendelian segregation ratios for abxaa or aaxab.")
with open(hybrid_pchis,"r") as f2:
    lines = f2.readlines()
    flen=len(lines)
    if (flen > 2):
        df1 = pd.read_csv(hybrid_pchis,header=None,sep='\t')
        df1.insert(4,'4',[l for l in range(len(df1))])
        df1.insert(5,'5',[k for k in range(len(df1))])
        df1.insert(6,'6',[k for k in range(len(df1))])
        list2=df1.values.tolist()
        for i in range(len(list2)):
            list2[i][4]=((int(list2[i][1])+int(list2[i][2]))+int(list2[i][3]))/4
            e1=list2[i][4]
            a=int(list2[i][1])
            b=int(list2[i][2])
            c=int(list2[i][3])
            list2[i][5]=stats.chisquare([a,b,c], f_exp = [e1,2*e1,e1]).pvalue
            list2[i][6]=stats.chisquare([a,b,c], f_exp = [e1,2*e1,e1]).statistic
        df1=pd.DataFrame(list2)
        df1.to_csv(z,sep='\t',index=False,header=False)
#    else:
#        print("No SSR genotypes are found that fit Mendelian segregation ratios for abxab.")


