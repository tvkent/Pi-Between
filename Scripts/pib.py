import pandas as pd
import numpy as np
# allo vs indica
for n in range(1,13):
	input1="/home/tvkent/projects/rice/angsd-wrapper/results/alloc"+str(n)+"_Diversity.mafs.gz"
	input3="/home/tvkent/projects/rice/angsd-wrapper/results/indic"+str(n)+"_Diversity.mafs.gz"
	out1="/home/tvkent/ai_"+str(n)+".txt"
	
	a = pd.read_table(input1, compression='gzip', usecols=['chromo','position','major','knownEM','nInd'])
	a=a[a.nInd == 4]
	a['a.maf']=a['knownEM']
	a['a.major']=a['major']
	a=a.drop('knownEM',1)
	a=a.drop('major',1)
	a['a.majorf']=1-a['a.maf']
	
	i = pd.read_table(input3, compression='gzip', usecols=['chromo','position','major','knownEM','nInd'])
	i=i[i.nInd >= 33]
	i['i.maf']=i['knownEM']
	i['i.major']=i['major']
	i=i.drop('knownEM',1)
	i=i.drop('major',1)
	i['i.majorf']=1-i['i.maf']
	
	ai = pd.merge(a,i,on=['chromo','position'],how='inner')
	ai=ai.drop(['nInd_x','nInd_y'], 1)
	ai['pib']= np.where(ai['a.major']==ai['i.major'],(ai['a.maf']*ai['i.majorf'] + ai['a.majorf']*ai['i.maf']),(ai['a.maf']*ai['i.maf'] + ai['a.majorf']*ai['i.majorf']))
	ai['position']=ai['position']/1000000
	ai['win']=np.floor(ai['position'])
	group=ai.groupby('win')
	aidiv=group['pib'].agg([np.sum,'count'])
	aidiv['div_perbp']=aidiv['sum']/aidiv['count']
	aidiv['chromo']=n
	aidiv.to_csv(out1, sep='\t')

# symp vs indica
for n in range(1,13):
        input2="/home/tvkent/projects/rice/angsd-wrapper/results/sympc"+str(n)+"_Diversity.mafs.gz"
        input3="/home/tvkent/projects/rice/angsd-wrapper/results/indic"+str(n)+"_Diversity.mafs.gz"
        out2="/home/tvkent/si_"+str(n)+".txt"
	
        s = pd.read_table(input2, compression='gzip', usecols=['chromo','position','major','knownEM','nInd'])
        s=s[s.nInd == 4]
        s['s.maf']=s['knownEM']
        s['s.major']=s['major']
        s=s.drop('major',1)
        s=s.drop('knownEM',1)
        s['s.majorf']=1-s['s.maf']
	
        i = pd.read_table(input3, compression='gzip', usecols=['chromo','position','major','knownEM','nInd'])
        i=i[i.nInd >= 33]
        i['i.maf']=i['knownEM']
        i['i.major']=i['major']
        i=i.drop('knownEM',1)
        i=i.drop('major',1)
        i['i.majorf']=1-i['i.maf']
	
        si = pd.merge(s,i,on=['chromo','position'],how='inner')
        si=si.drop(['nInd_x','nInd_y'], 1)
        si['pib']=np.where(si['s.major']==si['i.major'],(si['s.maf']*si['i.majorf'] + si['s.majorf']*si['i.maf']),(si['s.maf']*si['i.maf'] + si['s.majorf']*si['i.majorf']))
        si['position']=si['position']/1000000
        si['win']=np.floor(si['position'])
        group=si.groupby('win')
        sidiv=group['pib'].agg([np.sum,'count'])
        sidiv['div_perbp']=sidiv['sum']/sidiv['count']
        sidiv['chromo']=n
        sidiv.to_csv(out2, sep='\t')




