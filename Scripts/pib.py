import pandas as pd
import numpy as np

# allo vs indica
for n in range(1,13):#alloc9_intergenicDiversity.mafs indicach6_intergenicDiversity.filtered.mafs alloch11_intergenicDiversity.filtered.mafs
	input1="~/rice/alloc"+str(n)+"_5kbintergenic4foldDiversity.mafs"
	input3="~/rice/indicac"+str(n)+"_5kbintergenic4foldDiversity.mafs"
	out1="~/rice/ai_"+str(n)+"_1kb_5kbintergenic4fold.txt"
	
	a = pd.read_table(input1,usecols=['chromo','position','major','minor','knownEM','nInd'])
	a=a[a.nInd == 4]
	a['a.maf']=a['knownEM']
	a['a.major']=a['major']
	a['a.minor']=a['minor']
	a=a.drop('knownEM',1)
	a=a.drop('major',1)
	a['a.majorf']=1-a['a.maf']
	
	i = pd.read_table(input3, usecols=['chromo','position','major','minor','knownEM','nInd'])
	i=i[i.nInd >= 34]
	i['i.maf']=i['knownEM']
	i['i.major']=i['major']
	i['i.minor']=i['minor']
	i=i.drop('knownEM',1)
	i=i.drop('major',1)
	i['i.majorf']=1-i['i.maf']
	
	ai = pd.merge(a,i,on=['chromo','position'],how='inner')
	ai=ai.drop(['nInd_x','nInd_y'], 1)
	ai['keep'] = ai.apply(lambda x: ((x['a.major'] in x['i.major']) and (x['a.minor'] in x['i.minor'])) or ((x['a.major'] in x['i.minor']) and (x['a.minor'] in x['i.major'])), axis=1)
	ai=ai[ai['keep']==True]
	ai['pib']= np.where(ai['a.major']==ai['i.major'],(ai['a.maf']*ai['i.majorf'] + ai['a.majorf']*ai['i.maf']),(ai['a.maf']*ai['i.maf'] + ai['a.majorf']*ai['i.majorf']))
	ai['position']=ai['position']/1000
	ai['win']=np.floor(ai['position'])
	group=ai.groupby('win')
	aidiv=group['pib'].agg([np.sum,'count'])
	aidiv['div_perbp']=aidiv['sum']/aidiv['count']
	aidiv['chromo']=n
	aidiv.to_csv(out1, sep='\t')

# symp vs indica
for n in range(1,13):
	input2="~/rice/sympc"+str(n)+"_5kbintergenic4foldDiversity.mafs"
	input3="~/rice/indicac"+str(n)+"_5kbintergenic4foldDiversity.mafs"
	out2="~/rice/si_"+str(n)+"_1kb_5kbintergenic4fold.txt"
	
	s = pd.read_table(input2, usecols=['chromo','position','major','knownEM','nInd'])
	s=s[s.nInd == 4]
	s['s.maf']=s['knownEM']
	s['s.major']=s['major']
	s=s.drop('major',1)
	s=s.drop('knownEM',1)
	s['s.majorf']=1-s['s.maf']
	
	i = pd.read_table(input3,  usecols=['chromo','position','major','knownEM','nInd'])
	i=i[i.nInd >= 34]
	i['i.maf']=i['knownEM']
	i['i.major']=i['major']
	i=i.drop('knownEM',1)
	i=i.drop('major',1)
	i['i.majorf']=1-i['i.maf']
	
	si = pd.merge(s,i,on=['chromo','position'],how='inner')
	si=si.drop(['nInd_x','nInd_y'], 1)
	si['keep'] = si.apply(lambda x: ((x['s.major'] in x['i.major']) and (x['s.minor'] in x['i.minor'])) or ((x['s.major'] in x['i.minor']) and (x['s.minor'] in x['i.major'])), axis=1)
	si=si[si['keep']==True]
	si['pib']=np.where(si['s.major']==si['i.major'],(si['s.maf']*si['i.majorf'] + si['s.majorf']*si['i.maf']),(si['s.maf']*si['i.maf'] + si['s.majorf']*si['i.majorf']))
	si['position']=si['position']/1000
	si['win']=np.floor(si['position'])
	group=si.groupby('win')
	sidiv=group['pib'].agg([np.sum,'count'])
	sidiv['div_perbp']=sidiv['sum']/sidiv['count']
	sidiv['chromo']=n
	sidiv.to_csv(out2, sep='\t')




