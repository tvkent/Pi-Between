setwd('~/Rice_Results/PiB')
library(magrittr)
library(dplyr)
library(ggplot2)
library(data.table)

load("raw_genetic_maps.RData") #gets osat object = O. Sativa genetic map
rec<-select(osat,chr,cm,mb) %>% arrange(chr,mb,cm) %>% mutate(win=floor(mb)) %>%
  group_by(chr,win) %>% summarise(start_mb=min(mb),end_mb=max(mb),start_cm=min(cm),end_cm=max(cm),count=n()) %>%
  mutate(rec=(end_cm-start_cm)/(end_mb-start_mb)) %>% filter(count>1) %>% select(chr,win,rec) %>% setnames('chr','chromo')

x<-fread("alloGL.mafs",select=c("chromo","position","minor","unknownEM","nInd"),header=T)
allo<-filter(x, nInd==4) %>% setnames("unknownEM","maf") %>% mutate(majorf=1-maf) %>% select(chromo,position,minor,maf,majorf)
y<-fread("sympGL.mafs",select=c("chromo","position","minor","unknownEM","nInd"),header=T)
symp<-filter(y, nInd==4) %>% setnames("unknownEM", "maf") %>% mutate(majorf=1-maf) %>% select(chromo,position,minor,maf,majorf,nInd)
z<-fread("indiGL2.mafs",select=c("chromo","position","minor","unknownEM","nInd"),header=T)
indi<-filter(z, nInd==38) %>% setnames("unknownEM","maf") %>% mutate(majorf=1-maf) %>% select(chromo,position,minor,maf,majorf)
setkey(allo, by = 'chromo','position')
setkey(symp, by = 'chromo','position')
setkey(indi, by = 'chromo','position')
#allo vs symp
as<-allo[symp,] %>% na.omit() %>%
  mutate(pib=ifelse(minor==i.minor,(majorf*i.maf)+(i.majorf*maf),(majorf*i.majorf)+(i.maf*maf)),mb=position/1E6,win=floor(mb)) %>%
  group_by(chromo,win) %>% summarise(div=sum(pib),startmb=min(mb),endmb=max(mb),count=n()) %>% mutate(div_perbp=div/count)
adat<-merge(rec,as,by=c("chromo","win"))
ggplot(adat)+geom_point(aes(x=rec,y=div_perbp))+geom_smooth(aes(x=rec,y=div_perbp))+
  xlab("cM per Mb")+ylab(expression(paste("weighted ",pi[B],sep="")))+
  theme_bw()+ggtitle("Allopatric vs Sympatric")
cor.test(adat$rec,adat$div_perbp,method="spearman")
#allo vs indi
ai<-allo[indi,] %>% na.omit() %>%
  mutate(pib=ifelse(minor==i.minor,(majorf*i.maf)+(i.majorf*maf),(majorf*i.majorf)+(i.maf*maf)),mb=position/1E6,win=floor(mb)) %>%
  group_by(chromo,win) %>% summarise(div=sum(pib),startmb=min(mb),endmb=max(mb),count=n()) %>% mutate(div_perbp=div/count)
adat<-merge(rec,ai,by=c("chromo","win"))
ggplot(adat)+geom_point(aes(x=rec,y=div_perbp))+geom_smooth(aes(x=rec,y=div_perbp))+
  xlab("cM per Mb")+ylab(expression(paste("weighted ",pi[B],sep="")))+
  theme_bw()+ggtitle("Allopatric vs Indica")
cor.test(adat$rec,adat$div_perbp,method="spearman")
#symp vs indi
si<-symp[indi,] %>% na.omit() %>%
  mutate(pib=ifelse(minor==i.minor,(majorf*i.maf)+(i.majorf*maf),(majorf*i.majorf)+(i.maf*maf)),mb=position/1E6,win=floor(mb)) %>%
  group_by(chromo,win) %>% summarise(div=sum(pib),startmb=min(mb),endmb=max(mb),count=n()) %>% mutate(div_perbp=div/count)
adat<-merge(rec,si,by=c("chromo","win"))
ggplot(adat)+geom_point(aes(x=rec,y=div_perbp))+geom_smooth(aes(x=rec,y=div_perbp))+
  xlab("cM per Mb")+ylab(expression(paste("weighted ",pi[B],sep="")))+
  theme_bw()+ggtitle("Sympatric vs Indica")
cor.test(adat$rec,adat$div_perbp,method="spearman")
