setwd('~/Pi-Between/Results/')
library(magrittr)
library(dplyr)
library(ggplot2)
library(data.table)

load("raw_genetic_maps.RData") #gets osat object = O. Sativa genetic map
rec<-select(osat,chr,cm,mb) %>% arrange(chr,mb,cm) %>% mutate(win=floor(mb)) %>% group_by(chr,win) %>%
  summarise(start_mb=min(mb),end_mb=max(mb),start_cm=min(cm),end_cm=max(cm),count=n()) %>%
  mutate(rec=(end_cm-start_cm)/(end_mb-start_mb)) %>%filter(count>1) %>%
  select(chr,win,rec) %>% setnames('chr','chromo')

ai<-fread("ai.txt",header=T)
si<-fread('si.txt',header=T)
adat<-merge(rec,ai,by=c("chromo","win"))
ggplot(adat)+geom_point(aes(x=rec,y=div_perbp))+
  geom_smooth(aes(x=rec,y=div_perbp),method='lm')+xlab("cM per Mb")+
  ylab(expression(paste("weighted ",pi[B],sep="")))+
  theme_bw()+ggtitle("Allopatric vs Indica")
cor.test(adat$rec,adat$div_perbp,method="spearman")
adat<-merge(rec,si,by=c("chromo","win"))
ggplot(adat)+geom_point(aes(x=rec,y=div_perbp))+
  geom_smooth(aes(x=rec,y=div_perbp),method='lm')+xlab("cM per Mb")+
  ylab(expression(paste("weighted ",pi[B],sep="")))+
  theme_bw()+ggtitle("Sympatric vs Indica")
cor.test(adat$rec,adat$div_perbp,method="spearman")
