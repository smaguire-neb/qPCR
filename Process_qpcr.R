# melt curve summary
# analysis summary

library(dplyr)
library(readxl)
library(stringr)
setwd("~/Desktop/qpcr")

layout1<-read_excel("layout-plate1.xlsx")
base.dir<-"~/Desktop/qpcr/"

process.plate<-function(dir,plate.type = c("plate1","plate2")){
  dir.i<-paste0(base.dir,dir)
  cq.file<-list.files(dir.i,pattern="Quantification Cq Results -reformatted") %>%
    .[!str_detect(.,"~")] %>%
    paste0(base.dir,dir,"/",.)
  melt.file<-list.files(dir.i,pattern="Melt Curve Summary -reformatted")  %>%
    .[!str_detect(.,"~")] %>%
    paste0(base.dir,dir,"/",.)
  base.name<-list.files(dir.i,pattern="Quantification Cq Results") %>%
    .[!str_detect(.,"~")] %>%
    str_split(.,pattern=" - ",simplify = T) %>%
    .[1,1]
  cq.data<-read_excel(cq.file) %>%
    select(Well,Cq)
  melt.data<-read_excel(melt.file) %>%
    group_by(Well) %>%
    mutate(count=n(),
           melt.flag=ifelse(count == 1,F,T)) %>%
    mutate(melt.temp = suppressWarnings(as.numeric(`Melt Temp`))) %>%
    group_by(melt.flag,Well) %>%
    summarise(melt.temp=mean(melt.temp)) %>%
    as.data.frame()
  quant<-left_join(cq.data,melt.data)
  if(plate.type == "plate1"){
    layout <- read_excel("layout-plate1.xlsx") %>%
      mutate(Well = well) %>%
      select(miRname,Well,sequence)
  } else if(plate.type == "plate2") {
    layout <- read_excel("layout-plate2.xlsx")  %>%
      mutate(Well = well) %>%
      select(miRname,Well,sequence)
  }
  full_join(quant,layout)
}

results1<-process.plate("11-5-18","plate1")
results2<-process.plate("11-12-18-plate1-rep2","plate1")
r1<-
  results1 %>%
  filter( Cq < 37) %>%
  filter(!(is.na(melt.temp) & Cq > 35))

r2<-
  results2 %>%
  filter( Cq < 37) %>%
  filter(!(is.na(melt.temp) & Cq > 35))


left_join(r1,r2,by=c("Well","miRname","sequence")) %>%
  ggplot(aes(x=Cq.x,y=Cq.y))+geom_point() + xlab("Cq rep1") + 
  ylab("Cq rep 2")

qpcr.data.plate1<-bind_rows(results1,results2) %>%
  group_by(miRname) %>%
  mutate(mean.cq = mean(Cq),
         sd.cq=sd(Cq),
         cv=sd.cq/mean.cq) %>%
  filter(cv <= 0.05,
         mean.cq<37)
#Take a list of processed plates and convert them to qpcr data
#check the CVs. If they are bad eventually figure out some way to toss the outliers. 
#Toss the outliers if they fail the CV check then find the mean of all replicates, subtract the mean, toss the one with the highest abs. value. 


