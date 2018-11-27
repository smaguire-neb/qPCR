# melt curve summary
# analysis summary

library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(purrr)

setwd("~/Desktop/qpcr/Qpcr code/data")

process.plate<-function(data,plate.type){
  if(plate.type == "plate1"){
    layout <- read_excel("layout-plate1.xlsx") %>%
      mutate(Well = well) %>%
      select(miRname,Well,sequence)
  } else if(plate.type == "plate2") {
    layout <- read_excel("layout-plate2.xlsx")  %>%
      mutate(Well = well) %>%
      select(miRname,Well,sequence)
  }
  full_join(data,layout,by="Well") %>%
    mutate(Cq = as.numeric(as.character(Cq)),
           Melt.Temperature = as.numeric(as.character(Melt.Temperature))) %>%
    mutate(flag = ifelse((Cq < 37) | (!(is.na(Melt.Temperature) & Cq > 35)),F,T))
}

files<-list.files(pattern=".txt")
data<-
as_tibble(str_split(files,"-",simplify = T)) %>%
  select(-V3) %>%
  rename(plate.type = V1, replicate = V2) %>%
  mutate(file.names = files) %>%
  group_by(plate.type,replicate) %>%
  nest(.key="file.name") %>%
  mutate(raw.data = map(.$file.name, function(x) unlist(x) %>%
                      read.delim(.,skip=19))) %>%
  mutate(processed.data = map2(.$raw.data,.$plate.type, process.plate))

names(data$processed.data)<-paste(data$plate.type,data$replicate,sep="-")

long.data<-bind_rows(data$processed.data,.id="file")

long.data %>%
  mutate(plate = str_split(file,"-",simplify = T)[,1],
         rep = str_split(file,"-",simplify = T)[,2]) %>%
  filter(miRname == "UniSp3 IPC")


# data corrected by the IPC file.
corrected.long.data<-
mutate(long.data,plate = str_split(file,"-",simplify = T)[,1],
           rep = str_split(file,"-",simplify = T)[,2]) %>%
  filter(miRname == "UniSp3 IPC") %>%
  group_by(file) %>%
    mutate(ipc.replicate = mean(Cq),
           sd.replicate = sd(Cq),
           cv.replicate = (sd.replicate/ipc.replicate)*100) %>%
  group_by() %>%
  mutate(diff.mean = abs(Cq-ipc.replicate)) %>%
  group_by(file) %>%
  mutate(possible.outlier = ifelse(diff.mean == max(diff.mean),T,F)) %>%
  mutate(outlier.flag = ifelse(cv.replicate >1 & possible.outlier,T,F)) %>%
  filter(!(outlier.flag)) %>%
  group_by(file) %>%
  mutate(ipc.replicate = mean(Cq),
         sd.replicate = sd(Cq),
         cv.replicate = (sd.replicate/ipc.replicate)*100) %>%
  group_by(plate) %>%
  mutate(ipc.overall = mean(Cq)) %>%
  group_by(file) %>%
  summarise(ipc.cf = mean(ipc.replicate - ipc.overall)) %>%
  left_join(long.data,.) %>%
  mutate(corrected.cq = Cq-ipc.cf)

corrected.long.data %>%
  mutate(plate = str_split(file,"-",simplify = T)[,1],
         rep = str_split(file,"-",simplify = T)[,2]) %>%
  filter(!is.na(sequence)) %>%
  group_by(Well,sequence,miRname) %>%
  mutate(mean.diff=abs(corrected.cq-mean(corrected.cq)),
         outlier.test = max(mean.diff),
         outlier = ifelse(mean.diff == outlier.test,T,F)) %>%
  filter(!outlier) %>%
  summarise(mean.corrected.cq = mean(corrected.cq),
            sd.corrected.cq = sd(corrected.cq),
            cv.corrected.cq = (sd.corrected.cq/mean.corrected.cq)*100) %>%
  #716 genes left
  #filter(mean.corrected.cq < 37, cv.corrected.cq < 5) %>%
  filter(cv.corrected.cq < 5) %>%
  ggplot(aes(x=mean.corrected.cq))+geom_density()
  

#filter them so the just contain flag, cq, well, name and sequence
#make the data long
#group the data by name and sequence
#calculate the cv for each one
#eventually when I have 4 datapoints for each I would like to automatically remove
#the outlier datapoints. IE choose the set of 3 points that has the lowest CV. 

#Take a list of processed plates and convert them to qpcr data
#check the CVs. If they are bad eventually figure out some way to toss the outliers. 
#Toss the outliers if they fail the CV check then find the mean of all replicates, subtract the mean, toss the one with the highest abs. value. 


