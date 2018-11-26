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
         rep = str_split(file,"-",simplify = T)[,2])

left_join(data$processed.data[[2]],data$processed.data[[3]],by="Well") %>%
  ggplot(aes(x=Cq.x,y=Cq.y))+geom_point()


#Take a list of processed plates and convert them to qpcr data
#check the CVs. If they are bad eventually figure out some way to toss the outliers. 
#Toss the outliers if they fail the CV check then find the mean of all replicates, subtract the mean, toss the one with the highest abs. value. 


