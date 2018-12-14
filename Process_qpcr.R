# melt curve summary
# analysis summary

library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
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

corrected.summ.data<-
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
  filter(mean.corrected.cq < 37, cv.corrected.cq < 5)
  

#filter them so the just contain flag, cq, well, name and sequence
#make the data long
#group the data by name and sequence
#calculate the cv for each one
#eventually when I have 4 datapoints for each I would like to automatically remove
#the outlier datapoints. IE choose the set of 3 points that has the lowest CV. 

#Take a list of processed plates and convert them to qpcr data
#check the CVs. If they are bad eventually figure out some way to toss the outliers. 
#Toss the outliers if they fail the CV check then find the mean of all replicates, subtract the mean, toss the one with the highest abs. value. 

#splint.full.data<-read.delim("splint.full.tabular",header=F,skip=12)
splint.count.data<-
splint.full.data %>%
  group_by(V3) %>%
  summarise(splint.count=n()) %>%
  rename("miRname" = V3) %>%
  mutate(miRname = as.character(miRname))

#NN.full.data<-read.delim("NN.full.tabular",header=F,skip=12)
NN.count.data<-
  NN.full.data %>%
  group_by(V3) %>%
  summarise(NN.count=n()) %>%
  rename("miRname" = V3) %>%
  mutate(miRname = as.character(miRname))

layout1 <- read_excel("layout-plate1.xlsx")
layout2 <- read_excel("layout-plate2.xlsx")

mirNames<-
rbind(layout1,layout2) %>%
  filter(sequence != "",
         !is.na(sequence))

mirNames$miRname

length(intersect(as.character(splint.count.data$miRname),mirNames$miRname))
#464 present in splint dataset and miRnome dataset 
length(intersect(as.character(NN.count.data$miRname),mirNames$miRname))
#437

joined.data<-
inner_join(corrected.summ.data,splint.count.data) %>%
  inner_join(.,NN.count.data) %>%
  mutate(log.cq = log(mean.corrected.cq),
         log.splint.count = log(splint.count),
         log.NN.count = log(NN.count))

nrow(joined.data)

#690 detected by qpcr
# of those 452 were detected by Splint
# and 428 were detected by NN

454-428
nrow(inner_join(corrected.summ.data,splint.count.data))
nrow(inner_join(corrected.summ.data,NN.count.data))

splint.qpcr.seq<-
ggplot(joined.data,aes(x=log.cq,y=log.splint.count))+geom_point() + 
  geom_smooth(method="lm") +
  xlab("qPCR Cq (Log Scale)") +
  ylab("RNAseq Read Count (Log Scale)") +
  ggtitle("Splint Method")
NN.qpcr.seq<-
ggplot(joined.data,aes(x=log.cq,y=log.NN.count))+geom_point() + 
  geom_smooth(method="lm") +
  xlab("qPCR Cq (Log Scale)") +
  ylab("RNAseq Read Count (Log Scale)") +
  ggtitle("NEBNext")

summary(lm(log.cq~log.splint.count,data=joined.data))
summary(lm(log.cq~log.NN.count,data=joined.data))

save_plot("NNseqQ.png",NN.qpcr.seq)
save_plot("SplintSeqQ.png",splint.qpcr.seq)

# convert sequence names --------------------------------------------------

# This is just some one-off code to convert the sequences in the layout file to fasta format try and map them
# with the star pipeline to see if I can match up the names.

layout1 <- read_excel("layout-plate1.xlsx")
layout2 <- read_excel("layout-plate2.xlsx")

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

full.layout<-bind_rows(layout1,layout2) %>% 
  filter(!(is.na(sequence))) %>%
  mutate(dna.sequence=str_replace_all(sequence,"U","T")) %>%
  select(miRname,dna.sequence) %>%
  dplyr::rename("name"="miRname","seq"="dna.sequence")

writeFasta(full.layout,"miRnome.fasta")

qpcr.reads<-read.delim("qpcr_seqs2ReadsPerGene.out.tab",head=F)
nrow(qpcr.reads)
head(qpcr.reads)

qpcr.reads.filtered<-filter(qpcr.reads,!(V2 ==0 & V3 == 0 & V4 == 0))
sum(qpcr.reads.filtered[-1:-4,2])

full.layout

gtf<-rtracklayer::import("~/Desktop/STAR pipeline/miRNAannotation.gtf")
gtf_df=as_data_frame(gtf) #%>%
  #filter(type=="gene") %>%
  #dplyr::select(gene_id,gene_name)
test<-filter(gtf_df,seqnames == "chr5",
       end > 168768211,
       start < 168768211)
test$transcript_name
filter(qpcr.reads.filtered,V1=="ENSG00000207739.1")
ENSG00000207630.1
ENSG00000207703.1
ENSG00000207603.1

ENST00000384898.1
ENST00000384970.1
ENST00000384871.1