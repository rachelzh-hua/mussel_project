setwd("~/Documents/NuzhdinLab/mussel_project/analysis")

filenames <- list.files(path='counts/', full.names = TRUE)
all_data <- lapply(filenames, function(x)read.csv(x, header=TRUE, sep="\t"))
names(all_data) <- gsub('^.*[/]([^.]+)[.].*[.].*[.]counts$','\\1', filenames)

data1 <- list()
data2 <- list()

for (i in 1:length(all_data)){
  if(ncol(all_data[[i]])==19){
    #ps <- all_data[[i]]
    data1[[i]] <- all_data[[i]]
  } else{
    #ps2 <-  all_data[[i]]
    data2[[i]] <-all_data[[i]]
  }
}

names(data1) <- names(all_data)[1:7]
names(data2) <- names(all_data)

data1 <- data1[!sapply(data1, is.null)]
data2 <- data2[!sapply(data2, is.null)]


data<- data1[[1]]
data <- data[data$Class == "heterozygous" ,]
data <- data[
  data[,6]+data[,7]>20 & data[,8]+data[,9]>20 &
    data[,10]+data[,11]>20 & data[,12]+data[,13]>20 & data[,14]+data[,15]>20, 
]

data[,6] <- ifelse(data[,7]==0, data[,6]==NA, data[,6])
data[,7] <- ifelse(data[,6]==0, data[,7]==NA, data[,7])
data[,8] <- ifelse(data[,9]==0, data[,8]==NA, data[,8])
data[,9] <- ifelse(data[,8]==0, data[,9]==NA, data[,9])
data[,10] <- ifelse(data[,11]==0, data[,10]==NA, data[,10])
data[,11] <- ifelse(data[,10]==0, data[,11]==NA, data[,11])
data[,12] <- ifelse(data[,13]==0, data[,12]==NA, data[,12])
data[,13] <- ifelse(data[,12]==0, data[,13]==NA, data[,13])
data[,14] <- ifelse(data[,15]==0, data[,14]==NA, data[,14])
data[,15] <- ifelse(data[,14]==0, data[,15]==NA, data[,15])
data[,6:15][data[,6:15]==0] <- NA

data$D0_R <- data[,6]/(data[,7]+data[,6])
data$D0_A <- 1-data$D0_R
data$D11_R <- data[,8]/(data[,9]+data[,8])
data$D11_A <- 1-data$D11_R
data$D16_R <- data[,10]/(data[,11]+data[,10])
data$D16_A <- 1-data$D16_R
data$D23B_R <- data[,12]/(data[,13]+data[,12])
data$D23B_A <- 1- data$D23B_R
data$D23SM_R <- data[,14]/(data[,15]+data[,14])
data$D23SM_A <- 1 -data$D23SM_R

library(tidyverse)
library(dplyr)
library(tidyselect)

contig <- data[,c(1:4, 20:29)]
contig$FEMALE <- ifelse(data[,16]<10, 'ALT', 'REF')
contig$MALE <- ifelse(data[,18]<10, 'ALT', 'REF')
contig$D0_F <- ifelse(contig$FEMALE=='REF', contig[,5], contig[,6])
contig$D0_M <- ifelse(contig$MALE=='REF', contig[,5], contig[,6])
contig$D11_F <- ifelse(contig$FEMALE=='REF', contig[,7], contig[,8])
contig$D11_M <- ifelse(contig$MALE=='REF', contig[,7], contig[,8])
contig$D16_F <- ifelse(contig$FEMALE=='REF', contig[,9], contig[,10])
contig$D16_M <- ifelse(contig$MALE=='REF', contig[,9], contig[,10])
contig$D23B_F <- ifelse(contig$FEMALE=='REF', contig[,11], contig[,12])
contig$D23B_M <- ifelse(contig$MALE=='REF', contig[,11], contig[,12])
contig$D23SM_F <- ifelse(contig$FEMALE=='REF', contig[,13], contig[,14])
contig$D23SM_M <- ifelse(contig$MALE=='REF', contig[,13], contig[,14])
contig <- contig[,c(1:4, 17:26)]
contig <- contig %>% 
  group_by(X.CHROM)%>%
  summarise(D0_F=mean(D0_F),
            D0_M=mean(D0_M),
            D11_F=mean(D11_F),
            D11_M=mean(D11_M),
            D16_F=mean(D16_F),
            D16_M=mean(D16_M),
            D23B_F=mean(D23B_F),
            D23B_M=mean(D23B_M),
            D23SM_F=mean(D23SM_F),
            D23SM_M=mean(D23SM_M),)


