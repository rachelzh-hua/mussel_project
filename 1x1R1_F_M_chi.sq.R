setwd("~/Documents/NuzhdinLab/mussel_project/analysis")

filenames <- list.files(path='counts/', full.names = TRUE)
all_data <- lapply(filenames, function(x)read.csv(x, header=TRUE, sep="\t"))
names(all_data) <- gsub('^.*[/]([^.]+)[.].*[.].*[.]counts$','\\1', filenames)
data <- all_data[[1]]
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
data$FEMALE <- ifelse(data[,16]<10, 'ALT', 'REF')
data$MALE <- ifelse(data[,18]<10, 'ALT', 'REF')
data$D0_num_F_reads <- ifelse(data[,30]=='REF', data[,20]*(data[,7]+data[,6]), data[,21]*(data[,7]+data[,6]))
data$D0_num_M_reads <- ifelse(data[,31]=='REF', data[,20]*(data[,7]+data[,6]), data[,21]*(data[,7]+data[,6]))
data$D23B_num_F_reads <- ifelse(data[,30]=='REF', data[,26]*(data[,13]+data[,12]), data[,27]*(data[,13]+data[,12]))
data$D23B_num_M_reads <- ifelse(data[,31]=='REF', data[,26]*(data[,13]+data[,12]), data[,27]*(data[,13]+data[,12]))
data$D23SM_num_F_reads <- ifelse(data[,30]=='REF', data[,28]*(data[,15]+data[,14]), data[,29]*(data[,15]+data[,14]))
data$D23SM_num_M_reads <- ifelse(data[,31]=='REF', data[,28]*(data[,15]+data[,14]), data[,29]*(data[,15]+data[,14]))

data_23B <-data[,c(1,2,3,4,32:35)]
data_23SM <- data[,c(1,2,3,4,32,33,36,37)]

library(tidyverse)
library(dplyr)
library(tidyselect)

data_23B <- data_23B %>% group_by(X.CHROM) %>% 
  summarise(D0_num_F_reads=sum(D0_num_F_reads),
            D0_num_M_reads=sum(D0_num_M_reads),
            D23B_num_F_reads=sum(D23B_num_F_reads),
            D23B_num_M_reads=sum(D23B_num_M_reads))
data_23B <-na.omit(data_23B)
data_23B <-as.data.frame(data_23B)
data_23SM <- data_23SM %>% group_by(X.CHROM) %>% 
  summarise(D0_num_F_reads=sum(D0_num_F_reads),
            D0_num_M_reads=sum(D0_num_M_reads),
            D23SM_num_F_reads=sum(D23SM_num_F_reads),
            D23SM_num_M_reads=sum(D23SM_num_M_reads))
data_23SM <- na.omit(data_23SM)
data_23SM <- as.data.frame(data_23SM)

SNP_list <- list()

for(i in 1:nrow(data_23B)){
  SNP_list[[paste("SNP", i, sep = "_")]] <- matrix(nrow=2, ncol=2, 
                                                   dimnames = list(c('F', "M"), c("D0", "D23_B")))
  SNP_list[[paste("SNP", i, sep = "_")]][1,1] <- data_23B[i,2]
  SNP_list[[paste("SNP", i, sep = "_")]][1,2] <- data_23B[i,4]
  SNP_list[[paste("SNP", i, sep = "_")]][2,1] <- data_23B[i,3]
  SNP_list[[paste("SNP", i, sep = "_")]][2,2] <- data_23B[i,5]
  SNP_list[[paste("SNP", i, sep = "_")]]<-as.data.frame(SNP_list[[paste("SNP", i, sep = "_")]])
  
}

names(SNP_list) <- data_23B$X.CHROM

pvals <- c()
for(j in 1:length(SNP_list)){
  chisq.result <- chisq.test(SNP_list[[j]])
  pvals <- c(pvals,chisq.result$p.value)
}
data_23B$p.value <- pvals
library(qvalue)
res <- qvalue(data_23B$p.value)
data_23B$q.value <- res$qvalues
write.table(x=data_23B, file='select_marker/1x1R1_23B_contig_F_M_chi.sq.txt', quote=FALSE, row.names=TRUE, sep='\t')



SNP_list <- list()

for(i in 1:nrow(data_23SM)){
  SNP_list[[paste("SNP", i, sep = "_")]] <- matrix(nrow=2, ncol=2, 
                                                   dimnames = list(c('F', "M"), c("D0", "D23_SM")))
  SNP_list[[paste("SNP", i, sep = "_")]][1,1] <- data_23SM[i,2]
  SNP_list[[paste("SNP", i, sep = "_")]][1,2] <- data_23SM[i,4]
  SNP_list[[paste("SNP", i, sep = "_")]][2,1] <- data_23SM[i,3]
  SNP_list[[paste("SNP", i, sep = "_")]][2,2] <- data_23SM[i,5]
  SNP_list[[paste("SNP", i, sep = "_")]]<-as.data.frame(SNP_list[[paste("SNP", i, sep = "_")]])
  
}

names(SNP_list) <- data_23SM$X.CHROM

pvals <- c()
for(j in 1:length(SNP_list)){
  chisq.result <- chisq.test(SNP_list[[j]])
  pvals <- c(pvals,chisq.result$p.value)
}
data_23SM$p.value <- pvals
res_SM <- qvalue(data_23SM$p.value)
data_23SM$q.value <- res_SM$qvalues
write.table(x=data_23SM, file='select_marker/1x1R1_23SM_contig_F_M_chi.sq.txt', quote=FALSE, row.names=TRUE, sep='\t')





















