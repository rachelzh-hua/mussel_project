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

size <- data[,c(1:4, 12:15)]
size <- na.omit(size)
SNP_list <- list()
for(i in 1:nrow(size)){
  SNP_list[[paste("SNP", i, sep = "_")]] <- matrix(nrow=2, ncol=2, 
                                                   dimnames = list(c("ref", "alt"), c("D23B", "D23SM")))
  SNP_list[[paste("SNP", i, sep = "_")]][1,1] <- size[i,5]
  SNP_list[[paste("SNP", i, sep = "_")]][1,2] <- size[i,7]
  SNP_list[[paste("SNP", i, sep = "_")]][2,1] <- size[i,6]
  SNP_list[[paste("SNP", i, sep = "_")]][2,2] <- size[i,8]
  SNP_list[[paste("SNP", i, sep = "_")]]<-as.data.frame(SNP_list[[paste("SNP", i, sep = "_")]])
  
}
names(SNP_list) <- rownames(size)
pvals <- c()
for(j in 1:length(SNP_list)){
  chisq.result <- chisq.test(SNP_list[[j]])
  pvals <- c(pvals,chisq.result$p.value)
}

size$pvals <- pvals 


library(openxlsx)
OUT <- createWorkbook()
addWorksheet(OUT, paste(names(data1[1]),'all_size', sep='_'))
writeData(OUT, 
          paste(names(data1[1]),'all_size', sep='_'),
          size)
saveWorkbook(OUT, 
             '~/Documents/NuzhdinLab/mussel_project/analysis/size_QTL/all_size.xlsx',
             overwrite=T)



data<- data1[[5]]
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

size <- data[,c(1:4, 12:15)]
size <- na.omit(size)
SNP_list <- list()
for(i in 1:nrow(size)){
  SNP_list[[paste("SNP", i, sep = "_")]] <- matrix(nrow=2, ncol=2, 
                                                   dimnames = list(c("ref", "alt"), c("D23B", "D23SM")))
  SNP_list[[paste("SNP", i, sep = "_")]][1,1] <- size[i,5]
  SNP_list[[paste("SNP", i, sep = "_")]][1,2] <- size[i,7]
  SNP_list[[paste("SNP", i, sep = "_")]][2,1] <- size[i,6]
  SNP_list[[paste("SNP", i, sep = "_")]][2,2] <- size[i,8]
  SNP_list[[paste("SNP", i, sep = "_")]]<-as.data.frame(SNP_list[[paste("SNP", i, sep = "_")]])
  
}
names(SNP_list) <- rownames(size)
pvals <- c()
for(j in 1:length(SNP_list)){
  chisq.result <- chisq.test(SNP_list[[j]])
  pvals <- c(pvals,chisq.result$p.value)
}

size$pvals <- pvals 


addWorksheet(OUT, paste(names(data1[5]),'all_size', sep='_'))
writeData(OUT, 
          paste(names(data1[5]),'all_size', sep='_'),
          size)
saveWorkbook(OUT, 
             '~/Documents/NuzhdinLab/mussel_project/analysis/size_QTL/all_size.xlsx',
             overwrite=T)



