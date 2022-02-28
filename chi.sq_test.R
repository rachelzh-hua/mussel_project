setwd("~/Documents/NuzhdinLab/mussel_project/analysis")

filenames <- list.files(path='counts/', full.names = TRUE)
results <- strsplit(x=filenames, split='//', fixed=TRUE)

population <- sapply(
  X=results,
  FUN=function(x){x[2]}
)

population <- sapply(
  strsplit(x=population, split='.', fixed=T),
  FUN=function(x){x[1]}
)

chi.sq_table <- data.frame(
  replicate=character(),
  number_of_markers=integer(),
  total_number_het_markers=integer(),
  D0_D23B=integer(),
  D0_D23SM=integer()
)
#ROW_TO_ADD <- nrow(chi.sq_table) +1

data <- read.csv(file=filenames[7], header=TRUE, sep="\t")
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

data_23B <-data[,c(1,2,3,4,6,7,12,13)]
data_23B <-na.omit(data_23B)
data_23SM <- data[,c(1,2,3,4,6,7,14,15)]
data_23SM <- na.omit(data_23SM)

SNP_list <- list()

for(i in 1:nrow(data_23B)){
  SNP_list[[paste("SNP", i, sep = "_")]] <- matrix(nrow=2, ncol=2, 
                                                   dimnames = list(c("ref", "alt"), c("D0", "D23_B")))
  SNP_list[[paste("SNP", i, sep = "_")]][1,1] <- data_23B[i,5]
  SNP_list[[paste("SNP", i, sep = "_")]][1,2] <- data_23B[i,7]
  SNP_list[[paste("SNP", i, sep = "_")]][2,1] <- data_23B[i,6]
  SNP_list[[paste("SNP", i, sep = "_")]][2,2] <- data_23B[i,8]
  SNP_list[[paste("SNP", i, sep = "_")]]<-as.data.frame(SNP_list[[paste("SNP", i, sep = "_")]])
  
}

names(SNP_list) <- rownames(data_23B)

pvals <- c()
for(j in 1:length(SNP_list)){
  chisq.result <- chisq.test(SNP_list[[j]])
  pvals <- c(pvals,chisq.result$p.value)
}

#pvals_distribution <- hist(pvals, breaks=20)

data_23B_chisq <- data.frame(
  SNP_rownames = integer(length = nrow(data_23B)),
  X.CHROM = character(length = nrow(data_23B)),
  POS=integer(length = nrow(data_23B)),
  p.values = integer(length = nrow(data_23B))
)

data_23B_chisq$SNP_rownames <- rownames(data_23B)
data_23B_chisq$X.CHROM <- data_23B$X.CHROM
data_23B_chisq$POS <- data_23B$POS
data_23B_chisq$p.values <- pvals

#data_23B_chisq$p.adjust.BH <- p.adjust(p=data_23B_chisq$p.values, method='BH')
#length(which(data_23B_chisq$p.adjust.BH <= 0.05))
#data_23B_chisq$p.adjust.bonferroni <- p.adjust(p=data_23B_chisq$p.values, method='bonferroni')
#length(which(data_23B_chisq$p.adjust.bonferroni <= 0.05))

#data_23B_bonferroni <- data_23B_chisq[data_23B_chisq$p.adjust.bonferroni <= 0.05, ] 
#data_23B_bonferroni <- na.omit(data_23B_bonferroni) #634

SNP_list_SM <- list()

for(i in 1:nrow(data_23SM)){
  SNP_list_SM[[paste("SNP", i, sep = "_")]] <- matrix(nrow=2, ncol=2, 
                                                      dimnames = list(c("ref", "alt"), c("D0", "D23_B")))
  SNP_list_SM[[paste("SNP", i, sep = "_")]][1,1] <- data_23SM[i,5]
  SNP_list_SM[[paste("SNP", i, sep = "_")]][1,2] <- data_23SM[i,7]
  SNP_list_SM[[paste("SNP", i, sep = "_")]][2,1] <- data_23SM[i,6]
  SNP_list_SM[[paste("SNP", i, sep = "_")]][2,2] <- data_23SM[i,8]
  SNP_list_SM[[paste("SNP", i, sep = "_")]]<-as.data.frame(SNP_list_SM[[paste("SNP", i, sep = "_")]])
  
}

names(SNP_list_SM) <- rownames(data_23SM)

pvals_SM <- c()
for(j in 1:length(SNP_list_SM)){
  chisq.result_SM <- chisq.test(SNP_list_SM[[j]])
  pvals_SM <- c(pvals_SM,chisq.result_SM$p.value)
}

#pvals_distribution_SM <- hist(pvals, breaks=20)

data_23SM_chisq <- data.frame(
  SNP_rownames = integer(length = nrow(data_23SM)),
  X.CHROM = character(length = nrow(data_23SM)),
  POS = integer(length = nrow(data_23SM)),
  p.values = integer(length = nrow(data_23SM))
)

data_23SM_chisq$SNP_rownames <- rownames(data_23SM)
data_23SM_chisq$X.CHROM <- data_23SM$X.CHROM
data_23SM_chisq$POS <- data_23SM$POS
data_23SM_chisq$p.values <- pvals_SM
#length(which(data_23SM_chisq$p.values <= 0.05)) #1900

#data_23SM_chisq$p.adjust.BH <- p.adjust(p=data_23SM_chisq$p.values, method='BH')
#length(which(data_23SM_chisq$p.adjust.BH <= 0.05))  
#data_23SM_chisq$p.adjust.bonferroni <- p.adjust(p=data_23SM_chisq$p.values, method='bonferroni')
#length(which(data_23SM_chisq$p.adjust.bonferroni <= 0.05))   

#data_23SM_bonferroni <- data_23SM_chisq[data_23SM_chisq$p.adjust.bonferroni <= 0.05, ]
#data_23SM_bonferroni <- na.omit(data_23SM_bonferroni) #395

#library(dplyr)
#common_bonferroni <- intersect(rownames(data_23B_bonferroni), rownames(data_23SM_bonferroni))  
#common_bonferroni_table <- data_23SM_bonferroni[common_bonferroni,]  
#common_bonferroni_table <- na.omit(common_bonferroni_table)
#length(which(common_bonferroni_table$p.adjust.bonferroni <= 0.05))  

data_23B_pval <- data_23B_chisq[data_23B_chisq$p.values < 0.05/nrow(data), ]
data_23SM_pval <- data_23SM_chisq[data_23SM_chisq$p.values < 0.05/nrow(data), ]
common <- intersect(data_23B_pval$SNP_rownames, data_23SM_pval$SNP_rownames) 
common_pval <- data[common,]
p.value_B <- data_23B_pval[data_23B_pval$SNP_rownames %in% common, ]
common_pval$p.value_B <- p.value_B$p.values
p.value_SM <- data_23SM_pval[data_23SM_pval$SNP_rownames %in% common, ]
common_pval$p.value_SM <- p.value_SM$p.values

#chi.sq_table[ROW_TO_ADD+3, 'replicate'] <- population[1]
#chi.sq_table[ROW_TO_ADD+3, 'number_of_markers'] <- length(common_pval)
#chi.sq_table[ROW_TO_ADD+3, 'total_number_het_markers'] <- nrow(data)
#chi.sq_table[ROW_TO_ADD+3, 'D0_D23B'] <- length(which(data_23B_chisq$p.values < 0.05))
#chi.sq_table[ROW_TO_ADD+3, 'D0_D23SM'] <- length(which(data_23SM_chisq$p.values < 0.05))
chi.sq_table[nrow(chi.sq_table)+1,] <- c(population[5], 
                                        nrow(common_pval),
                                        nrow(data), 
 
                                    nrow(data_23B_pval),
                                    nrow(data_23SM_pval))


#library(openxlsx)
#OUT <- createWorkbook()
addWorksheet(OUT, paste(population[7],'significant_markers', sep='_'))
writeData(OUT, 
          paste(population[7],'significant_markers', sep='_'),
          common_pval)
saveWorkbook(OUT, 
             '~/Documents/NuzhdinLab/mussel_project/analysis/chisq.test/2.0_chi.sq_significant_markers.xlsx',
             overwrite=T)



data <- read.csv(file=filenames[8], header=TRUE, sep="\t")
data <- data[data$Class == "heterozygous" ,]
data <- data[
  data[,6]+data[,7]>20 & data[,8]+data[,9]>20 &
    data[,10]+data[,11]>20 & data[,12]+data[,13]>20, 
]

data[,6] <- ifelse(data[,7]==0, data[,6]==NA, data[,6])
data[,7] <- ifelse(data[,6]==0, data[,7]==NA, data[,7])
data[,8] <- ifelse(data[,9]==0, data[,8]==NA, data[,8])
data[,9] <- ifelse(data[,8]==0, data[,9]==NA, data[,9])
data[,10] <- ifelse(data[,11]==0, data[,10]==NA, data[,10])
data[,11] <- ifelse(data[,10]==0, data[,11]==NA, data[,11])
data[,12] <- ifelse(data[,13]==0, data[,12]==NA, data[,12])
data[,13] <- ifelse(data[,12]==0, data[,13]==NA, data[,13])
data[,6:13][data[,6:13]==0] <- NA

data_23B <-data[,c(1,2,3,4,6,7,12,13)]
data_23B <- na.omit(data_23B)

SNP_list <- list()

for(i in 1:nrow(data_23B)){
  SNP_list[[paste("SNP", i, sep = "_")]] <- matrix(nrow=2, ncol=2, 
                                                   dimnames = list(c("ref", "alt"), c("D0", "D23_B")))
  SNP_list[[paste("SNP", i, sep = "_")]][1,1] <- data_23B[i,5]
  SNP_list[[paste("SNP", i, sep = "_")]][1,2] <- data_23B[i,7]
  SNP_list[[paste("SNP", i, sep = "_")]][2,1] <- data_23B[i,6]
  SNP_list[[paste("SNP", i, sep = "_")]][2,2] <- data_23B[i,8]
  SNP_list[[paste("SNP", i, sep = "_")]]<-as.data.frame(SNP_list[[paste("SNP", i, sep = "_")]])
  
}

names(SNP_list) <- rownames(data_23B)

pvals <- c()
for(j in 1:length(SNP_list)){
  chisq.result <- chisq.test(SNP_list[[j]])
  pvals <- c(pvals,chisq.result$p.value)
}

#pvals_distribution <- hist(pvals, breaks=20)

data_23B_chisq <- data.frame(
  SNP_rownames = integer(length = nrow(data_23B)),
  X.CHROM = character(length = nrow(data_23B)),
  POS = integer(length = nrow(data_23B)),
  p.values = integer(length = nrow(data_23B))
)

data_23B_chisq$SNP_rownames <- rownames(data_23B)
data_23B_chisq$X.CHROM <- data_23B$X.CHROM
data_23B_chisq$POS <- data_23B$POS
data_23B_chisq$p.values <- pvals

#data_23B_chisq$p.adjust.BH <- p.adjust(p=data_23B_chisq$p.values, method='BH')
#length(which(data_23B_chisq$p.adjust.BH <= 0.05))  #14390
#data_23B_chisq$p.adjust.bonferroni <- p.adjust(p=data_23B_chisq$p.values, method='bonferroni')
#length(which(data_23B_chisq$p.adjust.bonferroni <= 0.05))  #7521 

#data_23B_bonferroni <- data_23B_chisq[data_23B_chisq$p.adjust.bonferroni <= 0.05, ] 
#data_23B_bonferroni <- na.omit(data_23B_bonferroni) #4342

data_23B_pval <- data_23B_chisq[data_23B_chisq$p.values < 0.05/nrow(data_23B), ]
chi.sq_table[nrow(chi.sq_table)+1,] <- c(population[8], 
                                         nrow(data_23B_pval),
                                         nrow(data), 
                                         nrow(data_23B_pval),
                                         NA
                                         )

common_pval <- data[data_23B_pval$SNP_rownames, ]
common_pval$p.value_B <- data_23B_pval$p.values


addWorksheet(OUT, paste(population[8],'significant_markers', sep='_'))
writeData(OUT, 
          paste(population[8],'significant_markers', sep='_'),
          common_pval)
saveWorkbook(OUT, 
             '~/Documents/NuzhdinLab/mussel_project/analysis/chisq.test/2.0_chi.sq_significant_markers.xlsx',
             overwrite=T)





write.table(x=chi.sq_table, 
            file='~/Documents/NuzhdinLab/mussel_project/analysis/chisq.test/chi.sq_table_p.values_over_total.txt', 
            quote=FALSE, row.names=TRUE, sep='\t')





