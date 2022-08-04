setwd("~/Documents/NuzhdinLab/mussel_project/analysis")

filenames <- list.files(path='contig_analysis/SNP/', full.names = TRUE)
all_data <- lapply(filenames, function(x)read.csv(x, header=TRUE, sep="\t"))
names(all_data) <- gsub('.*[//]([^_]+)([_])([^_]+)[_].*[_]chi.sq.txt$','\\1\\2\\3', filenames)

SNP_analysis <- data.frame(
  population = character(),
  number_of_SNP = integer(),
  number_of_sig_SNP = integer()
)

for(i in 1:length(all_data)){
  ROW_TO_ADD <- nrow(SNP_analysis) + 1
  
  SNP_analysis[ROW_TO_ADD, 'population'] <- names(all_data[i])
  SNP_analysis[ROW_TO_ADD, 'number_of_SNP'] <- nrow(all_data[[i]])
  SNP_analysis[ROW_TO_ADD, 'number_of_sig_SNP'] <- length(which(all_data[[i]]$q.value <= 0.05))
  
}
write.table(x=SNP_analysis, file='contig_analysis/SNP_analysis.txt', quote=FALSE, row.names=TRUE, sep='\t')

OUT <- createWorkbook()
for(i in 1:length(all_data)){
  addWorksheet(OUT, names(all_data[i]))
  writeData(OUT, names(all_data[i]), all_data[[i]])
}
saveWorkbook(OUT, 
             '~/Documents/NuzhdinLab/mussel_project/analysis/contig_analysis/all_SNP.xlsx',
             overwrite=T)

OUT1 <- createWorkbook()
for(i in 1:length(all_data)){
  addWorksheet(OUT1, names(all_data[i]))
  writeData(OUT1, names(all_data[i]), all_data[[i]][all_data[[i]]$q.value<=0.05, ])
}
saveWorkbook(OUT1, 
             'contig_analysis/significant_SNP.xlsx',
             overwrite=T)

data_1x1R1 <- merge(all_data[[2]], all_data[[1]], by='row.names', all.x =T)
names(data_1x1R1) <- gsub('([^.]+)[.]x$', '\\1_SM', names(data_1x1R1))
names(data_1x1R1) <- gsub('([^.]+)[.]y$', '\\1_B', names(data_1x1R1))
names(data_1x1R1) <- gsub('([^D0.reads$]+)[_].*', '\\1', names(data_1x1R1))
data_1x1R1 <- na.omit(data_1x1R1)

data_1x1R1 <- data_1x1R1[data_1x1R1$q.value_B<=0.05 & data_1x1R1$q.value_SM<=0.05, ]
data_1x1R1 <- data_1x1R1[,c(2:11, 18:21)]
OUT2 <- createWorkbook()
addWorksheet(OUT2, 'data_1x1R1_significant_SNPs')
writeData(OUT2, 
          1,
          data_1x1R1)
saveWorkbook(OUT2, 
             'contig_analysis/overlap_significant_SNPs.xlsx',
             overwrite=T)


data_1x1R2 <- all_data[[3]][all_data[[3]]$q.value<=0.05, ]
addWorksheet(OUT2, 'data_1x1R2_significant_SNPs')
writeData(OUT2, 2, data_1x1R2)
saveWorkbook(OUT2, 'contig_analysis/overlap_significant_SNPs.xlsx',overwrite=T)

data_1x1R3 <- all_data[[4]][all_data[[4]]$q.value<=0.05, ]
addWorksheet(OUT2, 'data_1x1R3_significant_SNPs')
writeData(OUT2, 3, data_1x1R3)
saveWorkbook(OUT2, 'contig_analysis/overlap_significant_SNPs.xlsx',overwrite=T)

data_1x3R1 <- merge(all_data[[5]], all_data[[6]], by='row.names', all.x =T)
names(data_1x3R1) <- gsub('([^.]+)[.]x$', '\\1_B', names(data_1x3R1))
names(data_1x3R1) <- gsub('([^.]+)[.]y$', '\\1_SM', names(data_1x3R1))
names(data_1x3R1) <- gsub('([^D0.reads$]+)[_].*', '\\1', names(data_1x3R1))
data_1x3R1 <- na.omit(data_1x3R1)
data_1x3R1 <- data_1x3R1[data_1x3R1$q.value_B<=0.05 & data_1x3R1$q.value_SM<=0.05, ]
data_1x3R1 <- data_1x3R1[,c(2:11, 18:21)]

addWorksheet(OUT2, 'data_1x3R1_significant_SNPs')
writeData(OUT2, 4, data_1x3R1)
saveWorkbook(OUT2, 'contig_analysis/overlap_significant_SNPs.xlsx',overwrite=T)

data_1x3R3 <- merge(all_data[[7]], all_data[[8]], by='row.names', all.x =T)
names(data_1x3R3) <- gsub('([^.]+)[.]x$', '\\1_B', names(data_1x3R3))
names(data_1x3R3) <- gsub('([^.]+)[.]y$', '\\1_SM', names(data_1x3R3))
names(data_1x3R3) <- gsub('([^D0.reads$]+)[_].*', '\\1', names(data_1x3R3))
data_1x3R3 <- na.omit(data_1x3R3)
data_1x3R3 <- data_1x3R3[data_1x3R3$q.value_B<=0.05 & data_1x3R3$q.value_SM<=0.05, ]
data_1x3R3 <- data_1x3R3[,c(2:11, 18:21)]

addWorksheet(OUT2, 'data_1x3R3_significant_SNPs')
writeData(OUT2, 5, data_1x3R3)
saveWorkbook(OUT2, 'contig_analysis/overlap_significant_SNPs.xlsx',overwrite=T)


data_2x1R1 <- merge(all_data[[9]], all_data[[10]], by='row.names', all.x =T)
names(data_2x1R1) <- gsub('([^.]+)[.]x$', '\\1_B', names(data_2x1R1))
names(data_2x1R1) <- gsub('([^.]+)[.]y$', '\\1_SM', names(data_2x1R1))
names(data_2x1R1) <- gsub('([^D0.reads$]+)[_].*', '\\1', names(data_2x1R1))
data_2x1R1 <- na.omit(data_2x1R1)
data_2x1R1 <- data_2x1R1[data_2x1R1$q.value_B<=0.05 & data_2x1R1$q.value_SM<=0.05, ]
data_2x1R1 <- data_2x1R1[,c(2:11, 18:21)]
addWorksheet(OUT2, 'data_2x1R1_significant_SNPs')
writeData(OUT2, 6, data_2x1R1)
saveWorkbook(OUT2, 'contig_analysis/overlap_significant_SNPs.xlsx',overwrite=T)

data_2x2R1 <- merge(all_data[[11]], all_data[[12]], by='row.names', all.x =T)
names(data_2x2R1) <- gsub('([^.]+)[.]x$', '\\1_B', names(data_2x2R1))
names(data_2x2R1) <- gsub('([^.]+)[.]y$', '\\1_SM', names(data_2x2R1))
names(data_2x2R1) <- gsub('([^D0.reads$]+)[_].*', '\\1', names(data_2x2R1))
data_2x2R1 <- na.omit(data_2x2R1)
data_2x2R1 <- data_2x2R1[data_2x2R1$q.value_B<=0.05 & data_2x2R1$q.value_SM<=0.05, ]
data_2x2R1 <- data_2x2R1[,c(2:11, 18:21)]
addWorksheet(OUT2, 'data_2x2R1_significant_SNPs')
writeData(OUT2, 7, data_2x2R1)
saveWorkbook(OUT2, 'contig_analysis/overlap_significant_SNPs.xlsx',overwrite=T)

data_2x3R2 <- all_data[[13]][all_data[[13]]$q.value<=0.05, ]
addWorksheet(OUT2, 'data_2x3R2_significant_SNPs')
writeData(OUT2, 8, data_2x3R2)
saveWorkbook(OUT2, 'contig_analysis/overlap_significant_SNPs.xlsx',overwrite=T)




