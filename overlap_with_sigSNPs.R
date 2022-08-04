setwd("~/Documents/NuzhdinLab/mussel_project/analysis")

library(readxl)
sheet_names <- excel_sheets('size_QTL/size_significantSNPs.xlsx')
size_SNPs <- lapply(sheet_names,
                   function(x) read_excel('size_QTL/size_significantSNPs.xlsx', sheet=x))
names(size_SNPs) <- sheet_names

sheet_names_2 <- excel_sheets('chisq.test/2.0_chi.sq_significant_markers.xlsx')
chisq_vSNPs <- lapply(sheet_names_2,
                    function(x) read_excel('chisq.test/2.0_chi.sq_significant_markers.xlsx', sheet=x))
names(chisq_vSNPs ) <- sheet_names_2
chisq_vSNPs <- chisq_vSNPs[-c(2,3,8)]

common <- intersect(chisq_vSNPs[[3]]$POS, size_SNPs[[3]]$POS)

data_1x1R1 <- merge(chisq_vSNPs[[1]], size_SNPs[[1]], by='POS')
data_1x3R1 <- merge(chisq_vSNPs[[2]], size_SNPs[[2]], by='POS')
data_1x3R3 <- merge(chisq_vSNPs[[3]], size_SNPs[[3]], by='POS')
data_2x1R1 <- merge(chisq_vSNPs[[4]], size_SNPs[[4]], by='POS')
data_2x2R1 <- merge(chisq_vSNPs[[5]], size_SNPs[[5]], by='POS')

overlap_1 <- lapply(ls(pattern='data[_]+'), function(x) get(x))
names(overlap_1) <- ls(pattern='data[_]+')

overlap_chi.sq <- data.frame(
  replicate = character(),
  number_of_size_SNP = integer(),
  number_of_sig_SNP_chi.sq = integer(),
  overlap = integer()
)

for(i in 1:length(overlap_1)){
  ROW_TO_ADD <- nrow(overlap_chi.sq) + 1
  
  overlap_chi.sq[ROW_TO_ADD, 'replicate'] <- names(overlap_1[i])
  overlap_chi.sq[ROW_TO_ADD, 'overlap'] <- nrow(overlap_1[[i]])
}

number_of_sig_SNP <- c()
for(j in 1:length(chisq_vSNPs)){
  number_of_sig_SNP <- c(number_of_sig_SNP, nrow(chisq_vSNPs[[j]]))
}
overlap_chi.sq$number_of_sig_SNP_chi.sq <- number_of_sig_SNP

number_of_size_SNP <- c()
for(h in 1:length(size_SNPs)){
  number_of_size_SNP <- c(number_of_size_SNP, nrow(size_SNPs[[h]]))
}
overlap_chi.sq$number_of_size_SNP <- number_of_size_SNP

write.table(x=overlap_chi.sq, file='size_QTL/size_overlap_chi.sq.txt', quote=FALSE, row.names=TRUE, sep='\t')

library(openxlsx)
OUT <- createWorkbook()
for(i in 1:length(overlap_1)){
  addWorksheet(OUT, names(overlap_1[i]))
  writeData(OUT, names(overlap_1[i]), overlap_1[[i]])
}
saveWorkbook(OUT, 
             '~/Documents/NuzhdinLab/mussel_project/analysis/size_QTL/size_overlap_chi.sq.xlsx',
             overwrite=T)




sheet_names_3 <- excel_sheets('select_marker/withContig/significantSNPs.xlsx')
uniform_vSNPs <- lapply(sheet_names_3,
                      function(x) read_excel('select_marker/withContig/significantSNPs.xlsx', sheet=x))
names(uniform_vSNPs ) <- sheet_names_3

common <- intersect(chisq_vSNPs[[3]]$POS, size_SNPs[[3]]$POS)

data_1x1R1 <- merge(rbind(uniform_vSNPs[[1]],uniform_vSNPs[[2]]), size_SNPs[[1]], by='POS')
data_1x1R1 <- data_1x1R1[data_1x1R1$X.CHROM.x==data_1x1R1$X.CHROM.y, ]

data_1x3R1 <- merge(rbind(uniform_vSNPs[[7]],uniform_vSNPs[[8]]), size_SNPs[[2]], by='POS')

data_1x3R3 <- merge(rbind(uniform_vSNPs[[9]],uniform_vSNPs[[10]]), size_SNPs[[3]], by='POS')
data_1x3R3 <- data_1x3R3[data_1x3R3$X.CHROM.x==data_1x3R3$X.CHROM.y, ]

data_2x1R1 <- merge(rbind(uniform_vSNPs[[11]],uniform_vSNPs[[12]]), size_SNPs[[4]], by='POS')

data_2x2R1 <- merge(rbind(uniform_vSNPs[[13]],uniform_vSNPs[[14]]), size_SNPs[[5]], by='POS')

overlap_1 <- lapply(ls(pattern='data[_]+'), function(x) get(x))
names(overlap_1) <- ls(pattern='data[_]+')

overlap_chi.sq <- data.frame(
  replicate = character(),
  number_of_size_SNP = integer(),
  number_of_sig_SNP_chi.sq = integer(),
  overlap = integer()
)

for(i in 1:length(overlap_1)){
  ROW_TO_ADD <- nrow(overlap_chi.sq) + 1
  
  overlap_chi.sq[ROW_TO_ADD, 'replicate'] <- names(overlap_1[i])
  overlap_chi.sq[ROW_TO_ADD, 'overlap'] <- nrow(overlap_1[[i]])
}

number_of_sig_SNP <- c()
for(j in 1:length(chisq_vSNPs)){
  number_of_sig_SNP <- c(number_of_sig_SNP, nrow(chisq_vSNPs[[j]]))
}
overlap_chi.sq$number_of_sig_SNP_chi.sq <- number_of_sig_SNP

number_of_size_SNP <- c()
for(h in 1:length(size_SNPs)){
  number_of_size_SNP <- c(number_of_size_SNP, nrow(size_SNPs[[h]]))
}
overlap_chi.sq$number_of_size_SNP <- number_of_size_SNP

write.table(x=overlap_chi.sq, file='size_QTL/size_overlap_chi.sq.txt', quote=FALSE, row.names=TRUE, sep='\t')

library(openxlsx)
OUT <- createWorkbook()
for(i in 1:length(overlap_1)){
  addWorksheet(OUT, names(overlap_1[i]))
  writeData(OUT, names(overlap_1[i]), overlap_1[[i]])
}
saveWorkbook(OUT, 
             '~/Documents/NuzhdinLab/mussel_project/analysis/size_QTL/size_overlap_chi.sq.xlsx',
             overwrite=T)

