setwd("~/Documents/NuzhdinLab/mussel_project/analysis")
library(readxl)
sheet_names <- excel_sheets('size_QTL/all_size.xlsx')
sig_size <- lapply(sheet_names,
                      function(x) read_excel('size_QTL/all_size.xlsx', sheet=x))
names(sig_size) <- sheet_names

library(openxlsx)
OUT <- createWorkbook()

for(i in 1:length(sig_size)){
  size <- sig_size[[i]][sig_size[[i]]$pvals < 0.05/nrow(sig_size[[i]]), ]
  addWorksheet(OUT, paste(names(sig_size)[i],'significantSNPs', sep='_'))
  writeData(OUT, paste(names(sig_size)[i],'significantSNPs', sep='_'), size)
  saveWorkbook(OUT, '~/Documents/NuzhdinLab/mussel_project/analysis/size_QTL/size_significantSNPs.xlsx', overwrite=T)
  
}

sizeQTL_table <- data.frame(
  replicate=character(length = 5L),
  total_number_het_markers=integer(length = 5L),
  number_sizeQTL=integer(length = 5L)
)

sheet_names2 <- excel_sheets('size_QTL/size_significantSNPs.xlsx')
sizeQTL <- lapply(sheet_names2,
                   function(x) read_excel('size_QTL/size_significantSNPs.xlsx', sheet=x))
names(sizeQTL) <- sheet_names2

sizeQTL_table$replicate <- gsub('^(.*)[_].*[_]size$','\\1', sheet_names)

total_number <- c()
for(j in 1:length(sig_size)){
  total_number <- c(total_number, nrow(sig_size[[j]]))
}
sizeQTL_table$total_number_het_markers <- total_number

sig_number <- c()
for(h in 1:length(sizeQTL)){
  sig_number <- c(sig_number, nrow(sizeQTL[[h]]))
}

sizeQTL_table$number_sizeQTL <- sig_number

write.table(
  sizeQTL_table,
  file = '~/Documents/NuzhdinLab/mussel_project/analysis/size_QTL/sizeQTL_table.txt',
  quote=FALSE, row.names=TRUE, sep='\t'
)

sizeQTL_table$percentage <- sizeQTL_table[,3]/sizeQTL_table[,2]
