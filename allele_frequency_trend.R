setwd("~/Documents/NuzhdinLab/mussel_project/analysis")

filenames <- list.files(path='counts/', full.names = TRUE)
all_data <- lapply(filenames, function(x)read.csv(x, header=TRUE, sep="\t"))
names(all_data) <- gsub('^.*[/]([^.]+)[.].*[.].*[.]counts$','\\1', filenames)

for (i in 1:length(all_data)){
  all_data[[i]] <- all_data[[i]][all_data[[i]][,5]=='heterozygous',]
  
  if(ncol(all_data[[i]])==19){
    all_data[[i]] <- all_data[[i]][all_data[[i]][,6]+all_data[[i]][,7] >20 &
                                     all_data[[i]][,8]+all_data[[i]][,9] >20 &
                                     all_data[[i]][,10]+all_data[[i]][,11] >20 &
                                     all_data[[i]][,12]+all_data[[i]][,13] >20 &
                                     all_data[[i]][,14]+all_data[[i]][,15] >20, ]
    
    all_data[[i]]$D0 <- all_data[[i]][,7]/(all_data[[i]][,7] + all_data[[i]][,6])
    all_data[[i]]$D11 <- all_data[[i]][,9]/(all_data[[i]][,9] + all_data[[i]][,8])-all_data[[i]]$D0 
    all_data[[i]]$D16 <- all_data[[i]][,11]/(all_data[[i]][,11] + all_data[[i]][,10])-all_data[[i]]$D0 
    all_data[[i]]$D23B <- all_data[[i]][,13]/(all_data[[i]][,13] + all_data[[i]][,12])-all_data[[i]]$D0
    all_data[[i]]$D23SM <- all_data[[i]][,15]/(all_data[[i]][,15] + all_data[[i]][,14])-all_data[[i]]$D0
    
    all_data[[i]] <- all_data[[i]][,20:24]
  } else{
    all_data[[i]] <- all_data[[i]][all_data[[i]][,6]+all_data[[i]][,7] >20 &
                                     all_data[[i]][,8]+all_data[[i]][,9] >20 &
                                     all_data[[i]][,10]+all_data[[i]][,11] >20 &
                                     all_data[[i]][,12]+all_data[[i]][,13] >20,] 
    
    all_data[[i]]$D0 <- all_data[[i]][,7]/(all_data[[i]][,7] + all_data[[i]][,6])
    all_data[[i]]$D11 <- all_data[[i]][,9]/(all_data[[i]][,9] + all_data[[i]][,8])-all_data[[i]]$D0 
    all_data[[i]]$D16 <- all_data[[i]][,11]/(all_data[[i]][,11] + all_data[[i]][,10])-all_data[[i]]$D0 
    all_data[[i]]$D23B <- all_data[[i]][,13]/(all_data[[i]][,13] + all_data[[i]][,12])-all_data[[i]]$D0
    
    all_data[[i]] <- all_data[[i]][,18:21]
  }
}

D11 <- c()
D16 <- c()
D23B <- c()
D23SM <- c()

for(j in 1:length(all_data)){
  D11 <- c(D11, mean(all_data[[j]][,2],na.rm=T))
  D16 <- c(D16, mean(all_data[[j]][,3], na.rm =T))
  D23B <- c(D23B, mean(all_data[[j]][,4], na.rm =T))
  if(ncol(all_data[[j]])==5){
    D23SM <- c(D23SM, mean(all_data[[j]][,5], na.rm =T))
  }else{
    D23SM <- c(D23SM, NA)
  }
}

seq_divergence <- data.frame(
  D0 <- 0,
  D11 <- D11,
  D16 <- D16,
  D23B <- D23B,
  D23SM <- D23SM
)

rownames(seq_divergence) <- names(all_data)
colnames(seq_divergence) <- c('D0', 'D11', 'D16', 'D23B', 'D23SM')

write.table(
  seq_divergence,
  file='~/Documents/NuzhdinLab/mussel_project/analysis/heterozygosity/seq_divergence_table.txt',
  quote=FALSE, row.names=TRUE, sep='\t'
)