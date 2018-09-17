median_centering <- function(input.data,subsetN=100,subsetSize=100,percentERpos=0.5,
                             idVar="PATIENT_ID"){
  gene_names <- c("MDM2","EGFR","PGR","GRB7","RRM2","KRT5","BIRC5","CDCA1","PTTG1","BLVRA",  
                  "MKI67","NAT1","MYBL2","MAPT","CXXC5","FOXA1","ORC6L","MLPH","BCL2","MELK",   
                  "MYC","GPR160","CDC20","ANLN","UBE2T","CEP55","TMEM45B","TYMS","ERBB2","EXO1",   
                  "SFRP1","ESR1","FGFR4","MIA","MMP11","CCNB1","CCNE1","PHGDH","KNTC2","BAG1",   
                  "KRT14","UBE2C","CENPF","CDC6","KRT17","CDH3","ACTR3B","FOXC1","KIF2C","SLC39A6")
  med_exp <- data.frame(genes=gene_names)
  
  for(i in 1:subsetN){
    if(!is.na(percentERpos)){
      positive <- sample_n(input.data[input.data$`ER Status`=="Positive",], percentERpos*subsetSize, replace=TRUE)
      negative <- sample_n(input.data[input.data$`ER Status`=="Negative",], (1-percentERpos)*subsetSize, replace=TRUE)
      dat <- as.data.frame(rbind(positive,negative))
      med_exp[,i+1] <- apply(dat[,colnames(dat) %in% gene_names],2,median)
    } else{
      dat <- sample_n(input.data, subsetSize, replace=TRUE)
      med_exp[,i+1] <- apply(dat[,colnames(dat) %in% gene_names],2,median)
    }
  }
  
  avg_med_exp <- apply(med_exp[,2:ncol(med_exp)],1,mean)
  
  for_pam50 <- as.data.frame(t(input.data[,colnames(input.data) %in% gene_names]))
  colnames(for_pam50) <- input.data[,c(idVar)]
  
  for(i in 1:nrow(for_pam50)){
    for_pam50[i,] <- apply(for_pam50[i,],1,function(x) as.numeric(x)-avg_med_exp[[i]])
  }
  
  for_pam50 <- as.data.frame(cbind(gene=row.names(for_pam50),
                                   for_pam50))
  return(for_pam50)
  
}

