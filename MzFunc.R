#Libraries
{
  library(Peptides) 
  library(stringr)
}

MzFunc <- function(mainsequence,carbaminomethylation=F,methoxid=F,deamidation=F,acetylation=F,propylation=T,methylAll=T,C13D3=F,CD3=F,Addition=NA){
  if((C13D3==T & CD3==F & methylAll==F) | (C13D3==F & CD3==T & methylAll==F)){
    print("Because heavy labelling selected, MethylAll is set to TRUE")
    methylAll <- T
  }
  if(C13D3==T & CD3==T)
    return("ERROR - Cannot do both CD3 and C13D3 labelling")
  
  mainsequence <- toupper(mainsequence)
  
  molweight <- mw(mainsequence,T)
  
  if(propylation==T){
    molweight <- molweight+56.02621+str_count(mainsequence,"K")*56.02621
  }
  if(carbaminomethylation==T){
    molweight <- molweight+57.02146*str_count(mainsequence,"C")
  }
  if(methoxid==T){
    molweight <- molweight+15.99491*str_count(mainsequence,"M")
  }
  if(deamidation==T){
    molweight <- molweight+0.98402*(str_count(mainsequence,"N")+str_count(mainsequence,"Q"))
  }
  
  aasequence <- mainsequence
  
  if(methylAll==T | acetylation==T){
    KR <- c(rep("K",str_count(substr(mainsequence,1,nchar(mainsequence)),"K")),
            rep("R",str_count(substr(mainsequence,1,nchar(mainsequence)),"R")))
    
    KR <- as.list(KR)
    if(methylAll==T & acetylation==F)
      KR <- lapply(KR, function(x) if(x=="K") x=c("k0","k1","k2","k3") else if(x=="R") x=c("r0","r1","r2"))
    if(methylAll==F & acetylation==T)
      KR <- lapply(KR, function(x) if(x=="K") x=c("k0","ka") else if(x=="R") x=c("r0","ra"))
    if(methylAll==T & acetylation==T)
      KR <- lapply(KR, function(x) if(x=="K") x=c("k0","k1","k2","k3","ka") else if(x=="R") x=c("r0","r1","r2","ra"))
    KR <- expand.grid(KR)
    KRm <- ifelse(KR=="k0"| KR=="r0", 0,
                  ifelse(KR=="k1" | KR=="r1", 14.01565,
                         ifelse(KR=="k2",ifelse(propylation==T,-27.99491,28.03130),
                                ifelse(KR=="k3",ifelse(propylation==T,-13.97926,42.04695),
                                       ifelse(KR=="r2",28.03130,
                                              ifelse(KR=="ka" | KR=="ra",ifelse(propylation==T & KR=="ka",-14.01565,42.01056),NA))))))
    
    KRn <- apply(KR,1,paste,collapse="")
    aasequence <- KRn
    molweight <- rowSums(KRm)+molweight
  }
  
    if(C13D3==T | CD3==T){
      methGroupCount <- KR
      for(i in 1:ncol(methGroupCount)){
        methGroupCount[,i] <- substr(methGroupCount[,i],2,2)
      }
      methGroupCount <- apply(methGroupCount,2,as.numeric)
      methGroupCount <- rowSums(methGroupCount,na.rm=T)
      
      res <- matrix(NA,0,2)
      for(i in 1:length(methGroupCount)){
        temp <- matrix(NA,methGroupCount[i],2)
        if(nrow(temp)==0){
          res <- rbind(res,c(paste(aasequence[i],"0H"),molweight[i]))
        } else {
          for(j in 0:methGroupCount[i]){
            res <- rbind(res,c(paste0(aasequence[i]," ",j,"H"),molweight[i]+ifelse(C13D3==T,j*4.022185,j*3.01883)))
          }
        }
      }
      aasequence <- res[,1]; molweight <- as.numeric(res[,2])
    }
  
  if(is.na(Addition)==F){
    molweight <- molweight+Addition
  }
  mono <- molweight
  plus1 <- (molweight+1*1.007276)/1
  plus2 <- (molweight+2*1.007276)/2
  plus3 <- (molweight+3*1.007276)/3
  plus4 <- (molweight+4*1.007276)/4
  plus5 <- (molweight+5*1.007276)/5
  
  names <- c("Sequence","Monoisotopic","Mz[H+]","Mz[2H+]","Mz[3H+]","Mz[4H+]","Mz[5H+]")
  
  masses <- c(Carboaminodmethylation= "+57.02146",Oxidation= "+15.99491",Acetylation="42.01056",Deamidation="-0.98402",Propylation= "+56.02621",Methyl= "+14.01565",`13CD3`= "+4.022185", CD3= "+3.01883")
  text <- masses[which(c(carbaminomethylation,methoxid,acetylation,deamidation,propylation,methylAll,C13D3,CD3)==T)]
  if(is.na(Addition)==F) text <- c(text,Addition=as.character(Addition))
  
  results <- list(Masses=data.frame(Sequence=aasequence,Monoisoptopic=mono,`Mz[H+]`=plus1,`Mz[2H+]`=plus2,`Mz[3H+]`=plus3,`Mz[4H+]`=plus4,`Mz[5H+]`=plus5),Modifications=text,sequence=mainsequence)
  #print(text)
  return(results)
}

