library(mzR)
library(ggplot2)
library(ggrepel)
library(cowplot)

comparison <- function(data, datahead, A, B, tolerance=0.5, massShift=4){
  peaksA <- peaks(data, which(datahead$seqNum==A))
  peaksA <- as.matrix(peaksA[order(peaksA[,2],decreasing=T),])
  peaksB <- peaks(data, which(datahead$seqNum==B))
  peaksB <- as.matrix(peaksB[order(peaksB[,2],decreasing=T),])
  
  lowA <- peaksA[which(peaksA[,1]<((datahead[which(datahead$seqNum==A),"precursorDa"]/2)-40)),,drop=F]
  highA <- peaksA[which(peaksA[,1]>=((datahead[which(datahead$seqNum==A),"precursorDa"]/2)-40)),,drop=F]
  lowB <- peaksB[which(peaksB[,1]<((datahead[which(datahead$seqNum==B),"precursorDa"]/2)-40)),,drop=F]
  highB <- peaksB[which(peaksB[,1]>=((datahead[which(datahead$seqNum==B),"precursorDa"]/2)-40)),,drop=F]
  
  if(nrow(lowA)==0 | nrow(highA)==0 | nrow(lowB)==0 | nrow(highB)==0){
    return()
  }
  
  lowA <- lowA[1:min(20,nrow(lowA)),,drop=F]
  highA <- highA[1:min(20,nrow(highA)),,drop=F]
  lowB <- lowB[1:min(20,nrow(lowB)),,drop=F]
  highB <- highB[1:min(20,nrow(highB)),,drop=F]
  
  if(nrow(lowA)==0 | nrow(highA)==0 | nrow(lowB)==0 | nrow(highB)==0){
    return()
  }
  
  count <- 0
  for(i in 1:nrow(lowA)){
    if(length(which(abs(lowA[i,1]-lowB[,1])<=tolerance))!=0)
      count <- count+1
  }
  if(count>(datahead[which(datahead$seqNum==A),"precursorDa"]/400)){
    firstB <- 0
    for(i in 1:nrow(highA)){
      if(length(which(abs(highA[i,1]-highB[,1])<=tolerance))!=0)
        firstB <- firstB+1
    }
    secondB <- 0
    highB[,1] <- highB[,1] - massShift
    for(i in 1:nrow(highA)){
      if(length(which(abs(highA[i,1]-highB[,1])<=tolerance))!=0)
        secondB <- secondB+1
    }
    if(secondB > firstB)
      return(c(A,B))
  }
  #return("No match")
}

sequencepairs <- function(data, rtInterval=120, massShift=4){
  datahead <- header(data)
  datahead$precursorDa <- datahead$precursorMZ * datahead$precursorCharge - (datahead$precursorCharge * 1.007276)
  ms2 <- subset(datahead,msLevel==2)
  
  potentialSpectra <- data.frame(A_scan=integer(0), B_scan=integer(0))
  
  for(i in 1:nrow(ms2)){
    temp <- ms2[which(ms2$retentionTime>(ms2$retentionTime[i]-rtInterval) &
                        ms2$retentionTime<(ms2$retentionTime[i]+rtInterval)),]
    
    temp <- temp[which(abs(temp$precursorDa-ms2$precursorDa[i])<(massShift+0.1) &
                         abs(temp$precursorDa-ms2$precursorDa[i])>(massShift-0.1)),]
    
    scans <- data.frame(scanA=integer(0),scanB=integer(0))
    if(nrow(temp)>0){
      for(j in 1:nrow(temp)){
        hits <- comparison(data,datahead,
                           ifelse(ms2[which(ms2$seqNum==ms2$seqNum[i]),"precursorDa"]<temp$precursorDa[j],
                                  ms2$seqNum[i],temp$seqNum[j]),
                           ifelse(ms2[which(ms2$seqNum==ms2$seqNum[i]),"precursorDa"]>temp$precursorDa[j],
                                  ms2$seqNum[i],temp$seqNum[j]))
        if(length(hits)==2)
          scans <- rbind(scans,hits)
      }
    }
    if(nrow(scans)>0){
      colnames(scans) <- colnames(potentialSpectra)
      potentialSpectra <- rbind(potentialSpectra,scans) 
    }
    print(paste(i, " of ", nrow(ms2), " (",round(i/nrow(ms2)*100,2),"%)"," complete"))
    print(paste0("||",paste(rep("=",round(i/nrow(ms2)*100,0)),collapse=""),paste(rep("_",100-round(i/nrow(ms2)*100,0)),collapse=""),"||"))
  }
  return(potentialSpectra)
}

spectraplot <- function(data, scan){
  datahead <- header(data)
  datahead$precursorDa <- datahead$precursorMZ * datahead$precursorCharge - (datahead$precursorCharge * 1.007276)
  peaks <- as.matrix(peaks(data,which(datahead$seqNum==scan)))
  peaks[,2] <- peaks[,2]/max(peaks[,2])*100
  peaks <- peaks[order(peaks[,2],decreasing = T),]
  peakslow <- peaks[which(peaks[,1]<(datahead[which(datahead$seqNum==scan),"precursorDa"]/2-40)),,drop=F]
  peakshigh <- peaks[which(peaks[,1]>(datahead[which(datahead$seqNum==scan),"precursorDa"]/2-40)),,drop=F]
  
  peaks <- rbind(peakslow[1:min(20,nrow(peakslow)),],peakshigh[1:min(20,nrow(peakshigh)),])
  peaks <- as.data.frame(peaks)
  peakslabel <- peaks[!duplicated(round(peaks[,1],1)),]
  
  plot <- ggplot() + 
    geom_linerange(data=peaks,aes(x=V1,ymin=0,ymax=V2)) + 
    geom_text_repel(data=peakslabel,aes(x=V1,y=V2+2,label=as.character(round(V1,1))),color="blue") +
    geom_vline(aes(xintercept=datahead[which(datahead$seqNum==scan),"precursorDa"]/2-40),color="red",linetype="dashed") +
    ggtitle(paste("Scan = ",scan),subtitle = paste("PrecursorDa = ", round(datahead[which(datahead$seqNum==scan),"precursorDa"],2)))
  
  return(plot)
}

sideBySidePlot <- function(data,scan1,scan2){
  plot1 <- spectraplot(data,scan1)
  plot2 <- spectraplot(data,scan2)
  return(plot_grid(plot1,plot2))
}
