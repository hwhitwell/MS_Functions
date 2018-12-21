library(mzR)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(compiler)
library(Hmisc)
library(ROCR)
library(pROC)

proton <- 1.00727646687991

#Compairs two spectra based on the number of matching ions
comparison <- function(data, datahead, A, B, tolerance=0.5, massShift=4, binwidth=100, numPerBin=10){
  peaksA <- peaks(data, which(datahead$acquisitionNum==A))
  peaksA <- as.matrix(peaksA[order(peaksA[,2],decreasing=T),])
  peaksB <- peaks(data, which(datahead$acquisitionNum==B))
  peaksB <- as.matrix(peaksB[order(peaksB[,2],decreasing=T),])
  
  x <- min(peaksA[,1])
  while(x<max(peaksA[,1])){
    range <- which(peaksA[,1]>x & peaksA[,1]<(x+binwidth))
    if(length(range)!=0){
      if(length(range)>numPerBin){
        temp <- peaksA[range,][1:numPerBin,]
      } else {
        temp <- peaksA[range,]}
      peaksA <- peaksA[-range,]
      peaksA <- rbind(peaksA,temp)
    }
    x <- x+binwidth
  }
  
  x <- min(peaksB[,1])
  while(x<max(peaksB[,1])){
    range <- which(peaksB[,1]>x & peaksB[,1]<(x+binwidth))
    if(length(range)!=0){
      if(length(range)>numPerBin){
        temp <- peaksB[range,][1:numPerBin,]
      } else {
        temp <- peaksB[range,]}
      peaksB <- peaksB[-range,]
      peaksB <- rbind(peaksB,temp)
    }
    x <- x+binwidth
  }
  
  lowA <- peaksA[which(peaksA[,1]<((datahead[which(datahead$acquisitionNum==A),"precursorMZ"]/2)-40)),,drop=F]
  highA <- peaksA[which(peaksA[,1]>=(datahead[which(datahead$acquisitionNum==A),"precursorMZ"]/2)),,drop=F]
  lowB <- peaksB[which(peaksB[,1]<((datahead[which(datahead$acquisitionNum==B),"precursorMZ"]/2)-40)),,drop=F]
  highB <- peaksB[which(peaksB[,1]>=(datahead[which(datahead$acquisitionNum==B),"precursorMZ"]/2)),,drop=F]
  
  if(nrow(lowA)==0 | nrow(highA)==0 | nrow(lowB)==0 | nrow(highB)==0){
    return()
  }
  
  count <- 0
  for(i in 1:nrow(lowA)){
    if(length(which(abs(lowA[i,1]-lowB[,1])<=tolerance))!=0)
      count <- count+1
  }
  if(count>(datahead[which(datahead$acquisitionNum==A),"precursorDa"]/400)){
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

#Compairs all spectra based on the number of matching ions (uses comparison function)
sequencepairs <- function(data, rtInterval=120, massShift=4, binwidth=100, numPerBin=10){
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
                           ifelse(ms2[which(ms2$acquisitionNum==ms2$acquisitionNum[i]),"precursorDa"]<temp$precursorDa[j],
                                  ms2$acquisitionNum[i],temp$acquisitionNum[j]),
                           ifelse(ms2[which(ms2$acquisitionNum==ms2$acquisitionNum[i]),"precursorDa"]>temp$precursorDa[j],
                                  ms2$acquisitionNum[i],temp$acquisitionNum[j]),
                           binwidth=binwidth,numPerBin=numPerBin)
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

#Generates a plot from MS peak list with peak filtering
spectraplot <- function(data, scan, binwidth=100, numPerBin=10, labels=T, numLabels=NA, precursorCharge=NA){
  datahead <- header(data)
  if(is.na(precursorCharge)){
    precursorCharge <- datahead$precursorCharge
  }
  datahead$precursorDa <- datahead$precursorMZ * precursorCharge - (precursorCharge * 1.007276)
  peaks <- as.matrix(peaks(data,which(datahead$acquisitionNum==scan)))
  peaks[,2] <- peaks[,2]/max(peaks[,2])*100
  
  peaks <- binning(peaks,binwidth=binwidth,numPerBin=numPerBin)
  
  plot <- ggplot() + 
    geom_linerange(data=peaks,aes(x=V1,ymin=0,ymax=V2)) + 
    geom_vline(aes(xintercept=(datahead[which(datahead$acquisitionNum==scan),"precursorDa"]+1.007276*2)/2),color="red",linetype="dashed") +
    ggtitle(paste("Scan = ",scan),subtitle = paste("PrecursorDa = ", round(datahead[which(datahead$acquisitionNum==scan),"precursorDa"],2))) +
    ylab("Normalised Intensity") +
    xlab("m/z")
  
  if(labels){
    peakslabel <- peaks[!duplicated(round(peaks[,1],1)),]
    if(!is.na(numLabels)){
      peakslabel <- peakslabel[order(peakslabel[,2],decreasing=T),][1:numLabels,]
    }
    plot <- plot + geom_text_repel(data=peakslabel,aes(x=V1,y=V2+2,label=as.character(round(V1,1))),color="blue")
  }
  return(plot)
}

#Puts two spectra next to each other (calls spectraplot twice)
sideBySidePlot <- function(data,scan1,scan2, binwidth=100, numPerBin=10, labels=T, numLabels=NA){
  plot1 <- spectraplot(data,scan1,binwidth=binwidth,numPerBin=numPerBin, labels=labels, numLabels=numLabels)
  plot2 <- spectraplot(data,scan2,binwidth=binwidth,numPerBin=numPerBin, labels=labels, numLabels=numLabels)
  return(plot_grid(plot1,plot2))
}

#Binning function - returns "numPerBin" peaks in each bin "binwidth" along spectra
binning <- function(peaks, binwidth, numPerBin){
  peaks <- as.data.frame(peaks)
  peaks <- peaks[order(peaks[,2],decreasing=T),]
  x <- min(peaks[,1])
  while(x<max(peaks[,1])){
    range <- which(peaks[,1]>x & peaks[,1]<(x+binwidth))
    if(length(range)!=0){
      temp <- peaks[range,]
      temp <- temp[1:min(nrow(temp),numPerBin),]
      peaks <- peaks[-range,]
      peaks <- rbind(peaks,temp)
    }
    x <- x+binwidth
  }
  return(peaks)
}

#Calculates regression of matched ions and matched complementary ions in two spectra, returning paisr that match the filter criteria
matchedRegression <- function(data, datahead, A, B, binwidth=100, numPerBin=10, tolerance=0.3,complementTolerance=0.5, plot=F, Filtered=T){
  #Generate matrix of mz and intensities for each scan
  peaksA <- peaks(data, which(datahead$acquisitionNum==A))
  peaksB <- peaks(data, which(datahead$acquisitionNum==B))
  
  #Filter peak lists by top x in binwidth y
  peaksA <- binning(peaksA,binwidth=binwidth,numPerBin=numPerBin)
  peaksB <- binning(peaksB,binwidth=binwidth,numPerBin=numPerBin)
  
  #Iterate through masses in peaksA, looking for matching peaks in B. Filtering by most intense peaks within the tolerance window.
  peaksA <- peaksA[order(peaksA[,2],decreasing = T),]
  peaksB <- peaksB[order(peaksB[,2],decreasing = T),]
  colnames(peaksA) <- c("MZ","Intensity")
  colnames(peaksB) <- c("MZ","Intensity")
  hitsA <- data.frame(MZ=numeric(0), Intensity=numeric(0))
  hitsB <- data.frame(MZ=numeric(0), Intensity=numeric(0))
  for(i in 1:nrow(peaksA)){
    tempA <- peaksA[which(abs((peaksA[i,1]-peaksA[,1]))<=tolerance),]
    if(nrow(tempA)!=0){
      tempA <- tempA[1,]
      tempB <- peaksB[which(abs((tempA[,1]-peaksB[,1]))<=tolerance),]
      if(nrow(tempB)!=0){
        hitsA <- rbind(hitsA,tempA)
        hitsB <- rbind(hitsB,tempB[1,])
      }
    }
  }
  if(nrow(hitsA)!=0 & nrow(hitsB)!=0){
    hits <- cbind(hitsA,hitsB[,2])
    if(nrow(hits)>datahead[which(datahead$acquisitionNum==A),"precursorDa"]/200 &
       nrow(hits)>3){
      
      colnames(hits) <- c("MZ","Intensity_A","Intensity_B")
      
      #calculate Peasons produce moment correlation coefficient, printing the plots and ultimatly returning the correlation score.
      correlation <- with(hits,cor.test(Intensity_A,Intensity_B))
      if(Filtered==T & correlation$estimate[[1]]<0.6){
        return(data.frame(A=numeric(0),B=numeric(0)))
      }
      
      if(plot==T){
        plot1 <- ggplot(hits, aes(x=Intensity_A/max(Intensity_A)*100,y=Intensity_B/max(Intensity_B)*100)) +
          geom_point(aes(colour=MZ)) +
          geom_smooth(method="lm") +
          xlab(paste0("normalised intensity (",A,")")) +
          ylab(paste0("normalised intensity (",B,")")) +
          ggtitle("Matched Ions",subtitle=paste0("PPMCC=",round(correlation$estimate,3)," P-value=",round(correlation$p.value,5)))
        
        print(plot1)}
      
      #Predict the mass of complementary masses (precursor mass[M+H] - matched ion mass) for A
      predictedComplementary <- hits[,1]
      predictedComplementary <- datahead[which(datahead$acquisitionNum==A),"precursorDa"] + 2*1.007276 - predictedComplementary
      #Look for the predicted ions in A within tolerance of 0.3Da
      complementaryHitsA <- data.frame("MZ"=numeric(0),"Intensity"=numeric(0))
      for(i in 1:length(predictedComplementary)){
        tempA <- peaksA[which(abs(predictedComplementary[i]-peaksA[,1])<=complementTolerance),]
        if(nrow(tempA)!=0){
          complementaryHitsA <- rbind(complementaryHitsA,tempA[1,])
        } 
      }
      
      if(Filtered==T & nrow(complementaryHitsA)<=2){
        return(data.frame(A=numeric(0),B=numeric(0)))
      }
      
      #Match these against B with the appropriate mass adjustment
      delta <- ifelse(datahead[which(datahead$acquisitionNum==A),"precursorDa"]<datahead[which(datahead$acquisitionNum==B),"precursorDa"],4,-4)
      complementaryA <- data.frame("MZ"=numeric(0),"Intensity"=numeric(0))
      complementaryB <- data.frame("MZ"=numeric(0),"Intensity"=numeric(0))
      for(i in 1:nrow(complementaryHitsA)){
        tempB <- peaksB[which(abs(peaksB[,1]-(complementaryHitsA[i,1]+delta))<=complementTolerance),]
        if(nrow(tempB)!=0){
          complementaryA <- rbind(complementaryA, complementaryHitsA[i,])
          complementaryB <- rbind(complementaryB, tempB[1,])
        }
      }
      
      if(Filtered==T & (nrow(complementaryA)<=2 & nrow(complementaryB)<=2)){
        return(data.frame(A=numeric(0),B=numeric(0)))
      }
      
      if(nrow(complementaryA)>2 & nrow(complementaryB)>2){
        results <- cbind(complementaryA,complementaryB)
        colnames(results) <- c("MZ_A","Intensity_A","MZ_B","Intensity_B")
        
        correlation2 <- with(results,cor.test(Intensity_A,Intensity_B))
        
        if(Filtered==T & (correlation2$estimate[[1]]<0.7 | correlation2$p.value>0.05 | is.na(correlation2$p.value)==T)){
          return(data.frame(A=numeric(0),B=numeric(0)))
        }
        
        complementaryResults <- c(correlation2$estimate[[1]],correlation2$p.value)
        if(plot==T){
          plot2 <- ggplot(results, aes(Intensity_A/max(Intensity_A)*100,Intensity_B/max(Intensity_B)*100)) +
            geom_point(aes(colour=MZ_A)) +
            geom_smooth(method="lm") +
            ylab(paste0("Normalised intensity (",B,")")) +
            xlab(paste0("Normalised intensity (",A,")")) +
            ggtitle("Complementary Matches",subtitle=paste0("PPMCC=",round(correlation2$estimate,3)," P-value=",round(correlation2$p.value,5))) 
          print(plot2)
        }
      } else {complementaryResults <- c(NA,NA,NA,NA)
      results <- data.frame("MZ_A"=NA,"Intensity_A"=NA,"MZ_B"=NA,"Intensity_B"=NA)}
      
      colnames(results) <- c("MZ_A","Intensity_A","MZ_B","Intensity_B")
      
      results <- data.frame("A"=A,"B"=B,
                            "PPMCC"=correlation$estimate[[1]],
                            "P-Val"=correlation$p.value,
                            "ComplementaryPPMCC"=complementaryResults[1],
                            "ComplementaryPVal"=complementaryResults[2],
                            "NumHits"=nrow(hits),
                            "HitIntensity_A"=mean(hits$Intensity_A),
                            "HitIntensity_B"=mean(hits$Intensity_B),
                            "ComplementaryHits"=nrow(results),
                            "ComplementaryIntensity_A"=mean(results$Intensity_A),
                            "ComplementaryIntensity_B"=mean(results$Intensity_B))
      return(results)
    } else return(data.frame(A=numeric(0),B=numeric(0)))
  } else return(data.frame(A=numeric(0),B=numeric(0)))
}
matchedRegression <- cmpfun(matchedRegression)

#Calculates precursor mass from mxXML header
precursorMass <- function(datahead){
  precusorDa <- datahead$precursorMZ * datahead$precursorCharge - (datahead$precursorCharge * 1.007276)
  return(precusorDa)
}

#Compares multiple  spectra using matchedRegression
regressionPairs <- function(data, rtInterval=120, massShift=4, binwidth=100, numPerBin=10, tolerance=0.3, complementTolerance=0.5, plot=F, suppressProgressBar=F, Filtered=T){
  datahead <- header(data)
  datahead$precursorDa <- precursorMass(datahead)
  ms2 <- subset(datahead,msLevel==2)
  
  results <- data.frame("A"=integer(0),
                        "B"=integer(0),
                        "PPMCC"=numeric(0),
                        "P-Val"=numeric(0),
                        "ComplementaryPPMCC"=numeric(0),
                        "ComplementaryPVal"=numeric(0),
                        "NumHits"=numeric(0),
                        "HitIntensity_A"=numeric(0),
                        "HitIntensity_B"=numeric(0),
                        "ComplementaryHits"=numeric(0),
                        "ComplementaryIntensity_A"=numeric(0),
                        "ComplementaryIntensity_B"=numeric(0))
  
  for(i in 1:nrow(ms2)){
    temp <- ms2[which(ms2$retentionTime>(ms2$retentionTime[i]-rtInterval) &
                        ms2$retentionTime<(ms2$retentionTime[i]+rtInterval)),]
    
    temp <- temp[which(abs(temp$precursorDa-ms2$precursorDa[i])<(massShift+2) &
                         abs(temp$precursorDa-ms2$precursorDa[i])>(massShift-2)),]
    
    tempresults <- data.frame("A"=integer(0),
                              "B"=integer(0),
                              "PPMCC"=numeric(0),
                              "P-Val"=numeric(0),
                              "ComplementaryPPMCC"=numeric(0),
                              "ComplementaryPVal"=numeric(0),
                              "NumHits"=numeric(0),
                              "HitIntensity_A"=numeric(0),
                              "HitIntensity_B"=numeric(0),
                              "ComplementaryHits"=numeric(0),
                              "ComplementaryIntensity_A"=numeric(0),
                              "ComplementaryIntensity_B"=numeric(0))
    if(nrow(temp)>0){
      for(j in 1:nrow(temp)){
        hits <- matchedRegression(data,datahead,ms2$acquisitionNum[i],temp$acquisitionNum[j],binwidth=binwidth,numPerBin=numPerBin,tolerance=tolerance, complementTolerance = complementTolerance,plot=plot, Filtered=Filtered)
        if(nrow(hits)>0)
          tempresults <- rbind(tempresults,hits)
      }
    }
    if(nrow(tempresults)>0){
      results <- rbind(results,tempresults) 
    }
    if(suppressProgressBar==F){
      print(paste(i, " of ", nrow(ms2), " (",round(i/nrow(ms2)*100,2),"%)"," complete"))
      print(paste0("||",paste(rep("=",round(i/nrow(ms2)*100,0)),collapse=""),paste(rep("_",100-round(i/nrow(ms2)*100,0)),collapse=""),"||")) 
    }
  }
  return(results)
}
regressionPairs <- cmpfun(regressionPairs)

#Called by matchedRegression.
matchedComplemetaryPeaks <- function(data, A, B, massShift=4, binwidth=100, numPerBin=10, tolerance=1, plot=F, Filtered){
  #Generate table of header data
  datahead <- header(data)
  datahead$precursorDa <- datahead$precursorMZ * datahead$precursorCharge - (datahead$precursorCharge * 1.007276)
  
  #Generate matrix of mz and intensities for each scan, ordered by normalised intensity
  peaksA <- peaks(data, which(datahead$acquisitionNum==A))
  peaksA[,2] <- peaksA[,2]/max(peaksA[,2])*100
  peaksB <- peaks(data, which(datahead$acquisitionNum==B))
  peaksB[,2] <- peaksB[,2]/max(peaksB[,2])*100
  
  #Filter peak lists by top x in binwidth y
  peaksA <- binning(peaksA,binwidth=binwidth,numPerBin=numPerBin)
  peaksB <- binning(peaksB,binwidth=binwidth,numPerBin=numPerBin)
  
  #Iterate through masses in peaksA, looking for matching peaks in B. Filtering by most intense peaks within the tolerance window.
  peaksA <- peaksA[order(peaksA[,2],decreasing = T),]
  peaksB <- peaksB[order(peaksB[,2],decreasing = T),]
  colnames(peaksA) <- c("MZ","Intensity")
  colnames(peaksB) <- c("MZ","Intensity")
  hitsA <- data.frame(MZ=numeric(0), Intensity=numeric(0))
  hitsB <- data.frame(MZ=numeric(0), Intensity=numeric(0))
  for(i in 1:nrow(peaksA)){
    tempA <- peaksA[which(abs((peaksA[i,1]-peaksA[,1]))<=tolerance),]
    if(nrow(tempA)!=0){
      tempA <- tempA[1,]
      tempB <- peaksB[which(abs((tempA[,1]-peaksB[,1]))<=tolerance),]
      if(nrow(tempB)!=0){
        hitsA <- rbind(hitsA,tempA)
        hitsB <- rbind(hitsB,tempB[1,])
      }
    }
  }
  if(nrow(hitsA)!=0 & nrow(hitsB)!=0){
    hits <- cbind(hitsA,hitsB)
    if(nrow(hits)>datahead[which(datahead$acquisitionNum==A),"precursorDa"]/200 &
       nrow(hits)>3){
      
      colnames(hits) <- c("MZ_A","Intensity_A","MZ_B","Intensity_B")
      
      #Predict the mass of complementary masses (precursor mass[M+H] - matched ion mass) for A
      predictedComplementary <- hits[,1]
      predictedComplementary <- datahead[which(datahead$acquisitionNum==A),"precursorDa"] + 1.007276 - predictedComplementary
      #Look for the predicted ions in A within tolerance of 0.3Da
      complementaryHitsA <- data.frame("MZ"=numeric(0),"Intensity"=numeric(0))
      for(i in 1:length(predictedComplementary)){
        tempA <- peaksA[which(abs(predictedComplementary[i]-peaksA[,1])<=tolerance),]
        if(nrow(tempA)!=0){
          complementaryHitsA <- rbind(complementaryHitsA,tempA[1,])
        }
      }
      #Match these against B with the appropriate mass adjustment
      delta <- ifelse(datahead[which(datahead$acquisitionNum==A),"precursorDa"]<datahead[which(datahead$acquisitionNum==B),"precursorDa"],4,-4)
      complementaryA <- data.frame("MZ"=numeric(0),"Intensity"=numeric(0))
      complementaryB <- data.frame("MZ"=numeric(0),"Intensity"=numeric(0))
      for(i in 1:nrow(complementaryHitsA)){
        tempB <- peaksB[which(abs(peaksB[,1]-(complementaryHitsA[i,1]+delta))<=tolerance),]
        if(nrow(tempB)!=0){
          complementaryA <- rbind(complementaryA, complementaryHitsA[i,])
          complementaryB <- rbind(complementaryB, tempB[1,])
        }
      }
      if(nrow(complementaryA)>2 & nrow(complementaryB)>2){
        results <- cbind(complementaryA,complementaryB)
        colnames(results) <- c("MZ_A","Intensity_A","MZ_B","Intensity_B")
        
        correlation <- with(results,cor.test(Intensity_A,Intensity_B))
        
        if(plot==T){
          plot <- ggplot(results, aes(Intensity_A,Intensity_B)) +
            geom_point(aes(colour=MZ_A)) +
            geom_smooth(method="lm") +
            ylab(paste0(B," Intensity")) +
            xlab(paste0(A," Intensity")) +
            ggtitle("",subtitle=paste0("PPMCC=",round(correlation$estimate,3)," P.Val=",round(correlation$p.value,5))) 
          print(plot)
        }
        return(data.frame("A"=A,"B"=B,"PPMCC"=correlation$estimate[[1]],"P-Val"=correlation$p.value)) 
      } else return(data.frame(A=numeric(0),B=numeric(0)))
    } else return(data.frame(A=numeric(0),B=numeric(0)))
  } else return(data.frame(A=numeric(0),B=numeric(0)))
}

#Identifies peaks in spectra.
peakpicker <- function(peaklist, threshold=2, normalisedColumn=2){
  peak <- ifelse(peaklist[,normalisedColumn]<threshold,0,1)
  if(peak[length(peak)]==1){
    peak <- c(peak,0)
  }
  begin <- vector()
  end <- vector()
  i=1
  while(i<length(peak)){
    if(peak[i]==1){
      begin <- c(begin,i)
      for(j in (i+1):length(peak)){
        if(peak[j]!=1){
          end <- c(end,(j-1))
          i <- j+1
          break()
        }
      }
    } else
      i <- i+1
  }
  peaklist <- data.frame(start=begin,end=end)
  return(peaklist)
}

#Bins peaks into an equally spaced matrix - for convolve function
peakBinning <- function(peakList,peakPickThreshold=2,binThreshold=0.004,plot=T){
  if(peakList[nrow(peakList),2]!=0){
    peakList <- rbind(peakList,c(peakList[nrow(peakList),1]+0.01,0))
  }
  pickPeaks <- peakpicker(peakList,threshold=peakPickThreshold)
  min <- seq(min(peakList[,1]),max(peakList[,1])-binThreshold,binThreshold)
  binned <- data.frame(min,
                       max=min+binThreshold,
                       Boolian=rep(F,length(min)),
                       Intensity=rep(0,length(min)))
  
  for(i in 1:nrow(peakList)){
    temp_min <- peakList[,1][pickPeaks$start[i]]
    temp_max <- peakList[,1][pickPeaks$end[i]]
    binned$Boolian[which(binned$max>=temp_min & binned$min<=temp_max)] <- T
  }
  rm(temp_min,temp_max)
  
  trueRows <- which(binned$Boolian==T)
  for(i in 1:length(trueRows)){
    temp_intensity <- peakList[,2][which(peakList[,1]>=binned[trueRows[i],"min"] &
                                           peakList[,1]<=binned[trueRows[i],"max"])]
    temp_intensity <- mean(temp_intensity)
    binned$Intensity[trueRows[i]] <- temp_intensity
  }
  
  rowsNA <- which(is.na(binned$Intensity))
  if(length(rowsNA)!=0){
    for(i in 1:length(rowsNA)){
      minMZ <- mean(binned$min[rowsNA[i]-1],binned$max[rowsNA[i]-1])
      minInt <- binned$Intensity[rowsNA[i]-1]
      maxInt <- 0
      row=rowsNA[i]
      while(maxInt==0){
        row = row+1
        maxInt <- binned$Intensity[row]
        if(is.na(maxInt)){
          maxInt <- 0
        }
        if(maxInt>0){
          maxMZ <- mean(binned$min[row],binned$max[row])
        }
      }
      temp <- data.frame(mz=c(minMZ,maxMZ),int=c(minInt,maxInt))
      temp_mz <- mean(binned$min[rowsNA[i]],binned$max[rowsNA[i]])
      temp_intensity <- lm(int~mz,temp)
      temp_intensity <- predict(temp_intensity,data.frame(mz=temp_mz))
      binned$Intensity[rowsNA[i]] <- temp_intensity
    } 
  }
  if(plot){
    print(ggplot(binned,aes(x=min,ymin=0,ymax=Intensity)) + geom_linerange())
  }
  return(binned)
}

#Identifies if a precursorMZ has a heavy pair in an MS1 spectra and if so performs convovle function.
convolveFunction <- function(raw,datahead,scan,heavyShift=4.02218, binThreshold=0.004, plot=T){
  scanType <- datahead[which(datahead$acquisitionNum==scan),"msLevel"]
  if(scanType!=2){
    results <- data.frame(MS2Scan=scan,
                          Type="Scan is not MS2",
                          Mz=precursorMz,
                          Charge=precursorCharge,
                          RT=datahead[which(datahead$acquisitionNum==scan),"retentionTime"]/60,
                          Intensity=datahead[which(datahead$acquisitionNum==scan),"precursorIntensity"],
                          MaxConvolve=NA,
                          ExpectedBinConvolve=NA,
                          sequentialConvolve=NA,
                          DeltaScore=NA,
                          Threshold=NA,
                          NumOfPeaks=NA,
                          Distance=NA,
                          DeltaDistance=NA,
                          glmScore=NA,
                          glmPredict=NA)
    return(results)
  }
  precursorScan <- datahead[which(datahead$acquisitionNum==scan),"precursorScanNum"]
  precursorMz <- datahead[which(datahead$acquisitionNum==scan),"precursorMZ"]
  precursorCharge <- datahead[which(datahead$acquisitionNum==scan),"precursorCharge"]
  peakList <- peaks(raw,which(datahead$acquisitionNum==precursorScan))
  
  #filter mass list to return the precursor isotope series with enough width to keep heavy/light series
  peakList <- peakList[which(peakList[,1]<precursorMz+4 & peakList[,1]>precursorMz-4),]
  
  if(nrow(peakList)==0){
    results <- data.frame(MS2Scan=scan,
                          Type="No peaks in range",
                          Mz=precursorMz,
                          Charge=precursorCharge,
                          RT=datahead[which(datahead$acquisitionNum==scan),"retentionTime"]/60,
                          Intensity=datahead[which(datahead$acquisitionNum==scan),"precursorIntensity"],
                          MaxConvolve=NA,
                          ExpectedBinConvolve=NA,
                          sequentialConvolve=NA,
                          DeltaScore=NA,
                          Threshold=NA,
                          NumOfPeaks=NA,
                          Distance=NA,
                          DeltaDistance=NA,
                          glmScore=NA,
                          glmPredict=NA)
    return(results)
  }
  
  #Normalised peakList
  peakList[,2] <- peakList[,2]/max(peakList[,2])*100
  
  #Calculate average isotopic m/z within each peak
  pickedPeaks <- peakpicker(peakList,7)
  averageMz <- matrix(NA,ncol=2,nrow=nrow(pickedPeaks))
  for(i in 1:nrow(averageMz)){
    averageMz[i,2] <- max(peakList[pickedPeaks[i,1]:pickedPeaks[i,2],2])
    averageMz[i,1] <- weighted.mean(peakList[pickedPeaks[i,1]:pickedPeaks[i,2],1],
                                    peakList[pickedPeaks[i,1]:pickedPeaks[i,2],2])
  }
  
  #identify precursor Mz and check if it is the first in the series
  averagePrecursorMz <- averageMz[which.min(abs(averageMz[,1]-precursorMz)),1]
  firstInSeries <- ifelse(which(averageMz[,1]==averagePrecursorMz)==1,T,F)
  while(firstInSeries==F){
    for(i in (which(averageMz[,1]==averagePrecursorMz)-1):1){
      if(averageMz[i,2]>averageMz[which(averageMz[,1]==averagePrecursorMz),2] &
         round((averagePrecursorMz-averageMz[i,1])*precursorCharge,0)==1){
        averagePrecursorMz <- averageMz[i,1]
      } else {
        break()
      }
    }
    firstInSeries <- T
  }
  
  #Check for a peak at the appropriate mass shift above or below.
  #If there is one eitherside, then it must be a multiply methylated ion.
  #If it is a multiply methylated ion, take selected precursor and the largest adjacent.
  #Generate heavy and light spectra for comparing with convolve function.
  matchedPair <- averageMz[which(abs(abs(averageMz[,1]-averagePrecursorMz)-heavyShift/precursorCharge)<0.01),,drop=F]
  if(nrow(matchedPair)>0){
    if(nrow(matchedPair)==1){
      precursorType <- ifelse(precursorMz>matchedPair[,1],"Heavy","Light")
    } else if(nrow(matchedPair)==2){
      if(matchedPair[1,1]<averagePrecursorMz &
         matchedPair[2,1]>averagePrecursorMz){
        precursorType <- ifelse(matchedPair[1,2]>matchedPair[2,2],"Heavy","Light")
      } else {
        precursorType <- "Unknown"
      }
    } else {
      precursorType <- "Unknown" 
    }
  } else {
    precursorType <- "Unknown"
  }
  rm(matchedPair)
  
  if(precursorType=="Unknown"){
    results <- data.frame(MS2Scan=scan,
                          Type="No heavy/light pair found",
                          Mz=precursorMz,
                          Charge=precursorCharge,
                          RT=datahead[which(datahead$acquisitionNum==scan),"retentionTime"]/60,
                          Intensity=datahead[which(datahead$acquisitionNum==scan),"precursorIntensity"],
                          MaxConvolve=NA,
                          ExpectedBinConvolve=NA,
                          sequentialConvolve=NA,
                          DeltaScore=NA,
                          Threshold=NA,
                          NumOfPeaks=NA,
                          Distance=NA,
                          DeltaDistance=NA,
                          glmScore=NA,
                          glmPredict=NA)
    return(results)
  }
  
  if(precursorType=="Light"){
    lightstart <- averagePrecursorMz-0.1
    lightend <- averagePrecursorMz + heavyShift/precursorCharge-0.1
    
    light <- peakList[which(peakList[,1]>=lightstart &
                              peakList[,1]<lightend),,drop=F]
    
    heavystart <- averagePrecursorMz+heavyShift/precursorCharge-0.1
    heavyend <- averagePrecursorMz+2*heavyShift/precursorCharge+0.1
    
    heavy <- peakList[which(peakList[,1]>=heavystart &
                              peakList[,1]<heavyend),,drop=F]
    
    expectedLightBin <- (averagePrecursorMz-lightstart)/binThreshold
    expectedHeavyBin <- (averagePrecursorMz + heavyShift/precursorCharge - heavystart)/binThreshold
  } else {
    if(precursorType=="Heavy"){
      lightstart <- averagePrecursorMz-heavyShift/precursorCharge-0.1
      lightend <- averagePrecursorMz-0.1
      
      light <- peakList[which(peakList[,1]>=lightstart &
                                peakList[,1]<lightend),,drop=F]
      
      heavystart <- averagePrecursorMz-0.1
      heavyend <- averagePrecursorMz + heavyShift/precursorCharge + 0.1
      
      heavy <- peakList[which(peakList[,1]>=heavystart &
                                peakList[,1]<heavyend),,drop=F]
      
      expectedLightBin <- (averagePrecursorMz - heavyShift/precursorCharge - lightstart)/binThreshold
      expectedHeavyBin <- (averagePrecursorMz - heavystart)/binThreshold
    }
  }
  
  if(nrow(light)==0 | nrow(heavy)==0 | sum(light[,2])==0 | sum(heavy[,2])==0){
    results <- data.frame(MS2Scan=scan,
                          Type="No heavy/light pair found",
                          Mz=precursorMz,
                          Charge=precursorCharge,
                          RT=datahead[which(datahead$acquisitionNum==scan),"retentionTime"]/60,
                          Intensity=datahead[which(datahead$acquisitionNum==scan),"precursorIntensity"],
                          MaxConvolve=NA,
                          ExpectedBinConvolve=NA,
                          sequentialConvolve=NA,
                          DeltaScore=NA,
                          Threshold=NA,
                          NumOfPeaks=NA,
                          Distance=NA,
                          DeltaDistance=NA,
                          glmScore=NA,
                          glmPredict=NA)
    return(results)
  }
  
  light[,2] <- light[,2]/max(light[,2])*100
  heavy[,2] <- heavy[,2]/max(heavy[,2])*100
  
  lightBinned <- peakBinning(light,peakPickThreshold=2,binThreshold = binThreshold,plot=F)
  lightBinned$Intensity <- lightBinned$Intensity/max(lightBinned$Intensity)*100
  heavyBinned <- peakBinning(heavy,peakPickThreshold=2,binThreshold = binThreshold,plot=F)
  heavyBinned$Intensity <- heavyBinned$Intensity/max(heavyBinned$Intensity)*100
  
  lightIntensity <- lightBinned$Intensity
  heavyIntensity <- heavyBinned$Intensity
  
  if(length(lightIntensity)<31 | length(heavyIntensity)<31){
    results <- data.frame(MS2Scan=scan,
                          Type="No heavy/light pair found",
                          Mz=precursorMz,
                          Charge=precursorCharge,
                          RT=datahead[which(datahead$acquisitionNum==scan),"retentionTime"]/60,
                          Intensity=datahead[which(datahead$acquisitionNum==scan),"precursorIntensity"],
                          MaxConvolve=NA,
                          ExpectedBinConvolve=NA,
                          sequentialConvolve=NA,
                          DeltaScore=NA,
                          Threshold=NA,
                          NumOfPeaks=NA,
                          Distance=NA,
                          DeltaDistance=NA,
                          glmScore=NA,
                          glmPredict=NA)
    return(results)
  }
  
  if(length(lightIntensity)<length(heavyIntensity)){
    lightIntensity <- c(lightIntensity,rep(0,length(heavyIntensity)-length(lightIntensity)))
  } else {
    heavyIntensity <- c(heavyIntensity,rep(0,length(lightIntensity)-length(heavyIntensity)))
  }
  
  convolveResults <- convolve(lightIntensity,heavyIntensity)
  
  if(plot){
    lightBin <- ggplot(lightBinned,aes(x=min,ymin=0,ymax=Intensity)) + geom_linerange() + ggtitle("Light Binned")
    heavyBin <- ggplot(heavyBinned,aes(x=min,ymin=0,ymax=Intensity)) + geom_linerange() + ggtitle("Heavy Binned")
    lightheavybined <- plot_grid(lightBin,heavyBin)
    y <- data.frame(x=1:length(convolveResults),y=convolveResults)
    convolved <- ggplot(y, aes(x,y)) +
      geom_point() +
      ggtitle("Convolve Function",subtitle=paste0("MS2Scan=",scan," Type=",precursorType, " PrecursorMz=",round(precursorMz,3),"\nscore=",round(max(convolveResults),0)," PrecursorCharge=",precursorCharge," RT=",round(datahead[which(datahead$acquisitionNum==scan),"retentionTime"]/60,2)))
    print(plot_grid(lightheavybined,convolved,nrow=2))
  }
  
  maxProductScore <- max(convolveResults)
  expectedBinScore <- max(convolveResults[c(1:30,(length(convolveResults)-30):length(convolveResults))])
  productScoreThreshold <- mean(convolveResults)+sd(convolveResults)
  convolveResults <- data.frame(MZ=1:length(convolveResults),Int=convolveResults)
  convolvePeaks <- peakpicker(convolveResults,productScoreThreshold)
  numPeak <- nrow(convolvePeaks)
  peakMaxBin <- matrix(NA,numPeak)
  for(i in 1:nrow(peakMaxBin)){
    temp <- convolveResults[convolvePeaks[i,1]:convolvePeaks[i,2],]
    peakMaxBin[i] <- temp[which.max(temp[,2]),1]
  }
  distance <- mean(binThreshold*(peakMaxBin[2:length(peakMaxBin)] - peakMaxBin[1:(length(peakMaxBin)-1)]))
  
  glmScore <- sum(c(maxProductScore,expectedBinScore,productScoreThreshold) * c(0.02642976 ,0.01477461,-0.05876565)) -1208.702  
  glmScore <- exp(glmScore)
  glmScore <- ifelse(glmScore==Inf, 1, glmScore/(1+glmScore))
  glmPredict <- ifelse(glmScore<0.5,0,1)
  
  if(glmPredict==1){
    precursorScans <- datahead[which(datahead$msLevel==1),"acquisitionNum"]
    precursorStart <- which(precursorScans==precursorScan) - 1
    count <- 0
    sequentialScore <- expectedBinScore
    #cutoff <- expectedBinScore-0.5*productScoreThreshold
    cutoff <- expectedBinScore-sd(convolveResults[,2])
    if(cutoff<productScoreThreshold){
      count <- 1
    } else {
      while(sequentialScore>=cutoff){
        count <- count+1
        sequentialPeaks <- peaks(raw,which(datahead$acquisitionNum==precursorScans[precursorStart]))
        sequentialLight <- sequentialPeaks[which(sequentialPeaks[,1]>=lightstart &
                                                   sequentialPeaks[,1]<lightend),,drop=F]
        sequentialHeavy <- sequentialPeaks[which(sequentialPeaks[,1]>=heavystart &
                                                   sequentialPeaks[,1]<heavyend),,drop=F]
        
        if(nrow(sequentialLight)==0 | nrow(sequentialHeavy)==0 | sum(sequentialLight[,2])==0 | sum(sequentialHeavy[,2])==0){
          sequentialScore <- 0
        } else {
          sequentialLight[,2] <- sequentialLight[,2]/max(sequentialLight[,2])*100
          sequentialHeavy[,2] <- sequentialHeavy[,2]/max(sequentialHeavy[,2])*100
          
          sequentialLight <- peakBinning(sequentialLight,peakPickThreshold=2,binThreshold = binThreshold,plot=F)
          sequentialHeavy <- peakBinning(sequentialHeavy,peakPickThreshold=2,binThreshold = binThreshold,plot=F)
          
          sequentialLight <- sequentialLight[,4]/max(sequentialLight[,4])*100
          sequentialHeavy <- sequentialHeavy[,4]/max(sequentialHeavy[,4])*100
          
          if(length(sequentialLight)<length(sequentialHeavy)){
            sequentialLight <- c(sequentialLight,rep(0,length(sequentialHeavy)-length(sequentialLight)))
          } else {
            sequentialHeavy <- c(sequentialHeavy ,rep(0,length(sequentialLight)-length(sequentialHeavy)))
          }
          
          sequentialConvolveResults <- convolve(sequentialLight,sequentialHeavy)
          sequentialScore <- ifelse(length(sequentialConvolveResults)<60,
                                    0,
                                    max(sequentialConvolveResults[c(1:30,(length(sequentialConvolveResults)-30):length(sequentialConvolveResults))]))
          
          precursorStart <- precursorStart-1
        }
      }
      precursorStart <- which(precursorScans==precursorScan) + 1
      sequentialScore <- expectedBinScore
      while(sequentialScore>=cutoff){
        count <- count+1
        sequentialPeaks <- peaks(raw,which(datahead$acquisitionNum==precursorScans[precursorStart]))
        sequentialLight <- sequentialPeaks[which(sequentialPeaks[,1]>=lightstart &
                                                   sequentialPeaks[,1]<lightend),,drop=F]
        sequentialHeavy <- sequentialPeaks[which(sequentialPeaks[,1]>=heavystart &
                                                   sequentialPeaks[,1]<heavyend),,drop=F]
        
        if(nrow(sequentialLight)==0 | nrow(sequentialHeavy)==0 | sum(sequentialLight[,2])==0 | sum(sequentialHeavy[,2])==0){
          sequentialScore <- 0
        } else {
          sequentialLight[,2] <- sequentialLight[,2]/max(sequentialLight[,2])*100
          sequentialHeavy[,2] <- sequentialHeavy[,2]/max(sequentialHeavy[,2])*100
          
          sequentialLight <- peakBinning(sequentialLight,peakPickThreshold=2,binThreshold = binThreshold,plot=F)
          sequentialHeavy <- peakBinning(sequentialHeavy,peakPickThreshold=2,binThreshold = binThreshold,plot=F)
          
          sequentialLight <- sequentialLight[,4]/max(sequentialLight[,4])*100
          sequentialHeavy <- sequentialHeavy[,4]/max(sequentialHeavy[,4])*100
          
          if(length(sequentialLight)<length(sequentialHeavy)){
            sequentialLight <- c(sequentialLight,rep(0,length(sequentialHeavy)-length(sequentialLight)))
          } else {
            sequentialHeavy <- c(sequentialHeavy ,rep(0,length(sequentialLight)-length(sequentialHeavy)))
          }
          
          sequentialConvolveResults <- convolve(sequentialLight,sequentialHeavy)
          sequentialScore <- ifelse(length(sequentialConvolveResults)<60,
                                    0,
                                    max(sequentialConvolveResults[c(1:30,(length(sequentialConvolveResults)-30):length(sequentialConvolveResults))]))
          
          precursorStart <- precursorStart+1
        }
      }
    }
    sequentialConvolve <- count-1
  } else {
    sequentialConvolve <- NA
  }
  
  results <- data.frame(MS2Scan=scan,
                        Type=precursorType,
                        Mz=precursorMz,
                        Charge=precursorCharge,
                        RT=datahead[which(datahead$acquisitionNum==scan),"retentionTime"]/60,
                        Intensity=datahead[which(datahead$acquisitionNum==scan),"precursorIntensity"],
                        MaxConvolve=maxProductScore,
                        ExpectedBinConvolve=expectedBinScore,
                        sequentialConvolve=sequentialConvolve,
                        DeltaScore=maxProductScore-expectedBinScore,
                        Threshold=productScoreThreshold,
                        NumOfPeaks=numPeak,
                        Distance=distance,
                        DeltaDistance=abs(distance-1/precursorCharge),
                        glmScore=glmScore,
                        glmPredict=glmPredict)
  
  return(results)
}

####OLEG's SCRIPT
auc_spec_sens <- function(data,Spec) #Calculation of AUC, Spec, Sens for data with 1st column outcome and other column - predictors
{colnames(data)[1]="output"
mylogit=glm(output~.,data=data,family="binomial") #ERROR OCCURING HERE
temp <- tryCatch(predict(mylogit,type="response",se=TRUE),error=function(e) 1) #This puts result to 0 if prediction fails
if (is.numeric(temp)==T){
  res <- list(auc=0,spec=0,sens=0)
} else{
  prGLM <- predict(mylogit,type="response",se=TRUE)
  Tmp=cbind.data.frame(data[,1],prGLM$fit)
  auc=somers2(Tmp[,2],Tmp[,1])[1]
  pred1<-prediction(Tmp[,2],Tmp[,1])
  perf1<-performance(pred1,"tpr","fpr")
  ind=which((1-perf1@x.values[[1]])>=Spec)
  num=ind[length(ind)]
  sens=perf1@y.values[[1]][num]
  ind=which(perf1@y.values[[1]]==sens)
  spec=max(1-perf1@x.values[[1]][ind])
  res=list(auc=auc,spec=spec,sens=sens)
}
return(res)
}


####OLEG's SCRIPT
best_comb_linear <- function(Data, y, x, num_ind, type, Spec)
  #looking for parameters for best linear combination: Data - datafile, y - number of column with outcome, x - columns with predictors, x_max - maximal number of variables in the final model, type - "AUC" or "Sens" (in this case "Spec" value will be used)
{sens_max=0
spec_max=0
auc_max=0
if(type=="All"){
  sens_max <- vector()
  spec_max <- vector()
  auc_max <- vector()
}
k<-0; ind_max<-rep(0, num_ind)

if(type=="All"){
  ind_max <- vector()
}

i<-num_ind	#looking through all combinations of "i" variables
tmp=combn(x,i)
for (j in 1:length(tmp[1,]))
{data=data.frame(Data[,c(y,tmp[,j])])
ROC_calc=auc_spec_sens(data,Spec)

if (type=="AUC")
  if (ROC_calc$auc>auc_max)
  {auc_max=ROC_calc$auc
  sens_max=ROC_calc$sens
  spec_max=ROC_calc$spec
  ind_max=tmp[,j]
  }
if (type=="Sens")
  if (ROC_calc$sens>sens_max)
  {auc_max=ROC_calc$auc
  sens_max=ROC_calc$sens
  spec_max=ROC_calc$spec
  ind_max=tmp[,j]
  }	
if(type=="All"){
  auc_max <- c(auc_max,ROC_calc$auc)
  sens_max <- c(sens_max,ROC_calc$sens)
  spec_max <- c(spec_max,ROC_calc$spec)
  ind_max <- c(ind_max,paste(tmp[,j],collapse=","))
}
}

res<-list(auc_max=auc_max, sens_max=sens_max, spec_max=spec_max, ind_max=ind_max)
return(res)
}
