#Return a XIC; tolerance=ppm, RTRange is vector of min/max RT, MZFilter is min/max mz values from zoom scan
library(caTools)
library(ggplot2)
library(ggrepel)
source("MSPairSelector_functions_alt.R")

XIC <- function(data, mass, tolerance=100, datahead=NA, RTRange=NA, MZFilter=NA, saveListOfMasses=F, threshold=2, peakArea=T, MS2=T){
  
  if((!is.na(RTRange) & length(RTRange)%%2!=0) |
     (!is.na(MZFilter) & length(MZFilter)%%2!=0)){
    return("For RTRange or MZFilter, min/max values must be provided in a 2 value vector")
  }
  
  if(is.na(datahead)){
    datahead <- header(data)
  }
  
  if(!is.na(RTRange)){
    minScan <- min(datahead[datahead$retentionTime>=min(RTRange),"seqNum"])
    maxScan <- max(datahead[datahead$retentionTime<=max(RTRange),"seqNum"])
    datahead <- datahead[datahead$seqNum>=minScan &
                           datahead$seqNum<=maxScan,]
  }
  
  dataheadMS2 <- datahead[datahead$msLevel==2,]
  datahead <- datahead[datahead$msLevel==1,]
  
  if(!is.na(MZFilter)){
    datahead <- datahead[datahead$lowMZ>=min(MZFilter) &
                           datahead$highMZ<=max(MZFilter),]
  }
  
  intensities <- vector(length=nrow(datahead))
  time <- vector(length=nrow(datahead))
  
  for(i in 1:nrow(datahead)){
    intensity <- peaks(data,datahead$seqNum[i])
    intensity <- intensity[intensity[,1]>mass*(1-(tolerance/1000000)) &
                             intensity[,1]<mass*(1+(tolerance/1000000)) &
                             intensity[,2]>0,,drop=F]
    if(nrow(intensity)>0){
      intensities[i] <- ifelse(nrow(intensity)==1,intensity[,2],trapz(intensity[,1],intensity[,2]))
      time[i] <- datahead[i,"retentionTime"]
    }
    rm(intensity)
  }
  XIC <- data.frame(RT=time,Intensity=intensities)
  XIC[is.na(XIC$Intensity),"Intensity"] <- 0
  XIC <- XIC[XIC$RT!=0,]
  #XIC <- XIC[XIC$Intensity=>0,]
  XIC$RT <- XIC$RT/60
  
  if(peakArea){
    peaks <- peakpicker(XIC,threshold=max(XIC$Intensity)*(threshold/100))
    
    integral <- peaks
    integral$Count <- peaks[,2]-peaks[,1]
    integral <- subset(integral,Count>0)
    RT <- vector(length=nrow(integral))
    Intensity <- vector(length=nrow(integral))
    integralValues <- vector(length=(nrow(integral)))
    for(i in 1:nrow(integral)){
      temp <- XIC[integral[i,1]:integral[i,2],]
      row <- which.max(temp[,"Intensity"])
      RT[i] <- temp[row,"RT"]
      Intensity[i] <- temp[row,"Intensity"]
      integralValues[i] <- trapz(temp$RT,temp$Intensity)
    }
    integral$RT <- RT
    integral$Intensity <- Intensity/max(Intensity)*100
    integral$Values <- integralValues 
  }
  
  XIC$Intensity <- XIC$Intensity/max(XIC$Intensity)*100
  
  if(saveListOfMasses){
    values <<-XIC
  }
  
  if(peakArea){
    polygons <- data.frame(RT=vector(),Intensity=vector(),PeakNumber=vector())
    for(i in 1:nrow(peaks)){
      temp <- XIC[peaks[i,1]:peaks[i,2],]
      temp$PeakNumber=LETTERS[i]
      polygons <- rbind(polygons,temp)
    } 
  }
  
  dataheadMS2 <- dataheadMS2[dataheadMS2$precursorMZ>=mass*(1-(tolerance/1000000)) &
                               dataheadMS2$precursorMZ<=mass*(1+(tolerance/1000000)),]
  
  
  plot <- ggplot(XIC, aes(x=RT, y=Intensity)) +
    geom_path() +
    annotate(geom="text", label=paste0("MZ=",mass, "; ppm = ", tolerance), x=-Inf, y=Inf, hjust=-0.01,vjust=1)
  
  if(peakArea){
    plot <- plot + geom_polygon(data=polygons, aes(x=RT,y=Intensity,fill=PeakNumber, alpha=0.75), show.legend=F) +
             geom_text_repel(data=integral, aes(x=RT, y=Intensity, label=round(Values,3)), nudge_x=0.1,nudge_y=1)
  } 
  
  if(MS2){
    plot <- plot + geom_segment(data=dataheadMS2, aes(y=-Inf,yend=80,x=retentionTime/60, xend=retentionTime/60),colour="red",linetype="dashed") +
      geom_text_repel(data=dataheadMS2, aes(label=seqNum, x=retentionTime/60, y=80), size=4, colour="red")
  }
  
  return(plot)
}
