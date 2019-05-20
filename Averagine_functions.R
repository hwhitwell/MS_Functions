library(enviPat)
library(pracma)
data(isotopes)
proton <- 1.007276467

averagineElements <- list(C=12, H=1.0078250321, O=15.9949146196, N=14.0030740049, S=31.9720710015)

#Used to generate averagineModel showing only the C13 isotope distributions for a mass calculated from a given Mz and charge.
averagineModel <- function(precursorMz, precursorCharge, tolerance=10, plot=T){
  AveragineMass <- averagineElements$C*4.9384 + averagineElements$H* 7.7583 + averagineElements$O*1.4773 + averagineElements$N*1.3577 + averagineElements$S*0.0417
  
  NumOfAveragines <- (precursorMz*precursorCharge-proton*precursorCharge)/AveragineMass
  
  C <- round(NumOfAveragines * 4.9384,0)
  H <- round(NumOfAveragines * 7.7583,0)
  O <- round(NumOfAveragines * 1.4773,0)
  N <- round(NumOfAveragines * 1.3577,0)
  S <- round(NumOfAveragines * 0.0417,0)
  
  Mw <- averagineElements$C*C + averagineElements$H*H + averagineElements$O*O + averagineElements$N*N + averagineElements$S*S
  
  formula <- ifelse(S==0, paste0("C",C,"H",H,"O",O,"N",N), paste0("C",C,"H",H,"O",O,"N",N,"S",S))
  
  #Using the numbers of averagines and known isotopic distributions, calculate the averagine spectra for the precursor mass.
  #Discard the peaks for other isotopes, then return the top 5 most intense (0 C13 to 4 C13)
  averaginemodel <- as.data.frame(isopattern(isotopes, formula)[[1]])
  #Group masses by matching those to the nearest 0.1Da and summing the abundances.
  #Remove any abundances that are less than the 5C13 containing isotope peak
  carbs <- ifelse(max(averaginemodel$`13C`)>=5,5,max(averaginemodel$`13C`))
  cutoff <- ifelse(S==0,
                   averaginemodel[averaginemodel$`1H`==H & averaginemodel$`12C`==(C-carbs) &
                                    averaginemodel$`14N`==N & averaginemodel$`16O`==O,"abundance"],
                   averaginemodel[averaginemodel$`1H`==H & averaginemodel$`12C`==(C-5) &
                                    averaginemodel$`14N`==N & averaginemodel$`16O`==O &
                                    averaginemodel$`32S`==S,"abundance"])
  colnames(averaginemodel)[1] <- "Mz"
  averaginemodel <- averaginemodel %>%
    mutate(Mz=round(Mz,1)) %>%
    mutate(Mz=factor(Mz)) %>%
    group_by(Mz) %>%
    summarise(abundance=sum(abundance)) %>%
    ungroup() %>%
    filter(abundance>cutoff) %>%
    mutate(Mz=as.numeric(as.character(Mz)))
  
  # averaginemodel <- averaginemodel[order(averaginemodel$abundance,decreasing=T),]
  # averaginemodel <- averaginemodel[averaginemodel$`15N`==0,]
  # averaginemodel <- averaginemodel[averaginemodel$`18O`==0,]
  # averaginemodel <- averaginemodel[averaginemodel$`17O`==0,]
  # averaginemodel <- averaginemodel[averaginemodel$`2H`==0,]
  # if(S!=0){
  #   averaginemodel <- averaginemodel[averaginemodel$`34S`==0,]
  #   averaginemodel <- averaginemodel[averaginemodel$`33S`==0,]
  #   averaginemodel <- averaginemodel[averaginemodel$`36S`==0,]
  # }
  averaginemodel <- as.data.frame(averaginemodel)
  averaginemodel$abundance <- averaginemodel$abundance/max(averaginemodel$abundance)*100
  averaginemodel <- averaginemodel[1:5,]
  
  #Optinal plot of the distribution
  if(plot){
    plot <- ggplot(averaginemodel, aes(x=1:nrow(averaginemodel), ymin=0, ymax=abundance)) +
      geom_linerange() +
      ggtitle("Averagine Model Isotopic Distribution")
    print(plot)
  }
  
  #Return the distributions.
  return(averaginemodel)
}

#Calculates a score for all possible combinations of peaks when compared to an averagine model - used to predict the correct precursorMz
determinePrecursor <- function(data, datahead=F, scan, tolerance=10, plot=T, precursorCharge=NA, precursorMz=NA, precursorScan=F){
  
  #Get data header if not provided by user
  if(!is.data.frame(datahead)){
    datahead <- header(data)
  }
  
  #Get peptide mz, charge and precrusors scan number from header
  dataheadRow <- which(datahead$acquisitionNum==scan)
  precursorMz <- ifelse(is.na(precursorMz),datahead[dataheadRow,"precursorMZ"],precursorMz)
  precursorCharge <- ifelse(is.na(precursorCharge),datahead[dataheadRow, "precursorCharge"],precursorCharge)
  if(precursorCharge==0) return(NULL)
  precursorScan <- ifelse(precursorScan,scan,datahead[dataheadRow, "precursorScanNum"])
  
  #Peak list for the MS spectra
  peakList <- mzR::peaks(data, which(datahead$acquisitionNum==precursorScan))
  
  #Isolate the peaks that are within 3.5Da of the precursor
  peakList <- peakList[which(peakList[,1]<(precursorMz+3.5/precursorCharge) & peakList[,1]>(precursorMz-3.5/precursorCharge)),]
  
  #
  if(nrow(peakList)==0){
    return(NA)
  }
  
  #Normalise peak list intensities
  peakList[,2] <- peakList[,2]/max(peakList[,2])*100
  
  #Pick peaks
  #cut off is half the averagine intensity for the 4C13 isotopic peak
  averaginePeptide <- averagineModel(precursorMz, precursorCharge, tolerance, F)
  pickedPeaks <- peakpicker(peakList, averaginePeptide[5,2])
  
  #Determine the averageMz for each peak
  averageMz <- matrix(NA,ncol=2,nrow=nrow(pickedPeaks))
  for(i in 1:nrow(averageMz)){
    averageMz[i,2] <- max(peakList[pickedPeaks[i,1]:pickedPeaks[i,2],2])
    averageMz[i,1] <- weighted.mean(peakList[pickedPeaks[i,1]:pickedPeaks[i,2],1],
                                    peakList[pickedPeaks[i,1]:pickedPeaks[i,2],2])
  }
  
  averageMz <- as.data.frame(averageMz)
  colnames(averageMz) <- c("Mz", "Intensity")
  
  #Differnce between C12 and C13
  carbonDiff <- (isotopes[which(isotopes$element=="[13]C"),"mass"]-isotopes[which(isotopes$element=="[12]C"),"mass"])
  #Select peaks that may belong to the isotope cluster, calculating them for -3Da and +2Da from precursorMz.
  theoreticalMasses <- seq(precursorMz-3*carbonDiff/precursorCharge,precursorMz+2*carbonDiff/precursorCharge,carbonDiff/precursorCharge)
  matchedMasses <- data.frame(Mz=theoreticalMasses,Intensity=NA)
  for(i in 1:nrow(matchedMasses)){
    #Find the closest match between the expected and the average Mz
    closest <- averageMz[which.min(abs(averageMz$Mz-matchedMasses[i,"Mz"])),]
    #check it is within tolerance of the expected mass, if it is, return the intensity, if it is not, return NA
    matchedMasses[i, "Intensity"] <- ifelse(abs(closest$Mz-matchedMasses[i,"Mz"])/matchedMasses[i,"Mz"]*1000000<=tolerance,
                                            closest$Intensity,NA)
  }
  
  #Test there is enough peaks
  matchedMasses$peakNum <- 1:nrow(matchedMasses)
  matchedMassesNaOmit <- na.omit(matchedMasses)
  if(nrow(matchedMassesNaOmit)<=2){
    return(NA)
  }
  
  #Calculate all possible combinations of peaks
  
  combinations <- combn(matchedMassesNaOmit$peakNum, m=3, simplify = F)  
  if(nrow(matchedMassesNaOmit)>3){
    for(i in 4:nrow(matchedMassesNaOmit)){
      temp <- combn(matchedMassesNaOmit$peakNum, m=i, simplify = F)
      combinations <- c(combinations, temp)
    }
  }
  
  calculateScore <- function(x){
    peaks <- matchedMasses[match(x,matchedMasses$peakNum),]
    
    #Remormalise peak intensities to the precursorMass
    peaks$Intensity <- peaks$Intensity/max(peaks$Intensity)*100
    averaginePeptide <- averagineModel(peaks$Mz[1], precursorCharge, plot=F)
    averaginePeaks <- x-x[1]+1
    if(length(which(averaginePeaks<=5))!=length(averaginePeaks)){
      remove <- which(averaginePeaks>5)
      averaginePeaks <- averaginePeaks[-remove]
      peaks <- peaks[-remove,]
    }
    averaginePeptide <- averaginePeptide[averaginePeaks,]
    
    #The score is the weighted average difference between the averagine peak and the experimental peak (excluding the 100% peak), weighting by abundance
    
    # score <- ifelse(nrow(peaks)<=1,NA,
    #                 sum(abs(peaks$Intensity[-1]-averaginePeptide$abundance[2:nrow(peaks)])*peaks$Intensity[-1])/sum(peaks$Intensity[-1]))
    # 
    
    # The score is the weighted average difference between the averagine peak and the experimental peak, weighting by abundance
    
    score <- ifelse(nrow(peaks)<=2,NA,
                    sum(abs(peaks$Intensity-averaginePeptide$abundance)*peaks$Intensity)/sum(peaks$Intensity))
    return(score)
  }
  
  scores <- lapply(combinations, function(x) calculateScore(x))
  matchedMz <- lapply(combinations, function(x) matchedMasses[which(matchedMasses$peakNum==x[1]),"Mz"])
  results <- data.frame(Mz=unlist(matchedMz),Score=unlist(scores))
  
  #test that not all the result scores are NA
  test <- na.omit(results)
  if(nrow(test)==0){
    return(NA)
  }
  
  
  if(plot){
    precursorMz <- results$Mz[which.min(results$Score)]
    bestScore <- results$Score[which.min(results$Score)]
    bestPeaks <- combinations[[which.min(results$Score)]]
    avraginePeaks <- bestPeaks-bestPeaks[1]+1
    if(length(which(avraginePeaks<=5))!=length(avraginePeaks)){
      remove <- which(avraginePeaks>5)
      avraginePeaks <- avraginePeaks[-remove]
    }
    
    isoPeaks <- matchedMasses[which(matchedMasses$Mz==precursorMz):nrow(matchedMasses),]
    isoPeaks <- isoPeaks[avraginePeaks,]
    isoPeaks$Intensity <- isoPeaks$Intensity/max(isoPeaks$Intensity,na.rm = T)*100
    averaginePeptide <- averagineModel(precursorMz, precursorCharge, tolerance, F)
    temp <- data.frame(Mz=c(avraginePeaks,1:nrow(averaginePeptide)),Intensity=c(isoPeaks$Intensity,averaginePeptide$abundance),Data=c(rep("Exp",nrow(isoPeaks)),rep("Averagine",nrow(averaginePeptide))))
    temp <- ggplot(temp, aes(x=Mz, ymin=0, ymax=Intensity, colour=Data, size=Data)) +
      geom_linerange() +
      ggtitle("Precursor and Averagine Isotopic Distributions", subtitle=paste0("PrecursorMz = ", round(precursorMz,3), " with score of ",round(bestScore,3)," using peaks ", paste(avraginePeaks,collapse=",")))
    print(temp)
  }
  
  return(results)
}

#For a given MS2, look in the precursor scan for a heavy/light pair.
MS1ShiftMatch <- function(data,datahead=F,MS2ScanNum,Delta=4.022185078, tolerance=100){
  #1 Check precursor with determinePrecursor - if the score is <8.6, use this precursor.
  #2 Look for mz +/- the appropriate mass shift and isolate the MS spectra for these peaks.
  #3 Check to see if the identified peak is the first peak in an isotope cluster of the same charge as the original precursor.
  #3 Check the score is the minimum for this peak and check assuming the charge is +/- 1
  #4 If the m/z is correct distance and it is the correct charge/isotope peak then keep!
  
  if(!is.data.frame(datahead)){
    datahead <- header(data)
  }
  
  precursorMz <- datahead[datahead$acquisitionNum==MS2ScanNum,"precursorMZ"]
  precursorCharge <- datahead[datahead$acquisitionNum==MS2ScanNum,"precursorCharge"]
  if(precursorCharge==0) return(NULL)
  precursorScan <- datahead[datahead$acquisitionNum==MS2ScanNum,"precursorScanNum"]
  
  #Check precursor Mz (assume charge is correct).
  scores <- determinePrecursor(data,datahead,MS2ScanNum,plot=F)
  if(!is.na(scores)){
    scores <- na.omit(scores)
    if(nrow(scores)>0){
      scores <- scores[which.min(scores$Score),]
      if(scores$Score<8.6){
        precursorMz <- scores$Mz
      }
    }
  }
  
  
  #Isolate all peaks in precusor scan +/- the mass shift.
  Delta <- Delta/precursorCharge
  peakList <- mzR::peaks(data,which(datahead$acquisitionNum==precursorScan))
  peakList <- peakList[peakList[,1]<=(precursorMz+Delta+4.1) & peakList[,1]>=(precursorMz-Delta-4.1),]
  if(nrow(peakList)==0) return(NULL)
  
  #Normalise peak list
  peakList[,2] <- peakList[,2]/max(peakList[,2])*100
  
  #Identify peaks and determine average mass of each peak and adjust the precursorMz to the average 
  pickedPeaks <- peakpicker(peakList,threshold = 5)
  averageMz <- matrix(NA,ncol=2,nrow=nrow(pickedPeaks))
  for(i in 1:nrow(averageMz)){
    averageMz[i,2] <- max(peakList[pickedPeaks[i,1]:pickedPeaks[i,2],2])
    averageMz[i,1] <- weighted.mean(peakList[pickedPeaks[i,1]:pickedPeaks[i,2],1],
                                    peakList[pickedPeaks[i,1]:pickedPeaks[i,2],2])
  }
  
  averageMz <- as.data.frame(averageMz)
  colnames(averageMz) <- c("Mz", "Intensity")
  
  precursorMz <- averageMz[which.min(abs(averageMz$Mz-precursorMz)),"Mz"]
  
  #Look if there is a peak +/- the appropriate delta
  above <- abs(precursorMz+Delta-averageMz$Mz)
  above <- averageMz[above<=(tolerance*precursorMz/1000000),]
  above <- above[which.max(above$Intensity),"Mz"]
  
  below <- abs(precursorMz-Delta-averageMz$Mz)
  below <- averageMz[below<=(tolerance*precursorMz/1000000),]
  below <- below[which.max(below$Intensity),"Mz"]
  
  #If there is a peak with the appropriate mass shift, check it is the beginning of an isotope cluster. If it is - it's a heavy/light pair!

  abovePair <- F
  if(length(above)>0){
    abovePrecursor <- determinePrecursor(data,datahead,scan=precursorScan,tolerance=10,plot=F,precursorCharge=precursorCharge,precursorMz=above,precursorScan=T)
    if(is.data.frame(abovePrecursor)){
      abovePrecursor <- abovePrecursor[which.min(abovePrecursor$Score),]
      } else {
        abovePrecursor <- FALSE
      }
    if(abovePrecursor!=F){
      abovePair <- ifelse(abovePrecursor$Mz==above & abovePrecursor$Score<8.6,T,F)
    }
  }
  
  belowPair <- F
  if(length(below)>0){
    belowPrecursor <- determinePrecursor(data,datahead,scan=precursorScan,tolerance=10,plot=F,precursorCharge=precursorCharge,precursorMz=below,precursorScan=T)
    if(is.data.frame(belowPrecursor)){
      belowPrecursor <- belowPrecursor[which.min(belowPrecursor$Score),]
    } else {
      belowPrecursor <- FALSE
    }
    if(belowPrecursor!=F){
      belowPair <- ifelse(belowPrecursor$Mz==below & belowPrecursor$Score<8.6,T,F)
    }
  }
  
  returnResults <- F
  
  if(abovePair & belowPair){
    results <- data.frame(MS1=precursorScan,MS2=MS2ScanNum,Mz1=precursorMz,Mz2=c(above,below))
    returnResults <- T
  }
  
  if(abovePair & !belowPair){
    results <- data.frame(MS1=precursorScan,MS2=MS2ScanNum,Mz1=precursorMz,Mz2=above)
    returnResults <- T
  }
  
  if(belowPair & !abovePair){
    results <- data.frame(MS1=precursorScan,MS2=MS2ScanNum,Mz1=precursorMz,Mz2=below)
    returnResults <- T
  }
  
  if(returnResults){
    return(results)
  }
  
}

#Calculate the retention time (at maximum intensity), maximum intensity and FWHM for a given precursor mass.
RTandFWHM <- function(data,datahead,MS1ScanNum,Mz,tolerance=10){
  MS1Scans <- datahead[datahead$msLevel==1,"seqNum"]
  tolerance <- tolerance/1000000
  
  #Calculate the average mz and thus the intensity (area under the peak of all matching peaks)
  peakList <- mzR::peaks(data,which(datahead$seqNum==MS1ScanNum))
  peakList <- peakList[(peakList[,1]>=(Mz-Mz*tolerance) & peakList[,1]<=(Mz+Mz*tolerance)),]
  if(nrow(peakList)==1){
    intensity <- peakList[,2]
  } else {
    intensity <- trapz(peakList[,1],peakList[,2])
  }
  
  intensityVector <- intensity
  rtVector <- datahead[datahead$seqNum==MS1ScanNum,"retentionTime"]
  
  #Calculate the intensity for the following scan, if that is not upward, then go the otherway.
  scanIndex <- which(MS1Scans==MS1ScanNum)
  currentScan <- MS1Scans[scanIndex+1]
  
  peakList <- mzR::peaks(data,which(datahead$seqNum==currentScan))
  peakList <- peakList[(peakList[,1]>=(Mz-Mz*tolerance) & peakList[,1]<=(Mz+Mz*tolerance)),]
  if(nrow(peakList)==1){
    intensity <- peakList[,2]
  } else {
    intensity <- trapz(peakList[,1],peakList[,2])
  }
  
  intensityVector <- c(intensityVector,intensity)
  rtVector <- c(rtVector,datahead[datahead$seqNum==currentScan,"retentionTime"])

  direction <- ifelse(intensityVector[2]>intensityVector[1],1,-1)
  
  #Repeat until the calculated intensity is less than 50% of the maximum intensity
    step <- ifelse(direction==1,2,-1)
  while(intensity>50){
    currentScan <- MS1Scans[scanIndex+step]
    
    peakList <- mzR::peaks(data,which(datahead$seqNum==currentScan))
    peakList <- peakList[(peakList[,1]>=(Mz-Mz*tolerance) & peakList[,1]<=(Mz+Mz*tolerance)),]
    if(nrow(peakList)==1){
      intensity <- peakList[,2]
    } else {
      intensity <- trapz(peakList[,1],peakList[,2])
    }
    
    intensityVector <- c(intensityVector,intensity)
    rtVector <- c(rtVector,datahead[datahead$seqNum==currentScan,"retentionTime"])
    
    intensity <- intensity/max(intensityVector)*100
    step <- step + direction
  }
    
    #Repeat in the other direction
    direction <- direction*-1
    step <- ifelse(direction==1,2,-1)
    intensity <- 100
    while(intensity>50){
      currentScan <- MS1Scans[scanIndex+step]
      
      peakList <- mzR::peaks(data,which(datahead$seqNum==currentScan))
      peakList <- peakList[(peakList[,1]>=(Mz-Mz*tolerance) & peakList[,1]<=(Mz+Mz*tolerance)),]
      if(nrow(peakList)==1){
        intensity <- peakList[,2]
      } else {
        intensity <- trapz(peakList[,1],peakList[,2])
      }
      
      intensityVector <- c(intensityVector,intensity)
      rtVector <- c(rtVector,datahead[datahead$seqNum==currentScan,"retentionTime"])
      
      intensity <- intensity/max(intensityVector)*100
      step <- step + direction
    }
    
    #order by retention time (seconds)
    intensityVector <- intensityVector[order(rtVector)]
    rtVector <- rtVector[order(rtVector)]
    
    #check the ends of intensityVector only contain 1 value below 50% at either end of the distribution
    normalisedIntensityVector <- intensityVector/max(intensityVector)*100
    normalisedIntensityVector <- normalisedIntensityVector[(min(which(normalisedIntensityVector>50))-1):(max(which(normalisedIntensityVector>50))+1)]
    
    #Calculate the retnetion time at the lower and upper 50%
    lower <- lm(rtVector[1:2]~normalisedIntensityVector[1:2])
    lower <- as.numeric(lower$coefficients[1] + lower$coefficients[2]*50)
    
    upper <- lm(rtVector[(length(rtVector)-1):length(rtVector)]~normalisedIntensityVector[(length(normalisedIntensityVector)-1):length(normalisedIntensityVector)])
    upper <- as.numeric(upper$coefficients[1] + upper$coefficients[2]*50)
    
    #Return the maximum intensity and the FwHM
    results <- data.frame(maxIntensity=max(intensityVector),RT=rtVector[which.max(intensityVector)],FwHM=upper-lower)
    return(results)
}

#Compare a pair of Mz values.
compareRTandFWHM <- function(data,datahead,MS1ScanNum,Mz1,Mz2,tolerance=10){
  first <- RTandFWHM(data,datahead,MS1ScanNum,Mz1,tolerance)
  second <- RTandFWHM(data,datahead,MS1ScanNum,Mz2,tolerance)
  
  colnames(first) <- paste0("A",".",colnames(first))
  colnames(second) <- paste0("B",".",colnames(second))
  
  results <- cbind(first,second)
  results$FWHM.diff <- abs(results$A.FwHM-results$B.FwHM)
  results$RT.diff <- abs(results$A.RT-results$B.RT)
  return(results)
}
