library(mzR)
source("C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/RScripts/MSPairSelector_functions_alt.R")
source("C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/RScripts/DevelopmentScripts/Averagine_functions.R")

data <- openMSfile("C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/RScripts/RMarkDownFiles/20180427_SETD2_001.mzXML", backend="Ramp")
datahead <- header(data)

testMS2 <- 10209
test <- MS1ShiftMatch(data,datahead,testMS2,tolerance=10)

# To test an entire mzXML file:
#####
pairs <- data.frame(MS1=numeric(),Mz1=numeric(),Mz2=numeric())
MS2ScanNums <- datahead[datahead$msLevel==2,"seqNum"]
for(i in 1:length(MS2ScanNums)){
  temp <- MS1ShiftMatch(data,datahead,MS2ScanNums[i],tolerance=10)
  pairs <- rbind(pairs,temp)
}

write.csv(pairs,"C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/RScripts/DevelopmentScripts/MS1PairsTest.csv",row.names=F)

results <- data.frame(MS1=numeric(),MS2=numeric(),Mz1=numeric(),Mz2=numeric(),A.maxIntensity=numeric(),A.RT=numeric(),A.FWHM=numeric(),B.maxIntensity=numeric(),B.RT=numeric(),B.FWHM=numeric(),RT.diff=numeric())
for(i in 1:nrow(pairs)){
  temp <- cbind(pairs[i,],compareRTandFWHM(data,datahead,pairs[i,"MS1"],pairs[i,"Mz1"],pairs[i,"Mz2"]))
  results <- rbind(results,temp)
}

results <- results[order(results$RT.diff,results$FWHM.diff),]

write.csv(results, "C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/RScripts/DevelopmentScripts/MSPairsTestResults.csv", row.names=F)
