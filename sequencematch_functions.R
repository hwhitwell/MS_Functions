library(mzR)
library(ggplot2)
library(OrgMassSpecR)

PeptideSpectrum <- function (expt, theory, t = 0.4, b = 5, label = "", xlim = c(100, 
                                                             1500), supress = FALSE, binwidth=100,numperbin=10) 
{
  if (!is.data.frame(expt)) 
    stop("the expt argument must be a data frame")
  if (t < 0 | t > 100) 
    stop("the t argument must be between 0 and 100")
  if (b < 0 | b > 100) 
    stop("the b argument must be between 0 and 100")
  if (!is.character(label)) 
    stop("the label argument must be a character string")
  expt <- data.frame(mz = expt[, 1], intensity = expt[, 2])
  expt$normalized <- round((expt$intensity/max(expt$intensity)) * 
                             100)
  
  expt <- expt[order(expt[,2],decreasing=T),]
  x <- min(expt[,1])
  while(x<max(expt[,1])){
    range <- which(expt[,1]>x & expt[,1]<(x+binwidth))
    if(length(range)!=0){
      if(length(range)>numperbin){
        temp <- expt[range,][1:numperbin,]
      } else {
        temp <- expt[range,]}
      expt <- expt[-range,]
      expt <- rbind(expt,temp)
    }
    x <- x+binwidth
  }
  
  expt <- subset(expt, expt$normalized >= b)
  matches <- vector("list")
  for (i in 1:nrow(expt)) {
    tmp_matches <- data.frame(NULL)
    tmp_matches <- theory[expt$mz[i] >= theory$ms2mz - t & 
                            expt$mz[i] <= theory$ms2mz + t, ]
    num_tmp_matches <- nrow(tmp_matches)
    expt_mz <- rep(expt$mz[i], times = num_tmp_matches)
    expt_intensity <- rep(expt$normalized[i], times = num_tmp_matches)
    matches[[i]] <- data.frame(expt_mz, expt_intensity, tmp_matches)
  }
  identifications <- as.data.frame(do.call("rbind", matches))
  if (nrow(identifications) == 0) 
    stop("no peaks were identified")
  row.names(identifications) <- 1:nrow(identifications)
  identifications$error <- as.character(format(round(identifications$expt_mz - 
                                   identifications$ms2mz, digits = 4),scientific = F))
    identifications$error[!identifications$error<0] <- paste0("+",identifications$error[!identifications$error<0])
  num_identifications <- nrow(identifications)
  getLocation <- function(type) {
    tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
    tmp <- strsplit(tmp, split = "[[:alpha:]]")[[1]][2]
    return(as.numeric(tmp))
  }
  location <- sapply(identifications$ms2type, getLocation)
  getSeries <- function(type) {
    tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
    tmp <- strsplit(tmp, split = "[[:digit:]]")[[1]][1]
    return(tmp)
  }
  series <- sapply(identifications$ms2type, getSeries)
  color <- sapply(series, function(x) ifelse(x =="p"|x=="w"|x=="n","green", ifelse(x=="b" | x == 
                                               "c", "red", "blue")))
  plot(expt$mz, expt$normalized, type = "h", xlim, ylim = c(0, 
                                                            150), xlab = "m/z", ylab = "intensity (%)", yaxs = "i", 
       yaxt = "n")
  axis(2, at = c(0, 50, 100), labels = c(0, 50, 100))
  if (supress == FALSE) {
    x_range <- xlim[2] - xlim[1]
    y_position <- vector("numeric")
    if (num_identifications == 1) {
      y_position[i] <- identifications$expt_intensity[i]
    }
    else {
      for (i in 1:(num_identifications - 1)) {
        if ((identifications$expt_mz[i + 1] - identifications$expt_mz[i])/x_range < 
            0.025 & all.equal(identifications$expt_intensity[i], 
                              100) != TRUE) {
          y_position[i] <- identifications$expt_intensity[i] + 
            40
          lines(rep(identifications$expt_mz[i], 2), c(identifications$expt_intensity[i], 
                                                      identifications$expt_intensity[i] + 30), 
                lty = "dotted", col = color[i])
        }
        else y_position[i] <- identifications$expt_intensity[i]
      }
    }
    y_position[num_identifications] <- identifications$expt_intensity[num_identifications]
    text(identifications$expt_mz, y_position + 12, labels = paste0(identifications$ms2type," ",round(identifications$ms2mz,3)," ",identifications$error), 
         col = color, srt = 90, cex = 0.75, family="sans")
    lines(identifications$expt_mz, identifications$expt_intensity, col=color, type="h", lwd=2)
  }
  seq_vector <- strsplit(as.character(identifications$ms1seq)[1], 
                         split = "")[[1]]
  num_residues <- length(seq_vector)
  plot.window(xlim = c(1, 20), ylim = c(0, 10))
  text(c(1:num_residues), 9, labels = seq_vector)
  for (i in 1:length(series)) {
    if (series[i] == "b" | series[i] == "c") 
      lines(c(location[i] + 0.25, location[i] + 0.55), 
            c(8.5, 9.5), col = "red")
    else lines(c(num_residues - location[i] + 0.45, num_residues - 
                   location[i] + 0.75), c(8.5, 9.5), col = "blue")
  }
  text(1, 10, label, pos=4, cex=0.8,col="red")
  return(identifications)
}

PeptideSpectrumNoPlot <- function (expt, theory, t = 0.4, b = 5, label = "", xlim = c(100, 
                                                                                     1500), supress = FALSE, binwidth=100, numperbin=10) 
{
  if (!is.data.frame(expt)) 
    stop("the expt argument must be a data frame")
  if (t < 0 | t > 100) 
    stop("the t argument must be between 0 and 100")
  if (b < 0 | b > 100) 
    stop("the b argument must be between 0 and 100")
  if (!is.character(label)) 
    stop("the label argument must be a character string")
  expt <- data.frame(mz = expt[, 1], intensity = expt[, 2])
  
  expt <- expt[order(expt[,2],decreasing=T),]
  x <- min(expt[,1])
  while(x<max(expt[,1])){
    range <- which(expt[,1]>x & expt[,1]<(x+binwidth))
    if(length(range)!=0){
      if(length(range)>numperbin){
        temp <- expt[range,][1:numperbin,]
      } else {
        temp <- expt[range,]}
      expt <- expt[-range,]
      expt <- rbind(expt,temp)}
    x <- x+binwidth
  }
  
  expt$normalized <- round((expt$intensity/max(expt$intensity)) * 
                             100)
  expt <- subset(expt, expt$normalized >= b)
  matches <- vector("list")
  for (i in 1:nrow(expt)) {
    tmp_matches <- data.frame(NULL)
    tmp_matches <- theory[expt$mz[i] >= theory$ms2mz - t & 
                            expt$mz[i] <= theory$ms2mz + t, ]
    num_tmp_matches <- nrow(tmp_matches)
    expt_mz <- rep(expt$mz[i], times = num_tmp_matches)
    expt_intensity <- rep(expt$normalized[i], times = num_tmp_matches)
    matches[[i]] <- data.frame(expt_mz, expt_intensity, tmp_matches)
  }
  identifications <- as.data.frame(do.call("rbind", matches))
  if (nrow(identifications) == 0) 
    return("no peaks were identified")
  row.names(identifications) <- 1:nrow(identifications)
  identifications$error <- round(identifications$expt_mz - 
                                   identifications$ms2mz, digits = 3)
  num_identifications <- nrow(identifications)
  getLocation <- function(type) {
    tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
    tmp <- strsplit(tmp, split = "[[:alpha:]]")[[1]][2]
    return(as.numeric(tmp))
  }
  location <- sapply(identifications$ms2type, getLocation)
  getSeries <- function(type) {
    tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
    tmp <- strsplit(tmp, split = "[[:digit:]]")[[1]][1]
    return(tmp)
  }
  return(identifications)
}

FragmentPeptideAmended <- function (sequence, fragments = "by", IAA = TRUE, N15 = FALSE, 
          custom = list()) 
{
  results_list <- vector("list")
  for (sequence_number in 1:length(sequence)) {
    peptide_vector <- strsplit(sequence[sequence_number], 
                               split = "")[[1]]
    peptide_length <- length(peptide_vector)
    if (peptide_length < 2) 
      stop("sequence must contain two or more residues")
    C <- 12
    H <- 1.0078250321
    O <- 15.9949146221
    S <- 31.97207069
    N <- ifelse(N15 == TRUE, 15.0001088984, 14.0030740052)
    proton <- 1.007276466
    electron <- 0.00054857990943
    residueMass <- function(residue) {
      if (residue == "A") 
        mass = C * 3 + H * 5 + N + O
      if (residue == "R") 
        mass = C * 6 + H * 12 + N * 4 + O
      if (residue == "N") 
        mass = C * 4 + H * 6 + N * 2 + O * 2
      if (residue == "D") 
        mass = C * 4 + H * 5 + N + O * 3
      if (residue == "E") 
        mass = C * 5 + H * 7 + N + O * 3
      if (residue == "Q") 
        mass = C * 5 + H * 8 + N * 2 + O * 2
      if (residue == "G") 
        mass = C * 2 + H * 3 + N + O
      if (residue == "H") 
        mass = C * 6 + H * 7 + N * 3 + O
      if (residue == "I") 
        mass = C * 6 + H * 11 + N + O
      if (residue == "L") 
        mass = C * 6 + H * 11 + N + O
      if (residue == "K") 
        mass = C * 6 + H * 12 + N * 2 + O
      if (residue == "M") 
        mass = C * 5 + H * 9 + N + O + S
      if (residue == "F") 
        mass = C * 9 + H * 9 + N + O
      if (residue == "P") 
        mass = C * 5 + H * 7 + N + O
      if (residue == "S") 
        mass = C * 3 + H * 5 + N + O * 2
      if (residue == "T") 
        mass = C * 4 + H * 7 + N + O * 2
      if (residue == "W") 
        mass = C * 11 + H * 10 + N * 2 + O
      if (residue == "Y") 
        mass = C * 9 + H * 9 + N + O * 2
      if (residue == "V") 
        mass = C * 5 + H * 9 + N + O
      if (residue == "C" & IAA == FALSE) 
        mass = C * 3 + H * 5 + N + O + S
      if (residue == "C" & IAA == TRUE) 
        mass <- ifelse(N15 == FALSE, C * 5 + H * 8 + 
                         N * 2 + O * 2 + S, C * 5 + H * 8 + N + 14.0030740052 + 
                         O * 2 + S)
      if (length(custom) != 0) 
        for (i in 1:length(custom$code)) if (residue == 
                                             custom$code[i]) 
          mass = custom$mass[i]
        return(mass)
    }
    masses <- sapply(peptide_vector, residueMass)
    pm <- sum(masses)
    p1 <- round(pm + H * 2 + O + proton, digits = 5)
    p2 <- round((pm + H * 2 + O + (2 * proton))/2, digits = 5)
    p3 <- round((pm + H * 2 + O + (3 * proton))/3, digits = 5)
    p4 <- round((pm + H * 2 + O + (4 * proton))/4, digits = 5)
    if (fragments == "by") {
      b1 <- vector(mode = "numeric", length = 0)
      b2 <- vector(mode = "numeric", length = 0)
      b3 <- vector(mode = "numeric", length = 0)
      bs <- vector(mode = "character", length = 0)
      bi <- vector(mode = "integer", length = 0)
      y1 <- vector(mode = "numeric", length = 0)
      y2 <- vector(mode = "numeric", length = 0)
      y3 <- vector(mode = "numeric", length = 0)
      ys <- vector(mode = "character", length = 0)
      yi <- vector(mode = "integer", length = 0)
      for (i in 1:(peptide_length - 1)) {
        mass <- sum(masses[1:i])
        b1[i] <- round(mass + proton, digits = 5)
        b2[i] <- round((b1[i] + proton)/2, digits = 5)
        b3[i] <- round((b1[i] + 2*proton)/3, digits = 5)
        bs[i] <- paste(peptide_vector[1:i], collapse = "")
        bi[i] <- i
      }
      for (j in 2:peptide_length) {
        mass <- sum(masses[j:peptide_length])
        y1[j - 1] <- round(mass + H * 2 + O + proton, 
                           digits = 5)
        y2[j - 1] <- round((y1[j - 1] + proton)/2, digits = 5)
        y3[j - 1] <- round((y1[j - 1] + 2*proton)/3, digits = 5)
        ys[j - 1] <- paste(peptide_vector[j:peptide_length], 
                           collapse = "")
        yi[j - 1] <- peptide_length - j + 1
      }
      ms1seq <- rep(sequence[sequence_number], times = ((3 * 
                                                           (length(bi))) + (3 * (length(yi)))))
      ms1z1 <- rep(p1, times = ((3 * (length(bi))) + (3 * 
                                                        (length(yi)))))
      ms1z2 <- rep(p2, times = ((3 * (length(bi))) + (3 * 
                                                        (length(yi)))))
      ms1z3 <- rep(p3, times = ((3 * (length(bi))) + (3 * 
                                                        (length(yi)))))
      ms2seq <- c(rep(bs, times = 3), rep(ys, times = 3))
      b1.type <- paste("[b", bi, "]1+", sep = "")
      b2.type <- paste("[b", bi, "]2+", sep = "")
      b3.type <- paste("[b", bi, "]3+", sep = "")
      y1.type <- paste("[y", yi, "]1+", sep = "")
      y2.type <- paste("[y", yi, "]2+", sep = "")
      y3.type <- paste("[y", yi, "]3+", sep = "")
      ms2type <- c(b1.type, b2.type, b3.type, y1.type, y2.type, y3.type)
      ms2mz <- c(b1, b2, b3, y1, y2, y3)
    }
    if (fragments == "cz") {
      c1 <- vector(mode = "numeric", length = 0)
      c2 <- vector(mode = "numeric", length = 0)
      cs <- vector(mode = "character", length = 0)
      ci <- vector(mode = "integer", length = 0)
      z1 <- vector(mode = "numeric", length = 0)
      z2 <- vector(mode = "numeric", length = 0)
      zs <- vector(mode = "character", length = 0)
      zi <- vector(mode = "integer", length = 0)
      for (i in 1:(peptide_length - 1)) {
        mass <- sum(masses[1:i])
        c1[i] <- round(mass + 3 * H + N + proton, digits = 5)
        c2[i] <- round((c1[i] + proton)/2, digits = 5)
        cs[i] <- paste(peptide_vector[1:i], collapse = "")
        ci[i] <- i
      }
      for (j in 2:peptide_length) {
        mass <- sum(masses[j:peptide_length])
        z1[j - 1] <- round(mass + O - N, digits = 5)
        z2[j - 1] <- round((z1[j - 1] + proton)/2, digits = 5)
        zs[j - 1] <- paste(peptide_vector[j:peptide_length], 
                           collapse = "")
        zi[j - 1] <- peptide_length - j + 1
      }
      ms1seq <- rep(sequence[sequence_number], times = ((2 * 
                                                           (length(ci))) + (2 * (length(zi)))))
      ms1z1 <- rep(p1, times = ((2 * (length(ci))) + (2 * 
                                                        (length(zi)))))
      ms1z2 <- rep(p2, times = ((2 * (length(ci))) + (2 * 
                                                        (length(zi)))))
      ms1z3 <- rep(p3, times = ((2 * (length(ci))) + (2 * 
                                                        (length(zi)))))
      ms2seq <- c(rep(cs, times = 2), rep(zs, times = 2))
      c1.type <- paste("[c", ci, "]1+", sep = "")
      c2.type <- paste("[c", ci, "]2+", sep = "")
      z1.type <- paste("[z", zi, "]1+", sep = "")
      z2.type <- paste("[z", zi, "]2+", sep = "")
      ms2type <- c(c1.type, c2.type, z1.type, z2.type)
      ms2mz <- c(c1, c2, z1, z2)
    }
    water <- MonoisotopicMass(ListFormula("H2O"))
    ammonia <- MonoisotopicMass(ListFormula("NH2"))
    results_list[[sequence_number]] <- data.frame(ms1seq, 
                                                  ms1z1, ms1z2, ms1z3, ms2seq, ms2type, ms2mz)
    results_list[[sequence_number]] <- rbind(results_list[[sequence_number]],
                                             data.frame(ms1seq=rep(ms1seq[1],12),
                                                        ms1z1=rep(ms1z1[1],12),
                                                        ms1z2=rep(ms1z2[1],12),
                                                        ms1z3=rep(ms1z3[1],12),
                                                        ms2seq=rep(ms1seq[1],12),
                                                        ms2type=c("[p0]1+","[p0]2+","[p0]3+","[p0]4+",
                                                                  "[w0]1+","[w0]2+","[w0]3+","[w0]4+",
                                                                  "[n0]1+","[n0]2+","[n0]3+","[n0]4+"),
                                                        ms2mz=c(ms1z1[1],ms1z2[1],ms1z3[1],round((ms1z1[1]+3*proton)/4,5),
                                                                ms1z1[1]-water,ms1z2[1]-water/2,ms1z3[1]-water/3,round((ms1z1[1]+3*proton)/4-water/4,5),
                                                                ms1z1[1]-ammonia,ms1z2[1]-ammonia/2,ms1z3[1]-ammonia/3,round((ms1z1[1]+3*proton)/4-ammonia/4,5))))
  }
  return(as.data.frame(do.call("rbind", results_list)))
}


SpectrumFinder <- function(mzxml, minrt=0, maxrt=7200, peptide, tolerance=0.5, base=2, mzrange, modifiedaa=NA, aminoacidmass=NA, bestspectrum=NA, IAA=F, N15=F,deltacutoff=100,precursorCharge=NA, binwidth=100,numperbin=10){
  aacode <- c("R","H","K","D","E","S","T","N","Q","C","U","G","P","A","I","L","M","F","W","Y","V")
  if(is.na(modifiedaa)==F){aacode <- c(aacode, modifiedaa)}
  checkpeptide <- unlist(strsplit(peptide,""))
  if (TRUE %in% is.na(match(checkpeptide,aacode))){
    stop("Unspecified amino acid in peptide")
  }
  if(is.na(modifiedaa) != is.na(aminoacidmass)){
    stop("Modified amino acids specified without masses or vice-versa")
  }
  if(class(mzxml)[1]!="mzRramp"){
    stop("Data file need to be of class mzRamp. Read in mzXML using mzR::openMSfile()")
  }
  if(is.na(modifiedaa)==F){
    pepfragments <- FragmentPeptideAmended(peptide, "by", IAA,  N15, custom=list(code=modifiedaa,mass=aminoacidmass))
  } else
    pepfragments <- FragmentPeptideAmended(peptide, "by", IAA, N15)
  
  mzxmltemp <- header(mzxml)
  mzxmlhead <- subset(mzxmltemp, msLevel==2)
  if(is.na(precursorCharge)==T){
    mzxmlhead <- mzxmlhead[which(abs((mzxmlhead$precursorMZ*mzxmlhead$precursorCharge-(mzxmlhead$precursorCharge-1)*1.007276)-pepfragments[1,2])<=deltacutoff),]
  } else {
    mzxmlhead <- mzxmlhead[which(abs((mzxmlhead$precursorMZ*precursorCharge-(precursorCharge-1)*1.007276)-pepfragments[1,2])<=deltacutoff),]
  }
  ms2scans <- mzxmlhead$seqNum[which(mzxmlhead$retentionTime>=minrt & mzxmlhead$retentionTime<=maxrt)]
  
  scores <- vector()
  scan <- vector()
  delta <- vector()
  
  if(is.na(bestspectrum)==F){
    print("Pass 1")
    bestspectrum <- PeptideSpectrumNoPlot(as.data.frame(peaks(mzxml,which(mzxmltemp$seqNum==bestspectrum))), pepfragments, tolerance, base, xlim=mzrange,binwidth=binwidth,numperbin=numperbin)
    limit <- length(unique(bestspectrum$ms2type))
    index <- 0
    for(i in ms2scans){
      msmass <- as.data.frame(peaks(mzxml,which(mzxmltemp$seqNum==i)))
      msmass <- msmass[order(msmass[,2],decreasing = T),]
      if(nrow(msmass)>0){
        temp <- PeptideSpectrumNoPlot(msmass,pepfragments, tolerance, base, xlim=mzrange, binwidth=binwidth, numperbin=numperbin)
        if(temp!="no peaks were identified"){
          if(length(unique(temp$ms2type)) >= limit*0.5){
            index <- index + 1
            score <- round(sum(msmass[match(temp$expt_mz,msmass[,1]),2])/msmass[1,2]*100,0)
            scores <- c(scores, score)
            scan <- c(scan, i)
            if(is.na(precursorCharge)==T){
              deltatemp <- round(pepfragments[1,2]-(mzxmltemp[which(mzxmltemp$seqNum==i),"precursorMZ"]*mzxmltemp[which(mzxmltemp$seqNum==i),"precursorCharge"])+(mzxmltemp[which(mzxmltemp$seqNum==i),"precursorCharge"]-1)*1.007276567,2)
            } else {
              deltatemp <- round(pepfragments[1,2]-(mzxmltemp[which(mzxmltemp$seqNum==i),"precursorMZ"]*precursorCharge)+(precursorCharge-1)*1.007276567,2)
            }
            delta <- c(delta, deltatemp)
          } 
        }
      }
    }
    scan <- scan[order(abs(delta), decreasing = T)]
    scores <- scores[order(abs(delta), decreasing = T)]
    delta <- delta[order(abs(delta), decreasing = T)]
    results <- list(scores=scores,scan=scan,delta=delta)

    print(paste(length(scan),"peptides selected"))
    print("Pass 2")
    for(i in 1:length(results$scan)){
      if(abs(results$delta[i]) <= deltacutoff){
        msmass <- as.data.frame(peaks(mzxml,which(mzxmltemp$seqNum==results$scan[i])))
        PeptideSpectrum(msmass,pepfragments,tolerance,base,xlim=mzrange,
                      label=paste("Scan", results$scan[i],
                                  "\n Score",results$score[i],
                                  "\n Delta (E-O) =",results$delta[i],"Da"),
                      binwidth=binwidth, numperbin=numperbin)}
    }
    results <<- results
    return(results)
  } else {
    print("Pass 1")
    for(i in ms2scans){
      msmass <- as.data.frame(peaks(mzxml,which(mzxmltemp$seqNum==i)))
      msmass <- msmass[order(msmass[,2],decreasing = T),]
      if(nrow(msmass)>0){
        temp <- PeptideSpectrumNoPlot(msmass,pepfragments, tolerance, base, xlim=mzrange, binwidth=binwidth, numperbin=numperbin)
        if(temp!="no peaks were identified"){
          score <- round(sum(msmass[match(temp$expt_mz,msmass[,1]),2]/msmass[1,2]*100),0)
          scores <- c(scores, score)
          scan <- c(scan, i)
          if(is.na(precursorCharge)==T){
            deltatemp <- round(pepfragments[1,2]-(mzxmltemp[which(mzxmltemp$seqNum==i),"precursorMZ"]*mzxmltemp[which(mzxmltemp$seqNum==i),"precursorCharge"])+(mzxmltemp[which(mzxmltemp$seqNum==i),"precursorCharge"]-1)*1.007276567,2)
          } else {
            deltatemp <- round(pepfragments[1,2]-(mzxmltemp[which(mzxmltemp$seqNum==i),"precursorMZ"]*precursorCharge)+(precursorCharge-1)*1.007276567,2)
          }
          delta <- c(delta, deltatemp)
        }
      }
    }
    
    scan <- scan[order(abs(delta),decreasing = T)]
    scores <- scores[order(abs(delta),decreasing = T)]
    delta <- delta[order(abs(delta), decreasing = T)]
    results <- list(scores=scores,scan=scan,delta=delta)
    
    print(paste(length(scan),"peptides selected"))
    print("Pass 2")
    for(i in 1:length(results$scan)){
      if(abs(results$delta[i]) <= deltacutoff){
        msmass <- as.data.frame(peaks(mzxml,which(mzxmltemp$seqNum==results$scan[i])))
      PeptideSpectrum(msmass,pepfragments,tolerance,base,xlim=mzrange,
                      label=paste0("Scan =", results$scan[i],"; ",
                                  "Score =",results$score[i],"; ",
                                  "delta (E-O) =",results$delta[i],"Da"),
                      binwidth=binwidth, numperbin=numperbin)}
    }
    results <<-results
    return(results)
  }
}

spectrumcompare <- function(file, scan1, fragment1, scan2, fragment2, file2=NA, t=0.5, b=2.5, scannames=c(scan1,scan2), precursorCharge=NA, binwidth=100, numperbin=10){
  filehead <- header(file)
  if(is.na(file2)==F){
    filehead2 <- header(file2)
  }
  firstscan <- data.frame(peaks(file, which(filehead$seqNum==scan1)))
  if(is.na(file2)==T){
    secondscan <- data.frame(peaks(file, which(filehead$seqNum==scan2)))
  } else {
    secondscan <- data.frame(peaks(file2, which(filehead2$seqNum==scan2)))
  }
  
  firstscan <- firstscan[order(firstscan[,2],decreasing=T),]
  x <- min(firstscan[,1])
  while(x<max(firstscan[,1])){
    range <- which(firstscan[,1]>x & firstscan[,1]<(x+binwidth))
    if(length(range)!=0){
      if(length(range)>numperbin)
        temp <- firstscan[range,][1:numperbin,]
      else temp <- firstscan[range,]
      firstscan <- firstscan[-range,]
      firstscan <- rbind(firstscan,temp) 
    }
    x <- x+binwidth
  }
  
  secondscan <- secondscan[order(secondscan[,2],decreasing=T),]
  x <- min(secondscan[,1])
  while(x<max(secondscan[,1])){
    range <- which(secondscan[,1]>x & secondscan[,1]<(x+100))
    if(length(range)!=0){
      if(length(range)>10)
        temp <- secondscan[range,][1:10,]
      else temp <- secondscan[range,]
      secondscan <- secondscan[-range,]
      secondscan <- rbind(secondscan,temp) 
    }
    x <- x+100
  }
  
  firstscanmatch <- PeptideSpectrumNoPlot(firstscan, fragment1, t=t, b=b, binwidth=binwidth, numperbin=numperbin)
  secondscanmatch <- PeptideSpectrumNoPlot(secondscan, fragment2, t=t, b=b, binwidth=binwidth, numperbin=numperbin)
  negfirstscan <- firstscanmatch
  negfirstscan$expt_intensity <- firstscanmatch$expt_intensity*-1
  scan12 <- rbind(negfirstscan,secondscanmatch)
  scan12$colour <- ifelse(substr(scan12$ms2type,2,2)=="b",1,2)
  scan12$colour <- factor(scan12$colour, levels=c(1,2), labels=c("b-ion", "y-ion"))
  labeladjust <- ifelse(scan12$expt_intensity<0, -8, 8)
  
  matches <- vector()
  first_mz <- vector()
  second_mz <- vector()
  meanvalue <- vector()
  for(i in 1:nrow(firstscanmatch)){
    thesame <- which(firstscanmatch[i,"ms2type"]==secondscanmatch[,"ms2type"])
    if(length(thesame)!=0){
      matches <- c(matches, as.character(firstscanmatch[i,"ms2type"]))
      first_mz <- c(first_mz,ifelse(duplicated(firstscanmatch[,"ms2type"])[i]==T,mean(firstscanmatch[which(firstscanmatch[,"ms2type"]==firstscanmatch[i,"ms2type"]),1]),firstscanmatch[i,1]))
      second_mz <- c(second_mz,ifelse(length(thesame)>1, second_mz <- mean(secondscanmatch[thesame,1]),secondscanmatch[thesame,1]))
      meanvalue <- c(meanvalue,ifelse(duplicated(firstscanmatch[,"ms2type"])[i]==T | length(thesame)>1, T, F))
      }
    }
  difference <- first_mz - second_mz
  matches <- data.frame(type=matches, first_mz=first_mz,second_mz=second_mz, MzDifference=round(difference,4), meanvalue=meanvalue)
  matches <- matches[order(matches$type),]
  matches <- matches
  
  scans <- header(file)$seqNum
  
  if(is.na(precursorCharge)==T){
    mass1 <- header(file)[which(scans==scan1),"precursorMZ"] * header(file)[which(scans==scan1),"precursorCharge"] - header(file)[which(scans==scan1),"precursorCharge"] * 1.007276567
    if(is.na(file2)==T){
      mass2 <- header(file)[which(scans==scan2),"precursorMZ"] * header(file)[which(scans==scan2),"precursorCharge"] - header(file)[which(scans==scan2),"precursorCharge"] * 1.007276567
    } else {
      mass2 <- header(file2)[which(scans==scan2),"precursorMZ"] * header(file2)[which(scans==scan2),"precursorCharge"] - header(file2)[which(scans==scan2),"precursorCharge"] * 1.007276567
    }
  } else {
    mass1 <- header(file)[which(scans==scan1),"precursorMZ"] * precursorCharge - precursorCharge * 1.007276567
    if(is.na(file2)==T){
      mass2 <- header(file)[which(scans==scan2),"precursorMZ"] * precursorCharge - precursorCharge * 1.007276567
    } else {
      mass2 <- header(file2)[which(scans==scan2),"precursorMZ"] * precursorCharge - precursorCharge * 1.007276567
    }
  }
  
  subtitle <- paste("Left =", scannames[1], ", Right =", scannames[2], "\nDelta =", round(mass1-mass2,4),"Da")
  
  firstscan <- data.frame(mz=firstscan[,1],Intensity=firstscan[,2])
  firstscan$Intensity <- -1*firstscan$Intensity/max(firstscan$Intensity)*100
  secondscan <- data.frame(mz=secondscan[,1],Intensity=secondscan[,2])
  secondscan$Intensity <- secondscan$Intensity/max(secondscan$Intensity)*100
  
  plot <- ggplot(scan12, aes(y=expt_mz, x=0, xmin=0, xmax=expt_intensity)) +
    geom_errorbarh(data=firstscan, inherit.aes = F, aes(x=0, xmin=0, xmax=Intensity, y=mz), colour="grey") +
    geom_errorbarh(data=secondscan, inherit.aes = F, aes(x=0, xmin=0, xmax=Intensity, y=mz), colour="grey") +
    geom_errorbarh(aes(colour=colour),size=0.25,height=0) +
    geom_text(aes(colour=colour, label=paste0(ms2type, ", ",round(expt_mz,2)," mz"), x=expt_intensity + labeladjust), show.legend=F, size=2) +
    theme_classic() +
    geom_vline(xintercept=0) +
    xlab("Normalised Intensity") +
    ylab("mz") +
    ggtitle("", subtitle=subtitle) +
    xlim(c(min(scan12$expt_intensity)-15, max(scan12$expt_intensity)+15)) +
    scale_colour_manual(values=c("red","blue")) +
    theme(legend.title=element_blank())
  
  print(matches)
  return(plot)
}