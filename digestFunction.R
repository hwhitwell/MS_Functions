library(OrgMassSpecR)

digest <- function (sequence, enzyme = "trypsin", missed = 0, IAA = TRUE, 
          N15 = FALSE, custom = list(), minaa = 0, maxaa = Inf) 
{
  if(is.na(as.numeric(substr(sequence,6,6)))==F){
        sequence <- readLines(paste0("https://www.uniprot.org/uniprot/",sequence,".fasta"))
        header <- sequence[1]
        sequence <- sequence[-1]
        sequence <- paste0(sequence,collapse="")
        print(header)
      }
  
  seq_vector <- strsplit(sequence, split = "")[[1]]
  end_position <- length(seq_vector)
  if (enzyme == "trypsin") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] == 
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_string <- gsub("KP", "!P", seq_string)
    seq_string <- gsub("RP", "!P", seq_string)
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "trypsin.strict") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] == 
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "pepsin") {
    if (seq_vector[end_position] == "F" | seq_vector[end_position] == 
        "L" | seq_vector[end_position] == "W" | seq_vector[end_position] == 
        "Y" | seq_vector[end_position] == "A" | seq_vector[end_position] == 
        "E" | seq_vector[end_position] == "Q") {
      seq_vector[end_position] <- "!"
    }
    stop <- grep("F|L|W|Y|A|E|Q", seq_vector)
    start <- stop + 1
  }
  if (enzyme =="LysN"){
    seq_vector <- sequence
    seq_vector <- strsplit(seq_vector,"")[[1]]
    stop <- grep("K|m|d|t", seq_vector) - 1
    start <- grep("K|m|d|t",seq_vector)
  }
  if (enzyme == "LysC") {
    if (seq_vector[end_position] == "K") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "gluC.E") {
    if (seq_vector[end_position] == "E") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("E", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "gluC.ED") {
    if (seq_vector[end_position] == "E" | seq_vector[end_position] == 
        "D") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("E|D", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "ArgC") {
    if (seq_vector[end_position] == 
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("R", seq_vector)
    start <- stop + 1
  }
  if (enzyme != "trypsin" & enzyme != "trypsin.strict" & enzyme != 
      "pepsin" & enzyme != "LysN" & enzyme != "LysC" & enzyme !="gluC.E" & enzyme != "gluC.ED" & enzyme != "ArgC") 
    stop("undefined enzyme, defined enzymes are LysC, LysN, trypsin, trypsin.strict, gluC.E, gluC.ED, ArgC and pepsin")
  if (length(stop) == 0) 
    warning("sequence does not contain cleavage sites")
  if (missed > length(stop)) 
    stop("number of specified missed cleavages is greater than the maximum possible")
  cleave <- function(sequence, start, stop, misses) {
    peptide <- substring(sequence, start, stop)
    mc <- rep(misses, times = length(peptide))
    result <- data.frame(peptide, start, stop, mc, stringsAsFactors = FALSE)
    return(result)
  }
  start <- c(1, start)
  stop <- c(stop, end_position)
  results <- cleave(sequence, start, stop, 0)
  if (missed > 0) {
    for (i in 1:missed) {
      start_tmp <- start[1:(length(start) - i)]
      stop_tmp <- stop[(1 + i):length(stop)]
      peptide <- cleave(sequence, start_tmp, stop_tmp, 
                        i)
      results <- rbind(results, peptide)
    }
  }
  C <- 12
  H <- 1.0078250321
  O <- 15.9949146221
  S <- 31.97207069
  N <- ifelse(N15 == TRUE, 15.0001088984, 14.0030740052)
  proton <- 1.007276466
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
    if (residue =="m")
      mass = C * 7 + H * 14 + N * 2 + O
    if (residue =="d")
      mass = C * 8 + H * 16 + N * 2 + O
    if (residue =="t")
      mass = C * 9 + H * 18 + N * 2 + O
    if (residue == "C" & IAA == FALSE) 
      mass = C * 3 + H * 5 + N + O + S
    if (residue == "C" & IAA == TRUE) 
      mass <- ifelse(N15 == FALSE, C * 5 + H * 8 + N * 
                       2 + O * 2 + S, C * 5 + H * 8 + N + 14.0030740052 + 
                       O * 2 + S)
    if (length(custom) != 0) 
      for (i in 1:length(custom$code)) if (residue == custom$code[i]) 
        mass = custom$mass[i]
      return(mass)
  }
  mz <- vector("list", length = nrow(results))
  for (i in 1:nrow(results)) {
    peptide_vector <- strsplit(results$peptide[i], split = "")[[1]]
    peptide_mass <- sum(sapply(peptide_vector, residueMass))
    mz[[i]] <- round((peptide_mass + H * 2 + O + (c(1, 2, 
                                                    3, 4) * proton))/c(1, 2, 3, 4), digits = 4)
  }
  mz <- as.data.frame(do.call("rbind", mz))
  names(mz) <- c("mz1", "mz2", "mz3", "mz4")
  results <- cbind(results, mz)
  results <- results[order(nchar(results$peptide)),]
  results <- results[which(nchar(results$peptide)>minaa & nchar(results$peptide)<maxaa),]
  return(results)
}

# digest <- function(protein,enzyme="Trypsin",minaa=0,maxaa=Inf){
#   enzymes <- c("Trypsin","LysN")
#   if(enzyme %in% enzymes == F)
#     return(paste(enzyme,"not defined"))
#   
#   if(is.na(as.numeric(substr(protein,6,6)))==F){
#     protein <- readLines(paste0("https://www.uniprot.org/uniprot/",protein,".fasta"))
#     header <- protein[1]
#     protein <- protein[-1]
#     protein <- paste0(protein,collapse="")
#     print(header)
#   }
#   
#   if(enzyme=="Trypsin"){
#     protein <- gsub("K","K.",protein)
#     protein <- gsub("R","R.",protein) 
#   }
#   
#   if(enzyme=="LysN"){
#     protein <- gsub("K",".K",protein) 
#   }
#   
#   protein <- unlist(strsplit(protein,"[.]"))
#   protein <- data.frame(sequence=protein)
#   M <- apply(protein, 1, function(x) MonoisotopicMass(ConvertPeptide(x)))
#   M1H <- apply(protein, 1, function(x) MonoisotopicMass(ConvertPeptide(x),charge=1))
#   M2H <- apply(protein, 1, function(x) MonoisotopicMass(ConvertPeptide(x),charge=2))
#   M3H <- apply(protein, 1, function(x) MonoisotopicMass(ConvertPeptide(x),charge=3))
#   M4H <- apply(protein, 1, function(x) MonoisotopicMass(ConvertPeptide(x),charge=4))
#   protein <- cbind(protein, M, M1H, M2H, M3H, M4H)
#   protein$sequence <- as.character(protein$sequence)
#   protein <- protein[order(nchar(protein$sequence)),]
#   protein <- protein[which(nchar(protein$sequence)>minaa & nchar(protein$sequence)<maxaa),]
#   
#   return(protein)
# }
