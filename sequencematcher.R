# source("C:/users/hwhitwel.CE-HWHITWEL/Documents/GitHub/MS_Functions/sequencematch_functions.R")
source("C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/RScripts/sequencematch_functions.R")

spectrum <- openMSfile("C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/Results/20190130_Melis_Tryp/20190130_ME402Tryp_001.mzXML",backend="Ramp")

optimalspectrum <- NA
minrt <- range(header(spectrum)$retentionTime)[1] #min/max RT in seconds.
maxrt <- range(header(spectrum)$retentionTime)[2]
peptide <- "SQVFtTNHFQK" #synthetic peptide

tolerance <- 0.8 #da (MS2 match tolerance)
base <- 2 #Minimum normalised intensity
IAA <- F #Are cystines iodoacetylated? (T/F)
mzrange <- c(100,2000) #x-axis range - limits the x-axis range on output plots
modifiedaminoacids <- c("m","d","t","k","g","o","a","i","E")
#EZLink <- MonoisotopicMass(ListFormula("C37H67N3O15S"))
aminoacidmass <- 
  #Mono
                  # c(MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+1*MonoisotopicMass(ListFormula("CH3"),isotopes=list(C=13.0033548378, H=2.0141017780))-1*MonoisotopicMass(ListFormula("H")),
                   c(MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+1*MonoisotopicMass(ListFormula("CH2")) + 56.02621,
  #Di
                   MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+2*MonoisotopicMass(ListFormula("CH3"),isotopes=list(C=13.0033548378, H=2.0141017780))-2*MonoisotopicMass(ListFormula("H")),
                   # MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+1*MonoisotopicMass(ListFormula("CH2"))+1*MonoisotopicMass(ListFormula("CH3"),isotopes=list(C=13.0033548378, H=2.0141017780))-1*MonoisotopicMass(ListFormula("H")),
                   # MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+2*MonoisotopicMass(ListFormula("CH2")),
  #Tri
                   #MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+3*MonoisotopicMass(ListFormula("CH3"),isotopes = list(C=13.0033548378, H=2.0141017780))-3*MonoisotopicMass(ListFormula("H")),
                   MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+3*MonoisotopicMass(ListFormula("CH2")),
  #Other
                   MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+MonoisotopicMass(ListFormula("C6N3O6H")),
                   MonoisotopicMass(ConvertPeptide("G"))-MonoisotopicMass(ListFormula("H2O"))+MonoisotopicMass(ListFormula("C6N3O6H")),
                   MonoisotopicMass(ConvertPeptide("M"))-MonoisotopicMass(ListFormula("H2O"))+MonoisotopicMass(ListFormula("O")),
                   MonoisotopicMass(ConvertPeptide("A"))-MonoisotopicMass(ListFormula("H2O"))+2*MonoisotopicMass(ListFormula("CH2")),
                   MonoisotopicMass(ConvertPeptide("I"))-MonoisotopicMass(ListFormula("H2O"))+1*MonoisotopicMass(ListFormula("CH2")),
                   MonoisotopicMass(ConvertPeptide("E"))-MonoisotopicMass(ListFormula("H2O")) + 56.02621)
deltacutoff <- 54 #Max difference between theortical parent ion and oberved (Da)
precursorCharge <- 2 #Option, can specify precursor charge if MS alocated precusor charges are incorrect. Set to NA if not using.
binwidth <- 100
numperbin <- 10

# pdf("C:/users/hwhitwel.CE-HWHITWEL/Desktop/temp.pdf")
SpectrumFinder(mzxml=spectrum,
               minrt=minrt,
               maxrt=maxrt,
               peptide=peptide,
               tolerance=tolerance,
               base=base,
               mzrange=mzrange,
               modifiedaa=modifiedaminoacids,
               aminoacidmass=aminoacidmass,
               bestspectrum=optimalspectrum,
               IAA=IAA,
               N15=F,
               deltacutoff=deltacutoff,
               precursorCharge=precursorCharge,
               binwidth=binwidth,
               numperbin=numperbin)
# dev.off()

