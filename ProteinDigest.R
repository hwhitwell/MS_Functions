source("C:/users/hwhitwel.CE-HWHITWEL/OneDrive - Imperial College London/ProteinDigest/digestFunction.R")
#For K{me} use "m", for K{2me} use "d", for K{3me} use "t"

# sequence <- "P68431" #This can be a single letter AA sequence OR a Uniprot accession number
sequence <- "GEAKSQVFmTNHFQKSTVR"
enzyme <- "LysN" #"trypsin", "LysC", "LysN", "pepsin", "trypsin.strict", "gluC.E", "gluC.ED", "ArgC"
minaa <- 1
maxaa <- Inf
IAA <- F
N15 <- F
missed = 3
#acustom <- list(code=c("K","g"),mass=c(MonoisotopicMass(ConvertPeptide("K"))-MonoisotopicMass(ListFormula("H2O"))+MonoisotopicMass(ListFormula("C6N3O6H")),
                                      MonoisotopicMass(ConvertPeptide("G"))-MonoisotopicMass(ListFormula("H2O"))+MonoisotopicMass(ListFormula("C6N3O6H"))))
custom <- list()

View(digest(sequence,enzyme,missed, IAA, N15, custom, minaa=minaa,maxaa=maxaa))
