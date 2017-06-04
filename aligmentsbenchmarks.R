# source("https://bioconductor.org/biocLite.R")
# biocLite("msa")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
# biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
# biocLite("rbenchmark")

library(rbenchmark)
library(msa)

################### fasta file loading ###################
# hg19 or hg38 is big big big file.
# That's file is not working in MacBook `13 2015 Late...
# We need to select the necessary files.

mySequenceFile <- "/Users/sigmadream/Works/R/RV100/BBA0001.tfa"
# mySequenceFile <- system.file("examples", "exampleDNA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences

######### ClustalW, ClustalOmega, Muscle benchmark #######
myFirstAlignment <- msa(mySequences)
oneFirstAlignment <- msa(mySequences, "ClustalW")
oneFirstAlignment
twoFirstAlignment <- msa(mySequences, "ClustalOmega")
twoFirstAlignment
threeFirstAlignment <- msa(mySequences, "Muscle")
threeFirstAlignment

benchmark(myFirstAlignment <- msa(mySequences),
          oneFirstAlignment <- msa(mySequences, "ClustalW"), 
          twoFirstAlignment <- msa(mySequences, "ClustalOmega"), 
          threeFirstAlignment <- msa(mySequences, "Muscle"), order=NULL)

######### Benchmark results #######
#                                                   test replications elapsed relative
# 1                  myFirstAlignment <- msa(mySequences)          100  76.059    2.426
# 2     oneFirstAlignment <- msa(mySequences, "ClustalW")          100  78.468    2.503
# 3 twoFirstAlignment <- msa(mySequences, "ClustalOmega")          100 122.395    3.904
# 4     threeFirstAlignment <- msa(mySequences, "Muscle")          100  31.351    1.000

