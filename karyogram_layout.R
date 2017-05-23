## library(devtools)
## install_github(genomicsclass/ERBS)
library(ERBS)

data(HepG2)
library(GenomeInfoDb)
seqlevels(HepG2, force=TRUE) = paste0("chr", 1:22)

data(GM12878)
seqlevels(GM12878, force=TRUE) = paste0("chr", 1:22)

library(ggbio)
autoplot(HepG2, layout="karyogram", main="ESRRA binding on HepG2")
