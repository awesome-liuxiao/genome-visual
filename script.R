library(msa)

system.file("tex", "texshade.sty", package = "msa")
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences

myFirstAlignment <- msa(mySequences)
myFirstAlignment

msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

texi2pdf("myFirstAlignment.tex", clean=TRUE)