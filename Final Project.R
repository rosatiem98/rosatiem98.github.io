library(tidyverse)
library(seqinr)
library(Biostrings)
library(msa)
library(phangorn)
library(ShortRead)
memory.limit()
#Practice with sreads

reads=readFasta("plant_ITS (1).fasta")
sequences=sread(reads) 
substr(sequences,1,5)
table(substr(sequences,1,5))

myreads=reads[substr(sequences,1,5)=="TAACA"]

myreads

dict=DNAStringSet(substr(sequences,1,5))

hits=vcountPattern("TAACA", dict,
                   max.mismatch = 1, with.indels = TRUE)
sum(hits)

sread(reads[hits])


myreads2=reads[substr(sequences,1,5)=="GTCCA"]

myreads2

dict=DNAStringSet(substr(sequences,1,5))

hits2=vcountPattern("GTCCA", dict, 
                   max.mismatch = 1, with.indels = TRUE)

sum(hits2)
sread(reads[hits2])

#Alignmet practice/example
mySequencefile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequencefile)
mySequences

myFirstAlignment <- msa(mySequences)
myFirstAlignment

print(myFirstAlignment, show="complete")

msaPrettyPrint(myFirstAlignment, y=c(164,213), output="asis",
               showNames = "none", showLogo = "none", askForOverwrite = FALSE)
