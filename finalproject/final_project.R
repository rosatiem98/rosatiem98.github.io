library(tidyverse)
library(seqinr)
library(Biostrings)
library(msa)
library(phangorn)
library(ShortRead)
library(ape)
library(ggplot2)
library(ggtree)
library(ggmsa)
library(zoom)
library(dplyr)
memory.limit()
#Practice with sreads

reads=readFasta("plant_ITS (1).fasta")
readssequences=sread(reads) 
substr(readssequences,1,5)
table(substr(readssequences,1,5))

myreads=reads[substr(readssequences,1,5)=="TAACA"]

myreads

dict=DNAStringSet(substr(readssequences,1,5))

hits=vcountPattern("TAACA", dict,
                   max.mismatch = 1, with.indels = TRUE)
sum(hits)

sread(reads[hits])


myreads2=reads[substr(readssequences,1,5)=="GTCCA"]

myreads2

dict=DNAStringSet(substr(readssequences,1,5))

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

#Using own data
mySequencefile <- system.file("plant_ITS (1).fasta", package = "Biostrings")
file.exists("plant_ITS (1).fasta")


dna <- readDNAStringSet("./plant_ITS (1).fasta")
myfasta <- read.fasta(file = "./plant_ITS (1).fasta")

#subsetting by bp length

counts <- stopifnot(all(getLength(myfasta)))
counts <- getLength(myfasta)
df <- as.data.frame(counts)

dist <- table(df$counts)
df2 <- as.data.frame(dist)
dna
#save as a file
df3 <- dna[-c(4802,4806,4711,4760,6694,4888,4887,4889,4892,4893,4890,4891,4810,4709,4758,2,44,4811,4897,4898,4712)]

alignment <- ggmsa(df3[c(1:10)]) # error message

alignment

#Alignment with msa
alignment2 <- msa(df3[c(1:10)])

alignment2

saveRDS(alignment2,"./alignment2.RDS")

msa10 <- readRDS(file="./alignment2.RDS")


print(msa10, show="complete")



phang.align10=as.phyDat(msa10,type="DNA")

#model test
mt <- modelTest(phang.align)

print(mt)

dm <- dist.ml(phang.align, model="JC69")


saveRDS(dm,"./msa10_distance.RDS")

msa10_UPGMA <- upgma(dm)
msa10_NJ <- NJ(dm)

plot(msa10_UPGMA, main="UPGMA")

plot(msa10_NJ, main="NJ")

parsimony(msa10_UPGMA, phang.align)
parsimony(msa10_NJ, phang.align)


msa10_optim <- optim.parsimony(msa10_NJ, phang.align)

msa10_pratchet <- pratchet(phang.align)


plot(msa10_optim)



plot(msa10_pratchet)


#Buidlding a tree with Maximum Likelihood

fit = pml(msa10_NJ, data=phang.align)
print(fit)


fitJC <- optim.pml(fit,model="JC", rearrangement = "stochastic")

logLik(fitJC)

#bootstrapping
bs <- bootstrap.pml(fitJC, bs=100, optNi=TRUE)

plotBS(midpoint(fitJC$tree), bs, p=50, type="p")

write.tree(bs,file = "./bootstrap_msa10.tree")

#Zahn example

treeNJ <- NJ(dm)

saveRDS(treeNJ,"./msa10_treeNJ.RDS")

fit = pml(treeNJ, data=phang.align)

saveRDS(fit,"./msa10_fit_treeNJ.RDS")

fit$tree$tip.label <- msa10

BiocManager::install("ggtree")


#Large align 

alignment3 <- msa(df3[c(1:400)])
alignment3
saveRDS(alignment3,"./alignment3.RDS")

msa400 <- readRDS(file="./alignment3.RDS")

print(msa400, show="complete")

phang.align400=as.phyDat(msa400,type="DNA")

#model test
mt400 <- modelTest(phang.align400)

print(mt400)

saveRDS(mt400,file="./mt400_modeltest")

dm400 <- dist.ml(phang.align400, model="JC69")


saveRDS(dm400,"./msa10_distance.RDS")

msa400_UPGMA <- upgma(dm400)
msa400_NJ <- NJ(dm400)

plot(msa400_UPGMA, main="UPGMA")

plot(msa400_NJ, main="NJ")

parsimony(msa400_UPGMA, phang.align400)
parsimony(msa400_NJ, phang.align400)


msa400_optim <- optim.parsimony(msa400_NJ, phang.align400)

msa400_pratchet <- pratchet(phang.align400)


plot(msa400_optim)



plot(msa400_pratchet)


#Buidlding a tree with Maximum Likelihood

fit400 = pml(msa400_NJ, data=phang.align400)
print(fit)


fitJC400 <- optim.pml(fit,model="JC", rearrangement = "stochastic")

logLik(fitJC400)

#bootstrapping
bs400 <- bootstrap.pml(fitJC400, bs=100, optNi=TRUE)

plotBS(midpoint(fitJC400$tree), bs400, p=50, type="p") #fix names

write.tree(bs400,file = "./bootstrap_msa400.tree")

#alignment with extraction sequences
dna2 <- readAAStringSet("./Research sequences/DNA extraction sequences.fasta")
test <- readDNAStringSet("./Research sequences/1D96DAA000.fasta")

t <- "fasta$"
test2 <- readFasta("./Research sequences/",t)

extalignment2 <- ggmsa(dna2)
saveRDS(extalignment2,"./ggmsa_align")

library(ggmsa)
extalignment <- ggmsa(f)
library(zoom)
a <- extalignment2
zm(a)

#ggmsa alignment with large fasta file 
dna <- readDNAStringSet("./plant_ITS (1).fasta")
df3 <- dna[-c(4802,4806,4711,4760,6694,4888,4887,4889,4892,4893,4890,4891,4810,4709,4758,2,44,4811,4897,4898,4712)]

largealign <- df3[1:100]
largeggmsa <- ggmsa(largealign)
largealign
saveRDS(largealign,"./largealign_ggmsa.RDS")

readRDS("./largealign_ggmsa.RDS")

#100 sequnce alignment with msa
msa100 <- msa(df3[1:100])
msa100

saveRDS(msa100,"./msa100.RDS")
msa100 <- readRDS("./msa100.RDS")
phang.align100=as.phyDat(msa100,type="DNA")

#model test
mt100 <- modelTest(phang.align100)

print(mt100)

saveRDS(mt100,file="./mt100_modeltest")

dm100 <- dist.ml(phang.align100, model="JC")


saveRDS(dm100,"./msa100_distance.RDS")
dm100 <- readRDS("./msa100_distance.RDS")
msa100_UPGMA <- upgma(dm100)
msa100_NJ <- NJ(dm100)

plot(msa100_UPGMA, main="UPGMA")

plot(msa100_NJ, main="NJ")

parsimony(msa100_UPGMA, phang.align100)
parsimony(msa100_NJ, phang.align100)


msa100_optim <- optim.parsimony(msa100_NJ, phang.align100)

msa100_pratchet <- pratchet(phang.align100)


plot(msa100_optim)



plot(msa100_pratchet)


#Buidlding a tree with Maximum Likelihood

fit100 = pml(msa100_NJ, data=phang.align100)
print(fit100)


fitJC100 <- optim.pml(fit100,model="JC", rearrangement = "stochastic")
saveRDS(fitJC100,"./fitJC100.RDS")

fitJC100 <- readRDS("./fitJC100.RDS")

logLik(fitJC100)

#bootstrapping
bs100 <- bootstrap.pml(fitJC100, bs=100, optNi=TRUE)

plotBS(midpoint(fitJC100$tree), bs100, p=50, type="p") #fix names and save 

giantalign <- msa(df3)

#Edited file itself from MEGA format to format of "DNA Extraction sequences" then added to file
Gymtest <- readDNAStringSet("./Gym Alignment.fasta")

gymtesta <- msa(Gymtest)
gymtesta

gggym <- ggmsa(Gymtest)
gggym

#Combined Research sequences with Angio and Gymno (-1 angio-bad DNA)

example <- readDNAStringSet("./Research sequences/DNA extraction sequences.fasta")

ggmsaex <- ggmsa(example,start=250, end=300 )

ggmsaex
saveRDS(ggmsaex,"./ggmsa_resaerch_alignment_small.RDS")
msaex <- msa(example)

saveRDS(msaex,"./research_sequences_alignment.RDS")
phang.alignex=as.phyDat(msaex,type="DNA")

#model test
researchmsa <- modelTest(phang.alignex)

print(researchmsa)

saveRDS(researchmsa,file="./mtex_modeltest")

dmex <- dist.ml(phang.alignex, model="JC")


saveRDS(dmex,"./msaex_distance.RDS")
dmex <- readRDS("./msaex_distance.RDS")
msaex_UPGMA <- upgma(dmex)
msaex_NJ <- NJ(dmex)

plot(msaex_UPGMA, main="UPGMA")

plot(msaex_NJ, main="NJ")

parsimony(msaex_UPGMA, phang.align100)
parsimony(msaex_NJ, phang.align100)


msaex_optim <- optim.parsimony(msaex_UPGMA, phang.alignex)

msaex_pratchet <- pratchet(phang.alignex)


plot(msaex_optim)

plot(msaex_pratchet)


#Buidlding a tree with Maximum Likelihood

fitex = pml(msaex_UPGMA, data=phang.alignex)
print(fitex)


fitJCex <- optim.pml(fitex,model="JC", rearrangement = "stochastic")
saveRDS(fitJCex,"./finalproject/fitJCex.RDS")

fitJCex <- readRDS("./finalproject/fitJCex.RDS")

logLik(fitJCex)

#bootstrapping
bsrex <- bootstrap.pml(fitJCex, bs=100, optNi=TRUE)

saveRDS(bsrex,"./finalproject/boot_strp_research_seqs.RDS")
plotBS((fitJCex$tree), bsrex, p=50, type="p", bs.col = "red") #fix names and save 



?plotBS

#subset large dataset to 300-500 dataset

#read in fasta file 
myspecies <- readDNAStringSet("./finalproject/Interested Species ITS 2.fasta")

myspmsa <- msa(myspecies)
myspmsa

shortggmsa <- ggmsa(myspecies[196:206], start=100, end = 150)
shortggmsa

saveRDS(myspmsa,"./myspmsa_alignment.RDS")

myspmsa <- readRDS("./myspmsa_alignment.RDS")
#subset by excluding duplicaded

duplicated(myspecies)
un <- unique(myspecies)
unspecies <- msa(un)

align399 <- readRDS("./finalproject/subset_msa_align.RDS")

align399[1:10]

msa399 <- as.phyDat(align399,type="DNA")

mt399 <- modelTest(msa399)

mt399 <- readRDS("./modeltest_399.RDs")

dist399 <- dist.ml(align399, model="JC69")
dist399 <- readRDS("./dm399_distance.RDS")

msa399_UPGMA <- upgma(dist399)
saveRDS(msa399_UPGMA,"./finalproject/UPGMA_399.RDS")

msa399_NJ <- NJ(dist399)
saveRDS(msa399_NJ,"./finalproject/NJ_399.RDS")

plot(msa399_UPGMA)
plot(msa399_NJ)

parsimony(msa399_UPGMA,msa399)

parsimony(msa399_NJ,msa399)

msa399_optim <- optim.parsimony(msa399_UPGMA,msa399)

msa399_optim <- readRDS("./msa399_optim")

msa399_pratchet <- pratchet(msa399)
msa399_pratchet <- readRDS("./msa399_pranchet")

plot(msa399_optim)

plot(msa399_pratchet)

fit399 <- pml(msa399_UPGMA,msa399)

fit399JC <- optim.pml(fit399, model="JC", rearrangement = "stochastic")

fit399JC <- readRDS("./fit399JC.RDS")

logLik(fit399JC)

bs399 <- bootstrap.pml(fit399JC,bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
bs399 <- readRDS("./boot_strp_research_seqs399.RDS")
tre <- plotBS(midpoint(fit399JC$tree), bs399, p=50, type="p")

genus195 <- tre$tip.label[1:195] %>% str_split(" ") %>% map_chr(2)
species195 <- tre$tip.label[1:195] %>% str_split(" ") %>% map_chr(3)
names195 <- paste0(genus195," ",species195)

names195

genus197 <- tre$tip.label[197] %>% str_split(" ") %>% map_chr(2)
species197 <- tre$tip.label[197] %>% str_split(" ") %>% map_chr(3)
names197 <- paste0(genus197," ",species197)

genus315 <- tre$tip.label[199:315] %>% str_split(" ") %>% map_chr(2)
species315 <- tre$tip.label[199:315] %>% str_split(" ") %>% map_chr(3)
names315 <- paste0(genus315," ",species315)

genus321 <- tre$tip.label[320:321] %>% str_split(" ") %>% map_chr(2)
species321 <- tre$tip.label[320:321] %>% str_split(" ") %>% map_chr(3)
names321 <- paste0(genus321," ",species321)

genus324 <- tre$tip.label[323:324] %>% str_split(" ") %>% map_chr(2)
species324 <- tre$tip.label[323:324] %>% str_split(" ") %>% map_chr(3)
names324 <- paste0(genus324," ",species324)

genus328 <- tre$tip.label[328] %>% str_split("") %>% map_chr(1)
species328 <- tre$tip.label[328] %>% str_split("Â") %>% map_chr(2)
names328 <- paste0(genus328," ",species328)


genus330 <- tre$tip.label[330] %>% str_split(" ") %>% map_chr(2)
species330 <- tre$tip.label[330] %>% str_split(" ") %>% map_chr(3)
names330 <- paste0(genus330," ",species330)

genus396 <- tre$tip.label[336:396] %>% str_split(" ") %>% map_chr(2)
species396 <- tre$tip.label[336:396] %>% str_split(" ") %>% map_chr(3)
names396 <- paste0(genus396," ",species396)

genus399 <- tre$tip.label[399] %>% str_split(" ") %>% map_chr(2)
species399 <- tre$tip.label[399] %>% str_split(" ") %>% map_chr(3)
names399 <- paste0(genus399," ",species399)





tre$tip.label[1:195] <- names195
tre$tip.label[197] <- names197
tre$tip.label[199:315] <- names315
tre$tip.label[320:321] <- names321
tre$tip.label[323:324] <- names324
tre$tip.label[328] <- names328
tre$tip.label[330] <- names330
tre$tip.label[336:396] <- names396
tre$tip.label[399] <- names399


tre$tip.label



tre <- plot(tre, font=0.5, edge.width=c(1))

subtreeplot(tre)

tree=midpoint(tre)

rtree <- ggplot(tre) + geom_tree() + theme_tree() +
  geom_tiplab(align = TRUE, linesize = 0.2)

test_tree <- ggtree(tre, layout="fan", open.angle = 120)



names <- tre$tip.label
df <- data.frame(names)

library(tidyr)
