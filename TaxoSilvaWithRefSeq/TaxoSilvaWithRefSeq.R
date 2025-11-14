################################

#Taxonomic assignation of RefSeq 16S rRNA sequences from Hyphomicrobiales types species + RH_AL1 and JAJXBW01: comparie SILVA 1138.1 and SILVA 138.2
library(seqinr)
library(dada2)

pD="/Users/jean-baptisteleducq/Desktop/"
pS<-"/Users/jean-baptisteleducq/Dropbox/Lichenibacterium/For-Submission-to-SAM/Resubmissiom-oct-2025/TaxoSilvaWithRefSeq"

setwd(pS)
SEQ16S<-read.fasta("RefSeq16S.fasta") 

SEQ16Sc<-sapply(c(1:length(SEQ16S)),function(x){
	return(paste(unlist(SEQ16S[[x]]),collapse=""))
})

#### With SILVA v 138.1

refSILVA<-assignTaxonomy(SEQ16Sc, paste(pD,"silva_nr_v138_train_set.fa",sep=""))

rownames(refSILVA)<-names(SEQ16S)

setwd(pS)

write.table(refSILVA,"RefSeq16S_classif-by_SILVA138-1.txt")

#### With SILVA v 138.2

refSILVA<-assignTaxonomy(SEQ16Sc, paste(pD,"silva_nr99_v138.2_toSpecies_trainset.fa.gz",sep=""))

rownames(refSILVA)<-names(SEQ16S)

setwd(pS)

write.table(refSILVA,"RefSeq16S_classif-by_SILVA138-2.txt")
