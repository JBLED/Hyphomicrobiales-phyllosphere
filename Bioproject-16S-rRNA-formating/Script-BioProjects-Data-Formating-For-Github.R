

#####Retrieve metadata for each bioproject on https://www.ncbi.nlm.nih.gov/Traces/study/


########################################################################################################################################################################################################################
########################################################################################################################################################################################################################
########################################################################################################################################################################################################################
########################################################################################################################################################################################################################

#Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.20")

library(dada2)

#Number of threads that can be used

NC=8

########################################################################################################################################################################################################################
######BioProject PRJDB13521

pW<-"/Users/jean-baptiste/Desktop/PRJDB13521"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)

# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##Only keep bacteria samples
sample.names<-intersect(DESC[grep("-B",DESC$LibraryName),1],subset(DESC[,1],DESC[,4]==DESC[,6]))

#retrieve relevant run info

INFO<-sapply(sample.names,function(x){
	
	return(subset(DESC$LibraryName,DESC[,1]==x))
})
write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")

#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 
 
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-250
#Reverse read
TRx<-220
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))


########################################################################################################################################################################################################################
######BioProject PRJEB36420

pW<-"/Users/jean-baptiste/Desktop/PRJEB36420"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)

# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

#Path to samples
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE)) #forward reads
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE)) #reverse reads

# sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")

#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 
 
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-220
#Reverse read
TRx<-220
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

########################################################################################################################################################################################################################
######BioProject PRJEB40635

pW<-"/Users/jean-baptiste/Desktop/PRJEB40635"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)
# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

#Path to samples
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE)) #forward reads
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE)) #reverse reads

# sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")

#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-280
#Reverse read
TRx<-250
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

########################################################################################################################################################################################################################
######BioProject PRJEB43933

pW<-"/Users/jean-baptiste/Desktop/PRJEB43933"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)

# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

#Path to samples
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE)) #forward reads
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE)) #reverse reads

# sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

INFO<-sapply(sample.names,function(x){
	return(subset(DESC$SampleName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")

#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-250
#Reverse read
TRx<-250
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

########################################################################################################################################################################################################################
######BioProject PRJEB59055

pW<-"/Users/jean-baptiste/Desktop/PRJEB59055"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW


##Only keep samples that have the same number of reads in forward and reverse files!
#sample.names<-subset(DESC[,1],DESC[,6]==DESC[,4])
sample.names<-DESC[,1]

INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")

#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-240
#Reverse read
TRx<-240
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC, matchIDs=T) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))



########################################################################################################################################################################################################################
######BioProject PRJNA1001021

#This is a big dataset so only the 5000 first reads per sample were downloaded

pW<-"/Users/jean-baptiste/Desktop/PRJNA1001021-5k"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW


##Only keep samples that have the same number of reads in forward and reverse files!
sample.names<-subset(DESC[,1],DESC[,6]==DESC[,4])


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$SampleName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")


#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


## truncLen ###this was check first on a couple samples
#forward read
TFx<-240
#Reverse read
TRx<-240
# trimLeft
TLx<-20


getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")


########################################################################################################################################################################################################################
######BioProject PRJNA1134710

pW<-"/Users/jean-baptiste/Desktop/PRJNA1134710"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

#Only keep 16s (Remove ITS)

ix<-grep("16S",DESC[,12])

##Only keep samples that have the same number of reads in forward and reverse files!
sample.names<-subset(DESC[ix,1],DESC[ix,6]==DESC[ix,4])


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$SampleName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")




#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-300
#Reverse read
TRx<-300
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

########################################################################################################################################################################################################################
######BioProject PRJNA665460

pW<-"/Users/jean-baptiste/Desktop/PRJNA665460"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##Only keep samples that have the same number of reads in forward and reverse files!
#sample.names<-subset(DESC[,1],DESC[,6]==DESC[,4])
sample.names<-DESC[,1]


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")



#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-250
#Reverse read
TRx<-250
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))


########################################################################################################################################################################################################################
######BioProject PRJNA718425

pW<-"/Users/jean-baptiste/Desktop/PRJNA718425"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##Remove metagenomic samples
sample.names<-subset(DESC[,1],DESC[,15]=="GENOMIC")


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")



#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-280
#Reverse read
TRx<-280
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))




########################################################################################################################################################################################################################
######BioProject PRJNA867490

pW<-"/Users/jean-baptiste/Desktop/PRJNA867490"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##samples 
sample.names<-subset(DESC[,1],DESC[,6]==DESC[,4])


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")




#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-250
#Reverse read
TRx<-250
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))


########################################################################################################################################################################################################################
######BioProject PRJNA290145

pW<-"/Users/jean-baptiste/Desktop/PRJNA290145"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##samples 
sample.names<-DESC[,1]

INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")


#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-220
#Reverse read
TRx<-200
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))


###Forward and Reverse files were inverted in a certain pattern

z<-seq(1,length(sample.names)-1,length.out=(length(sample.names)/2))

out<-sapply(z,function(x){
	print(x)
	out1 <- filterAndTrim(fnFs[x], filtFs[x], fnRs[x+1], filtRs[x+1], truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE
	out2 <- filterAndTrim(fnFs[x+1], filtFs[x+1], fnRs[x], filtRs[x], truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE
	return(list(t(rbind(out1, out2))))
		
})

out<-matrix(unlist(out),ncol=2,byrow=T)


errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs[c(z,z+1)], filtFs[c(z,z+1)], dadaRs[c(z+1,z)], filtRs[c(z+1,z)], verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out[c(z,z+1),], sapply(dadaFs, getN)[c(z,z+1)], sapply(dadaRs, getN)[c(z+1,z)], sapply(mergers, getN)[c(z,z+1)], rowSums(seqtab.nochim)[c(z,z+1)])
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))



########################################################################################################################################################################################################################
######BioProject PRJNA558995

pW<-"/Users/jean-baptiste/Desktop/PRJNA558995"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##samples 
sample.names<-subset(DESC[,1],DESC[,6]==DESC[,4])


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")




#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-250
#Reverse read
TRx<-250
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))



########################################################################################################################################################################################################################
######BioProject PRJNA948675

pW<-"/Users/jean-baptiste/Desktop/PRJNA948675"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##samples (remove ITS and genomes)
sample.names<-DESC[grep("16S",DESC$LibraryName),1]


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")





#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-250
#Reverse read
TRx<-250
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))


########################################################################################################################################################################################################################
######BioProject PRJNA729807 (16S only)

pW<-"/Users/jean-baptiste/Desktop/PRJNA729807"

setwd(pW)

#List of Fastq files
RUN<-read.table("SRR.numbers",header=T)

#Description of runs
DESC<-read.csv("runinfo.csv",header=T)


# Forward and reverse fastq filenames have format: SRRname_1.fastq.gz and SRRname_2.fastq.gz

path<-pW

##samples (remove rpoB)
sample.names<-DESC[grep("16S",DESC$LibraryName),1]


INFO<-sapply(sample.names,function(x){
	return(subset(DESC$LibraryName,DESC[,1]==x))
})

write.table(cbind(sample.names,as.vector(INFO)),"RunInfo.txt")




#path to samples
fnFs<-paste(path,paste(sample.names,"_1.fastq.gz",sep=""),sep="/")
fnRs <-paste(path,paste(sample.names,"_2.fastq.gz",sep=""),sep="/")


#Check quality profile of 4 randomly chosen libraries if order to set trimming
z<-sample(c(1:length(fnFs)),4)

setwd(path)

pdf("quality_profiles_F.pdf")
plotQualityProfile(fnFs[z])
dev.off()

pdf("quality_profiles_R.pdf")
plotQualityProfile(fnRs[z])
dev.off()

#Your reads must still overlap after truncation in order to merge them later! Your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them. 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## truncLen 
#forward read
TFx<-280
#Reverse read
TRx<-250
# trimLeft
TLx<-20

getN <- function(x) sum(getUniques(x))

print(paste(date(),"; Trim left=", TLx,"; Trunc forward=",TFx,"; Trunc reverse=",TRx))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(TFx,TRx),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=TLx,compress=TRUE, multithread=NC) # On Windows set multithread=FALSE

errF <- learnErrors(filtFs, multithread=NC)
errR <- learnErrors(filtRs, multithread=NC)

setwd(path)
pdf(paste("plotErrors_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

	#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=NC, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=NC, pool="pseudo")

#Concatenate paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=F)#justConcatenate is Optional. Default FALSE. If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them. USE ONLY FOR RPOB, NOT FOR 16S
	#Construct sequence table

seqtab <- makeSequenceTable(mergers)
	#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=NC, verbose=TRUE)

	#Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track,paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

pdf(paste("Reads_Track_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".pdf",sep=""))
par(mar=c(6,4,1,1),bty="n",mfrow=c(1,2))
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
barplot(t(cbind(track[,6],(track[,5]-track[,6]),(track[,2]-track[,5]),(track[,1]-track[,2]))/track[,1]),las=2,cex.names=0.6,cex.axis=0.6,ylab="reads")
dev.off()

seqtab.nochim<-t(seqtab.nochim)

write.table(seqtab.nochim,paste("SeqTableNoChim_TF=",TFx,"_TR=",TRx,"_TL=",TLx,".txt",sep=""))

##########################
##########################
##########################
##########################
##########################
##########################


#16S Bioprojects
BP<-c("PRJDB13521_250_220",
"PRJEB36420_220_220",
"PRJEB40635_280_250",
"PRJEB43933_250_250",
"PRJEB59055_240_240",
"PRJNA1001021-5k_240_240",
"PRJNA1134710_300_300",
"PRJNA665460_250_250",
"PRJNA718425_280_280",
"PRJNA867490_250_250",
"PRJNA290145_220_200",
"PRJNA558995_250_250",
"PRJNA948675_250_250",
"PRJNA729807_280_250")

#########################
# rarefaction: produce a rarefied abundance table for each bioproject and for different thresholds (NS)

TLx=20
Kx=50 #Number of replicate resampling
Lx=20 #number of ressampling value to test in the range 1-total number of sequence in the sample

par(mar=c(4,4,1,1),bty="n")

#Thresholds for rarefaction (to adjust)
NS=c(1000,2000,3000,4000,5000,10000,20000)

#start of the loop (for each project)
sapply(BP,function(x){
	print(x)
	PBx<-unlist(strsplit(x,split='_'))
	TFx=as.numeric(PBx[2])
	TRx=as.numeric(PBx[3])
	
	pW<-paste("/Users/jean-baptiste/Desktop/", PBx[1],sep="")

	setwd(pW)
	
	#Raw abundance table
	seqtab.nochim<-read.table(paste("SeqTableNoChim_TF=",TFx,"_TR=",
		TRx,"_TL=",TLx,".txt",sep=""),header=T)

	#number of ASV and sequences in each sample
	Nx<-t(sapply(c(1:length(seqtab.nochim[1,])),function(x){
		Fx<-as.vector(seqtab.nochim[,x])
		return(c(length(subset(Fx,Fx!=0)),sum(Fx)))
	}))

	#Adjust the range of thresholds
	NSadj<-subset(NS,NS<=max(Nx[,2]))
	
	setwd(pW)
	
	pdf(paste(c("Rarefaction-",paste(PBx,collapse="_"),".pdf"),
		collapse=""),width=5,height=5)
	
	#Plot rarefaction curves
	plot(1,1,
		ylim=c(1,max(Nx[,1])),xlim=c(10,max(Nx[,2])),
		xlab="Sampled sequences",ylab="ASVs",
		cex.axis=0.8,las=1,log="xy",type="l",main=paste(PBx,collapse="_"),
		cex.main=0.8)

	sapply(c(1:length(seqtab.nochim[1,]))
		[order(Nx[,1],decreasing=T)],function(x){
		#print(x)
		Fx<-as.vector(seqtab.nochim[,x])
		Fx<-subset(Fx,Fx!=0)
		Fx<-unlist(sapply(c(1:length(Fx)),function(x){
			return(seq(x,x,length.out=Fx[x]))
		}))
		Rx<-round(seq(1,Nx[x,2],length.out=Lx))
		RRk<-apply(sapply(Rx,function(x){
			#print(x)
			RRx=x
			return(sapply(c(1:Kx),function(x){
				return(length(unique(sample(Fx,
				RRx,replace=F))))
			}))
		}),2,mean)
		points(Rx,RRk,type="l",col=rgb(0,0,0,0.1))
		points(Nx[x,2],Nx[x,1],col="black",pch=21,cex=0.5,bg="white")
	})	

	segments(NSadj,seq(1,1,length.out=length(NSadj)),	
		NSadj,seq(max(Nx[,1]),max(Nx[,1]),
		length.out=length(NSadj)),col="red",lty=2)

	#Number of samples left for each threshold
	NsampLeft<-c(length(Nx[,1]),sapply(NSadj,function(x){
		return(length(subset(Nx[,2], Nx[,2]>=x)))
	}))

	legend("topleft",cex=0.5, horiz=F,
		legend= paste(c("all:",paste("NS=", NSadj,":",sep="")),
		c(NsampLeft)),text.col= "red",box.col=NA,
		border=NA,ncol=1)

	dev.off()

	#Total Number of ASVs
	NR=length(seqtab.nochim[,1])

	#For each NS value, sample randomly NS sequences in each sample and produce a rarefied abundance table
	
	sapply(NSadj,function(x){
		NSx=x
		#Remove samples with no enought reads
		iS<-subset(c(1:length(seqtab.nochim[1,])), Nx[,2]>= NSx)
	 

		seqtab.nochimR<-sapply(iS,function(x){
			#print(x)
			Fx<-as.vector(seqtab.nochim[,x])
			Ix=subset(c(1:length(Fx)),Fx!=0)
			Fx<-subset(Fx,Fx!=0)
			Fx<-unlist(sapply(c(1:length(Fx)),function(x){
				return(seq(Ix[x],Ix[x],length.out=Fx[x]))
			}))

			Fs<-sample(Fx, NSx,replace=F)
			return(
				sapply(c(1: NR),function(x){
					return(length(subset(Fs,Fs==x)))
			}))
		})
		rownames(seqtab.nochimR)<-rownames(seqtab.nochim)
		colnames(seqtab.nochimR)<-colnames(seqtab.nochim)[iS]

		#Remove unsampled ASVs
		seqtab.nochimR <-subset(seqtab.nochimR,
			apply(seqtab.nochimR,1,sum)>=1)

		setwd(pW)

		write.table(seqtab.nochimR,paste(c(paste(PBx,collapse="_"),
			"-rarefied_NS=",
			NSx,".txt"),collapse=""))
			
		print(c(paste(PBx,collapse="_"),NSx,dim(seqtab.nochimR)))	
	})

})	

