

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
