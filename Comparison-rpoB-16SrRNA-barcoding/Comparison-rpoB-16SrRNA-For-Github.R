#########rpoB barcoding: retrive ASV table from PRJNA729807 = Leducq et al 2022 (mBio)

pS<-"/Users/jean-baptisteleducq/Dropbox/Lichenibacterium/Data-16s"
pG<-"/Users/jean-baptisteleducq/Dropbox/Lichenibacterium/Available-Genome-Ressources/Genomes"
library(seqinr)
library(parallel)
library(dada2)

#1 Retrieve rarefied table (NS=3000) from mBio paper (Dataset S1, table k-rpoB-ASV-details)

setwd(pS)

rpoBtable<-read.table("rpoB-PRJNA729807.txt",header=T) #rpoB ASV table


#2 Assign ASV sequences using rpoB sequences retrieved from genomes (rpoB_genomics.fasta) as refererence - compare with previous assignation 
setwd(pG)

REFrpoB<-read.fasta("rpoB_genomics.fasta") #rpoB sequence from reference genomes

#Only keep complete rpoB sequences (used to build the rpoB phylogeny)
KEEP<-names(read.fasta("rpoB_genomics_aligned.fas"))

IX<-sapply(c(1:length(REFrpoB)),function(x){
	nx=gsub(";","",names(REFrpoB)[x])
	return(length(intersect(KEEP,nx)))
})
#Index of reference sequences that can be keeped (complete)
KEEP<-subset(c(1:length(REFrpoB)),IX==1)
#Index of reference sequences that should be manually checked (see if there's a complete overlap with ASV sequences)
CHECK<-subset(c(1:length(REFrpoB)),IX==0)

#Reformate names of reference sequences

NewN<-matrix(unlist(strsplit(names(REFrpoB),split=";")),ncol=3,byrow=T)

NewN<-cbind("Bacteria","Pseudomonadota","Alphaproteobacteria","Hyphomicrobiales", NewN[,c(2:3)])

NewN<-cbind(NewN[,c(1:5)],t(sapply(c(1:length(NewN[,1])),function(x){
	sx<-c(unlist(strsplit(NewN[x,6],split="_")),"NA")
	return(ifelse(sx=="sp","NA",sx)[1:2])
})))

NewN<-sapply(c(1:length(NewN[,1])),function(x){
	return(paste(NewN[x,],collapse=";"))	
})

NewN<-gsub("_","", NewN)

REFrpoB_KEEP<-sapply(KEEP,function(x){
	return(list(unlist(REFrpoB[[x]])))
})

REFrpoB_CHECK<-sapply(CHECK,function(x){
	return(list(unlist(REFrpoB[[x]])))
})

setwd(pG)

write.fasta(REFrpoB_KEEP,names= NewN[KEEP],"rpoB_genomics-NewN.fasta")

write.fasta(REFrpoB_CHECK,names= NewN[CHECK],"rpoB_genomics-to-check.fasta")

####In MEGA: Manually select sequences from REFrpoB_CHECK that overlap ASV sequences, and transfer them to REFrpoB_KEEP (create a new file)


ASV<-rpoBtable[,2]

print(date())

TAXArpoB <- assignTaxonomy(ASV, paste(pG,"rpoB_genomics-NewN2.fasta",sep="/"), multithread=8,minBoot = 50)

print(date())

rownames(TAXArpoB)<-NULL

TAXArpoB <-ifelse(is.na(TAXArpoB)==T,"NA", TAXArpoB)

#Taxonomy from the previous paper
ORDtable<-ifelse(is.na(rpoBtable[,6])==T,
	"NA",rpoBtable[,6])
FAMtable<-ifelse(is.na(rpoBtable[,7])==T,
	"NA",rpoBtable[,7])
GENtable<-ifelse(is.na(rpoBtable[,8])==T,
	"NA",rpoBtable[,8])
	
#Compare both taxonomies

U1<-sort(unique(paste(ORDtable ,
	FAMtable, GENtable,sep="_")))
U2<-sort(unique(paste(TAXArpoB[,5], 
	TAXArpoB[,6],sep="_"))) 

COMP12<-sapply(U1,function(x){
	X=x
	z<-subset(c(1:length(ASV)), 
		paste(ORDtable,FAMtable, 
		GENtable,sep="_") ==X)
	return(sapply(U2,function(x){
		length(intersect(z,subset(c(1:length(ASV)),
		paste(TAXArpoB[,5], 
		TAXArpoB[,6],sep="_")==x)))
	}))
})

setwd(pS)

write.table(COMP12,
	"rpoBtax_comparison-Ogier-REF.txt")
	
#Because reference sequences only included phyllosphere-associated Hyphomicrobiales, we only performed reclassification for ASVs previously identified as Rhizobiales (former Hyphomicrobiales), hence excluding Caulobacterales, Sphingomonadales and Rhodospirillales. Within Rhizobiales, we excluded ASVs previously assigned to other families (Aurantimonadaceae, Bradyrhizobiaceae, Brucellaceae, Hyphomicrobiaceae, Phyllobacteriaceae and Rhizobiaceae)	

#For each ASV: If previous classification: 0; if new classification : 1

CLASSnew<-ifelse(ORDtable!="Rhizobiales",0,
	ifelse(FAMtable=="Aurantimonadaceae",0,
	ifelse(FAMtable =="Bradyrhizobiaceae",0,
	ifelse(FAMtable =="Brucellaceae",0,
	ifelse(FAMtable =="Hyphomicrobiaceae",0,
	ifelse(FAMtable =="Phyllobacteriaceae",0,
	ifelse(FAMtable =="Rhizobiaceae",0,1	
)))))))

ORDnew<-ifelse(CLASSnew==0, ORDtable,
	"Hyphomicrobiales")
ORDnew<-ifelse(ORDnew =="Rhizobiales", 
	"Hyphomicrobiales", ORDnew)
	
FAMnew<-ifelse(CLASSnew==0, FAMtable,
	TAXArpoB[,5])
	
GENnew<-ifelse(CLASSnew==0, GENtable,
	TAXArpoB[,6])		

#Index for each forest
iMSH<-grep("MSH",colnames(rpoBtable))
iSBL<-grep("SBL",colnames(rpoBtable))

#Recalculate abundance per sample

TAXnew<-paste(ORDnew, FAMnew, GENnew,sep="_")

ABnew<-sapply(unique(TAXnew),function(x){
	ABx<-subset(rpoBtable[,sort(c(iMSH,iSBL))], 
		TAXnew ==x)
	return(apply(ABx,2,sum))	
})

ABnewSBL<-sapply(unique(TAXnew),function(x){
	ABx<-subset(rpoBtable[,sort(c(iSBL))], 
		TAXnew ==x)
	return(apply(ABx,2,sum))	
})

ABnewMSH<-sapply(unique(TAXnew),function(x){
	ABx<-subset(rpoBtable[,sort(c(iMSH))], 
		TAXnew ==x)
	return(apply(ABx,2,sum))	
})

SUMall<-cbind(apply(ABnew,2,mean),
	apply(ABnewMSH,2,mean),
	apply(ABnewSBL,2,mean))/3000
	
colnames(SUMall)<-c("All","MSH","SBL")

SUMall<-subset(SUMall,apply(SUMall,1,max)>=0.01)

barplot(SUMall,col=c("blue","red","cyan",
	"green","grey","yellow2","red3","yellow4"),
	legend=T)

setwd(pS)

write.table(ABnew,"Abundance-per-taxa-rpoB-juilly25.txt")

write.table(cbind(ASV,TAXnew),"rpoB-new-tax-juilly25.txt",row.names=F)

#3 abundance per forest

TAXs<-c("RHAL1","Methylobacterium","Enterovirga",
	"Lichenifustis","Lichenibacterium","JAJXWB01")
	
COLs<-c("blue","red","pink2",
	"cyan3","green3","yellow2")

Fu<-c("MSH","SBL")

Fmin=0.0001

#dnsity (violin plot) parameters
DENS=400
SVP=0.05

#abundance of orders

ABo<-sapply(c("Hyphomicrobiales",
	"Caulobacterales"),function(x){
	X=x
	vx<-grep(X,colnames(ABnew))
	unlist(sapply(c(1:length(Fu)),function(x){
		ABx<-apply(ABnew[grep(Fu[x],
			rownames(ABnew)),vx]/3000,1,sum)
		z<-mean(ABx)
	}))	
})

#abundance of families

ABf<-sapply(c("Methylobacteriaceae",
	"Lichenihabitantaceae"),function(x){
	X=x
	vx<-grep(X,colnames(ABnew))
	unlist(sapply(c(1:length(Fu)),function(x){
		ABx<-apply(ABnew[grep(Fu[x],
			rownames(ABnew)),vx]/3000,1,sum)
		z<-mean(ABx)
	}))	
})

#abundance of genera

setwd(pS)


pdf("rpoB-abundance-July2025.pdf",
	height=5,width=4)

par(mar=c(4,4,1,1),bty="n")

plot(-10,-10,xlim=c(1,3),ylim=c(Fmin,1),
	log="y",las=1,cex.axis=0.8,xlab="",
	ylab="Relative abundance (rpoB barcoding, log scale)",cex.axis=0.6,xaxt="n")



ABg<-sapply(c(1:length(TAXs)),function(x){
	X=x
	vx<-grep(TAXs[X],colnames(ABnew))
	unlist(sapply(c(1:length(Fu)),function(x){
		ABx<-ABnew[grep(Fu[x],
			rownames(ABnew)),vx]/3000
		z<-mean(ABx)	
		ABx<-ifelse(ABx==0,Fmin,ABx)
		#ABx<-subset(ABx,ABx!=0)	
		points(seq(x,x,length.out=length(ABx))
			+sample(seq(-0.025,0.025,
			length.out=1001))[c(1:length(ABx))]
			+X/8,
			ABx,col=COLs[X],pch=19,cex=0.2)
		
		XX<-density(ABx,n= DENS,na.rm=T)$y
		XX=SVP*XX/max(XX)
		YY<-density(ABx,n= DENS,na.rm=T)$x
		
		polygon(
			c(XX,(-XX)[DENS:1])+X/8+x,
			c(YY,YY[DENS:1]),
			col=rgb(col2rgb(COLs[X])[1]/255,
				col2rgb(COLs[X])[2]/255,
				col2rgb(COLs[X])[3]/	
				255,0.3),border=COLs[X],lwd=0.5)	
				
		segments(x+X/8-SVP,median(ABx),x+X/8+SVP,
			median(ABx),lend=2,
			lwd=3,col="black")			
		
		return(z)	
	}))
})

dev.off()

colnames(ABg)<-TAXs

#4 Compare genus / famillies abundance with estimation from 16s barcoding 

pS<-"/Users/jean-baptisteleducq/Dropbox/Lichenibacterium/Data-16s"

library(vegan)

setwd(pS)
ABnew <-read.table("Abundance-per-taxa-rpoB-juilly25.txt",header=T) #rpoB abundance per taxon

SUM16S<-read.table("SUM_Tax-PRJNA729807_280_250-June2025.txt",header=T) #taxon abundance per sample for 16S
N16S<-read.table("Names-16S-PRJNA729807.txt",header=T) #correspondance between SRA and sample names for 16S


#Samples for which both 16S and rpoB were screened
SAMPcomp<-intersect(N16S[,1],
	gsub("[.]","-",rownames(ABnew)))
PLOT<-matrix(unlist(strsplit(SAMPcomp,split="-")),
	ncol=4,byrow=T)[,1]
TIME<-matrix(unlist(strsplit(SAMPcomp,
	split="-")),ncol=4,byrow=T)[,4]
TIME<-gsub("A",1, TIME)
TIME<-gsub("B",1, TIME)
TIME<-gsub("C",2, TIME)
TIME<-gsub("D",2, TIME)
TIME<-gsub("E",3, TIME)
TIME<-gsub("F",3, TIME)
TIME<-gsub("G",4, TIME)
TIME<-gsub("H",4, TIME)

#Taxa present in both 16S and rpoB screens

TAXs<-intersect(unique(unlist(strsplit(
	colnames(ABnew),split="_"))),
	unique(unlist(strsplit(
	colnames(SUM16S),split="[.]"))))

#For each shared taxa, look at correlation in abudance between rpoB nd 16S for samples with possible comparison

COORD<-sapply(SAMPcomp,function(x){
	c(subset(c(1:length(N16S[,1])),N16S[,1]==x),
		subset(c(1:length(rownames(ABnew))),
		gsub("[.]","-",rownames(ABnew))==x))
})

Fmin=0.0001

setwd(pS)

pdf("rpoB-vs-16S-July2025.pdf",
	height=5,width=5)

par(mar=c(4,4,1,1),bty="n")

plot(-10,-10,xlim=c(Fmin,1),ylim=c(Fmin,1),
	log="xy",las=1,cex.axis=0.6,
	xlab="Relative abundance (16S rRNA barcoding, log scale)",
	ylab="Relative abundance (rpoB barcoding, log scale)")

TAXs<-c("RHAL1","Methylobacterium","Enterovirga",
	"Lichenifustis","Lichenibacterium","JAJXWB01")
COLs<-c("blue","red","pink3",
	"cyan3","green2","yellow3")

TEST<-t(sapply(c(1:length(TAXs)),function(x){
	print(x)
	X= TAXs[x]
	a16S<-SUM16S[COORD[1,],
		c(grep(X,colnames(SUM16S)),
		grep(X,colnames(SUM16S)))][,1]
	aRPOB<-apply(as.data.frame(ABnew)[COORD[2,],
		c(grep(X,colnames(ABnew)),
		grep(X,colnames(ABnew)))]/6000,1,sum)
	points(a16S+ Fmin, aRPOB+ Fmin,
		pch=ifelse(PLOT=="SBL",21,22),
		col=COLs[x],cex=0.8,
		bg=rgb(col2rgb(COLs[x])[1]/255,
			col2rgb(COLs[x])[2]/255,
			col2rgb(COLs[x])[3]/255,0.3))
	
	print(X)
	print(anova(lm(aRPOB ~ a16S* PLOT )))
	#return(c(cor.test(a16S, aRPOB)$estimate,
	#	cor.test(a16S, aRPOB)$p.value))
	
	return(list(t(cbind(PLOT, SAMPcomp,X,a16S,aRPOB))))
}))

dev.off()

TEST<-matrix(unlist(TEST),ncol=5,byrow=T)

#Hellinger transformation of relative abundances
HELL<-as.data.frame(decostand(cbind(as.numeric(TEST[,4]),as.numeric(TEST[,5])) ,method="hellinger", MARGIN=2))

MOD16S<-HELL[,1]
MODrpoB<-HELL[,2]
MODplot<-TEST[,1]
MODsamp<-TEST[,2]
MODtax<-TEST[,3]


plot(as.data.frame(decostand(cbind(as.numeric(TEST[,4]),as.numeric(TEST[,5])) ,method="hellinger", MARGIN=2)))

anova(lm(MODrpoB ~ MOD16S* MODplot* MODtax))

anova(lm(MODrpoB ~ MODtax +MOD16S
	+MOD16S:MODtax+MODplot:MODtax ))

#redo a correlation test REP times by randomly permutating 16S abundance values - Keep associations of values within samples and plots

#create matrices of rpoB and 16S values: one for each plot and samples in lines and taxa in column

MATsbl<-subset(cbind(MODsamp, MODtax, MOD16S, MODrpoB),
	MODplot=="SBL")
MATmsh<-subset(cbind(MODsamp, MODtax, MOD16S, MODrpoB),
	MODplot=="MSH")
	

COLs2<-sapply(MODtax,
	function(x){
	return(subset(COLs,TAXs==x))	
})
	
		
	
MATsbl16S<-sapply(sort(unique(MATsbl[,1])),function(x){
	z<-as.numeric(subset(MATsbl[,3],MATsbl[,1]==x))
	taxz<-subset(MATsbl[,2],MATsbl[,1]==x)
	return(z[order(taxz)])
})
MATsblRPOB<-sapply(sort(unique(MATsbl[,1])),function(x){
	z<-as.numeric(subset(MATsbl[,4],MATsbl[,1]==x))
	taxz<-subset(MATsbl[,2],MATsbl[,1]==x)
	return(z[order(taxz)])
})
MATmsh16S<-sapply(sort(unique(MATmsh[,1])),function(x){
	z<-as.numeric(subset(MATmsh[,3],MATmsh[,1]==x))
	taxz<-subset(MATmsh[,2],MATmsh[,1]==x)
	return(z[order(taxz)])
})
MATmshRPOB<-sapply(sort(unique(MATmsh[,1])),function(x){
	z<-as.numeric(subset(MATmsh[,4],MATmsh[,1]==x))
	taxz<-subset(MATmsh[,2],MATmsh[,1]==x)
	return(z[order(taxz)])
})		

cor.test(MODrpoB, MOD16S)
CORreal<-cor.test(unlist(c(MATsblRPOB, MATmshRPOB)), 
	unlist(c(MATsbl16S, MATmsh16S)))$estimate

cSBL<-c(1:length(MATsblRPOB[1,]))
cMSH<-c(1:length(MATmshRPOB[1,]))

REP=10000
Fmin=0.1

setwd(pS)

pdf("rpoB-vs-16S-perm-July2025.pdf",
	height=5,width=10)

par(mar=c(4,4,1,1),bty="n",mfrow=c(1,2))

plot(log10(MOD16S+ Fmin),MODrpoB,col= COLs2,log="",
	pch=ifelse(MODplot =="SBL",21,22),
	las=1,cex.axis=0.8,xlab="16S rRNA (Hellinger transformation, log scale)",ylab="rpoB (Hellinger transformation)",xaxt="n")

xcat<-c(-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4)

axis(1,at=xcat,labels=round(10^xcat-Fmin,2),cex.axis=0.8)
	
legend("bottomright",cex=0.8,text.font=c(3,3,3,3,3,3,1,1), 
	horiz=F,
	legend= c(TAXs,"SBL","MSH"),
	text.col= "black",box.col=NA,
	border=NA,col= c(COLs,"black","black"),
	pch=c(21,21,21,21,21,21,21,22),
	lty=c(0,0,0,0,0,0,1,2),
	y.intersp=c(0.8),
	title="",title.col="black",ncol=1)	

ordiellipse(
	cbind(log10(MOD16S+Fmin),MODrpoB),kind="sd",conf=0.75,
	groups= paste(MODplot, MODtax,sep="_"),
	col=c(COLs[order(TAXs)], COLs[order(TAXs)]),
	alpha= 15,lty=c(2,2,2,2,2,2,2,1,1,1,1,1,1))
	

MAT16Sall<-log10(cbind(MATsbl16S, MATmsh16S)+Fmin)
cALL<-c(1:length(MAT16Sall[1,]))

LMrep1<-sapply(c(1:REP),function(x){
	#randomisation of samples among and within plots
	modx<-lm(c(as.vector(MATsblRPOB), 
		as.vector(MATmshRPOB)) ~ 
		as.vector(MAT16Sall[,sample(cALL)]))
	abline(modx,col=rgb(0.8,0.8,0.8,0.05))		
	#modx$coef
	cor.test(c(as.vector(MATsblRPOB), 
		as.vector(MATmshRPOB)) , 
		as.vector(MAT16Sall[,sample(cALL)]))$estimate
})

LMrep2<-sapply(c(1:REP),function(x){
	#randomisation of samples within plots
	modx<-lm(c(as.vector(MATsblRPOB), 
		as.vector(MATmshRPOB)) ~ 
		log10(c(as.vector(MATsbl16S[,sample(cSBL)]), 
		as.vector(MATmsh16S[,sample(cMSH)]))+Fmin))
	abline(modx,col=rgb(0.2,0.2,0.2,0.05))		
	#modx$coef
	cor.test(c(as.vector(MATsblRPOB), 
		as.vector(MATmshRPOB)), 
		log10(c(as.vector(MATsbl16S[,sample(cSBL)]), 
		as.vector(MATmsh16S[,sample(cMSH)]))+Fmin))$estimate
})


MODreal<-lm(MODrpoB ~ log10(MOD16S+Fmin))



abline(MODreal,col="red",lty=1,lwd=2)	

#Distribution of slope

hist(as.vector(LMrep1),breaks=100,
	col="grey",border=NA,
	las=1,cex.axis=0.8,main="",
	xlab="Correlation coefficient",
	xlim=c(0.50,0.62),ylim=c(0,700))
par(new=T)
hist(as.vector(LMrep2),breaks=50,
	col=rgb(0,0,0,0.5),border=NA,
	las=1,cex.axis=0.8,main="",
	xlab="",xaxt="n",
	xlim=c(0.50,0.62),ylim=c(0,700))	
points(cor.test(MODrpoB , 
	log10(MOD16S+Fmin))$estimate,0,pch=21,cex=2,
	col=rgb(1,0,0),lwd=2,bg="white")	

legend("topright",cex=0.8, 
	horiz=F,
	legend= c("observed","expected (permutations among sites)","expected (permutations within sites)"),
	text.col= "black",box.col=NA,
	border=NA,col= c("red","grey",rgb(0,0,0,0.5)),
	pch=c(21,22,22),
	pt.bg=c("white","grey",rgb(0,0,0,0.5)),
	lty=c(0,0,0),
	y.intersp=c(0.8),
	title="",title.col="black",ncol=1)

	
dev.off()

