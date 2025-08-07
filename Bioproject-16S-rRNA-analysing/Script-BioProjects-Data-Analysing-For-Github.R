################################

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

library(seqinr)
library(dada2)

pD="/Users/jean-baptiste/Desktop/"
pS<-"/Users/jean-baptiste/Dropbox/Lichenibacterium/Data-16s"
pG<-"/Users/jean-baptiste/Dropbox/Lichenibacterium/Available-Genome-Ressources/Genomes"


#16S references from genomes only, taxonomy refined from whole-genome phylogeny (see script Script-Genome-Analysis.R)

setwd(pG)
SEQsimp<-read.fasta("SEQ16s_genomics.fasta") 

#Families and Genera for each reference sequence
iREF<-c(1:length(SEQsimp))

TAX<-sapply(iREF,function(x){
	tax<-gsub("RH_AL","RHAL",names(SEQsimp)[x])
	tax<-unlist(strsplit(tax,split="[;]"))
	tax<-unlist(strsplit(tax,split="_"))
	return(paste(tax[c(2,3)],collapse="_"))
})


#Reference sequences pasted to check how dada2 classify them

 SEQsimp2<-sapply(c(1:length(SEQsimp)),function(x){
 	return(paste(unlist(SEQsimp[[x]]),collapse=""))
 })

refSILVA<-assignTaxonomy(SEQsimp2, paste(pD,"/silva_nr99_v138.2_toSpecies_trainset.fa.gz",sep=""), multithread=8)
rownames(refSILVA)<-names(SEQsimp)

setwd(pS)


write.table(refSILVA,"SEQ16s_genomics_classif-by_SILVA.txt")


#Retrieve all ASVs for NS=NSx

NSx=3000

ASV<-sapply(BP,function(x){
	print(x)
	PBx<-unlist(strsplit(x,split='_'))
	TFx=as.numeric(PBx[2])
	TRx=as.numeric(PBx[3])
	
	pW<-paste("/Users/jean-baptiste/Desktop/", PBx[1],sep="")

	setwd(pW)	
	
	seqtab.nochimR<-read.table(paste(c(paste(PBx,collapse="_"),
			"-rarefied_NS=",
			NSx,".txt"),collapse=""),header=T)
	
	return(rownames(seqtab.nochimR))
	#return(dim(seqtab.nochimR))
})

ASV<-unique(unlist(ASV))

####Taxonomic classification using Dada2 with SILVA v138.2 database

print(date())

taxaSILVA <- assignTaxonomy(ASV, paste(pD,"/silva_nr99_v138.2_toSpecies_trainset.fa.gz",sep=""), multithread=8)

print(date())

setwd(pS)

write.table(taxaSILVA,paste("ASV-taxo-SILVA_NS=",NSx,".txt",sep=""))

####Taxonomic classification using Dada2 with reference genomes

print(date())

taxaGEN <- assignTaxonomy(ASV, paste(pS,"SEQ16s_genomics-TAXA.fasta",sep="/"), multithread=8)

print(date())

setwd(pS)

write.table(taxaGEN,paste("ASV-taxo-GENOMES_NS=",NSx,".txt",sep=""))



############################
#refine taxonomy only for PAH ASV

library(seqinr)
library(dada2)

pD="/Users/jean-baptiste/Desktop/"
pS<-"/Users/jean-baptiste/Dropbox/Lichenibacterium/Data-16s"
pG<-"/Users/jean-baptiste/Dropbox/Lichenibacterium/Available-Genome-Ressources/Genomes"

#16S references from genomes only, taxonomy refined from whole-genome phylogeny (see script Script-Genome-Analysis.R)

setwd(pG)
SEQsimp<-read.fasta("SEQ16s_genomics.fasta") 


###########Retrieve ASV taxonomy based on SILVA and ref genomes using assignTaxonomy function
setwd(pS)

NSx=3000
taxaSILVA<-read.table(paste("ASV-taxo-SILVA_NS=",NSx,".txt",sep=""),header=T)
taxaGEN <-read.table(paste("ASV-taxo-GENOMES_NS=",NSx,".txt",sep=""),header=T)

ASV<-rownames(taxaSILVA)

rownames(taxaSILVA)<-NULL
rownames(taxaGEN)<-NULL

###Assign taxonomy of ASV classified by SILVA as Beijerinckiaceae by using 16S sequence from reference genomes as reference

#Place RHAL1 in Beijerinckiaceae and JAJXWB01 in Methylocystaceae: even if it appears as a divergent groups according to phylogenomics, it form a monophyletic group with it and could barely distinguished by 16S

TAX<-ifelse(TAX=="RHAL1_NA","Beijerinckiaceae_RHAL1",TAX)
TAX<-ifelse(TAX=="JAJXWB01_NA","Methylocystaceae_JAJXWB01",TAX)

#Format for classification
TAXf<-paste("Bacteria;Pseudomonadota;Alphaproteobacteria;Hyphomicrobiales",gsub("_",";",TAX),sep=";")

setwd(pS)

write.fasta(SEQsimp,names= TAXf,"SEQ16s_genomics-TAXA.fasta")

#Reassign only "Beijerinckiaceae" ASVs using 16S sequences from genomes 

iB<-c(grep("Beijerinckiaceae",taxaSILVA[,5]), grep("Phreatobacter",taxaSILVA[,6]),grep("Alsobacter",taxaSILVA[,6]))


###with assignTaxonomy
print(date())

taxaB <- assignTaxonomy(ASV[iB], paste(pS,"SEQ16s_genomics-TAXA.fasta",sep="/"), multithread=8)

taxaB <- assignTaxonomy(ASV[iB], paste(pS,"SEQ16s_genomics-TAXA.fasta",sep="/"),outputBootstraps=T, multithread=8)


taxaS<- assignTaxonomy(ASV[iB], paste(pD,"/silva_nr99_v138.2_toSpecies_trainset.fa.gz",sep=""),outputBootstraps=T, multithread=8)


taxaBboot<-taxaB$boot
rownames(taxaBboot)<-NULL

taxaBtax<-taxaB$tax
rownames(taxaBtax)<-NULL


taxaSboot<-taxaS$boot
rownames(taxaSboot)<-NULL

taxaStax<-taxaS$tax
rownames(taxaStax)<-NULL


print(date())

rownames(taxaB)<-NULL

#Compare taxonomy with SILVA and references

COMP<-sapply(sort(unique(apply(taxaSILVA[iB,c(5,6)],1,paste,collapse="_"))),function(x){
	iX=subset(iB,apply(taxaSILVA[iB,c(5,6)],1,paste,collapse="_")==x)
	sapply(sort(unique(apply(taxaBtax[,c(5,6)],1,paste,collapse="_"))),function(x){
		iY=subset(iB,apply(taxaBtax[,c(5,6)],1,paste,collapse="_")==x)
		return(length(intersect(iY,iX)))
	})
})

COMP<-t(t(summary(as.factor(apply(cbind(taxaSILVA[iB,c(5,6)],
	taxaBtax[,c(5,6)]),1,paste,collapse="_")),maxsum=1000)))


write.table(COMP,"Comparison-taxo-SILVA-refined.txt")

####Rules for assignation 

##Cleaning of SILVA results
#In genus: 
#replace "alphaI cluster" and "Neo-b11" by NA, 

GENs<-taxaSILVA[iB,6]
FAMs<-taxaSILVA[iB,5]

#replace Methylorubrum by Methylobacterium
GENs<-gsub("Methylorubrum","Methylobacterium", GENs)

#In family (mostly previously assigned to Beijerinckiaceae): 
# assign "alphaI cluster" and "Neo-b11" to NA

FAMs<-ifelse(GENs=="alphaI cluster",NA,FAMs)
GENs<-gsub("alphaI cluster",NA, GENs)

FAMs<-ifelse(GENs=="Neo-b11",NA,FAMs)
GENs<-gsub("Neo-b11",NA, GENs)

# keep Beijerinckia, Methylocapsa, Methylocella, Methyloferula, Methylorosula (absent from genome) and Methylovirgula in Beijerinckiaceae
 
# Alsobacter to Alsobacteraceae (previously to Hyphomicrobiales Incertae Sedis)
FAMs<-ifelse(GENs=="Alsobacter","Alsobacteraceae",FAMs)

# Bosea to Boseaceae
FAMs<-ifelse(GENs=="Bosea","Boseaceae",FAMs)

# Chelatococcus to Chelatococcaceae
FAMs<-ifelse(GENs=="Chelatococcus","Chelatococcaceae",FAMs)

# Methylobacterium, Methylorubrum, Enterovirga, Microvirga, Psychroglaciecola (absent from genomes) and Salinarimonas to Methylobacteriaceae
FAMs<-ifelse(GENs=="Methylobacterium","Methylobacteriaceae",FAMs)
FAMs<-ifelse(GENs=="Enterovirga","Methylobacteriaceae",FAMs)
FAMs<-ifelse(GENs=="Microvirga","Methylobacteriaceae",FAMs)
FAMs<-ifelse(GENs=="Psychroglaciecola","Methylobacteriaceae",FAMs)
FAMs<-ifelse(GENs=="Salinarimonas","Methylobacteriaceae",FAMs)

# Methylocystis and Methylosinus to Methylocystaceae
FAMs<-ifelse(GENs=="Methylocystis","Methylocystaceae",FAMs)
FAMs<-ifelse(GENs=="Methylosinus","Methylocystaceae",FAMs)

# Phreatobacter to Phreatobacteraceae (previously to Hyphomicrobiales Incertae Sedis)
FAMs<-ifelse(GENs=="Phreatobacter","Phreatobacteraceae",FAMs)

# Rhodoblastus to Rhodoblastaceae 
FAMs<-ifelse(GENs=="Rhodoblastus","Rhodoblastaceae",FAMs)

# Roseiarcus to Roseiarcaceae
FAMs<-ifelse(GENs=="Roseiarcus","Roseiarcaceae",FAMs)

# Lichenibacterium to Lichenihabitantaceae
FAMs<-ifelse(GENs=="Lichenibacterium","Lichenihabitantaceae",FAMs)

FAMg<-taxaBtax[,5]
GENg<-taxaBtax[,6]

sum(ifelse(FAMg==FAMs,1,0),na.rm=T) #2501 ASVs for which both methods give the same family, including**
sum(ifelse(GENg==GENs,1,0),na.rm=T) #**1632 ASVs for which both methods give the same genus
sum(ifelse(is.na(FAMg)==T,1,0)) #262 ASVs for which ref genomes-based database give no family, including*
sum(ifelse(is.na(FAMs)==T,1,0)) #329 ASVs for which SILVA database give no family, including*
sum(ifelse(ifelse(is.na(FAMs)==T,1,0)+ifelse(is.na(FAMg)==T,1,0)==2,1,0)) #*65 ASVs for which no method return a family
sum(ifelse(FAMg!=FAMs,1,0),na.rm=T) #62 ASVs for which both methods give different family

FAMg<-ifelse(is.na(FAMg)==T,"NA",FAMg)
FAMs<-ifelse(is.na(FAMs)==T,"NA", FAMs)
GENg<-ifelse(is.na(GENg)==T,"NA", GENg)
GENs<-ifelse(is.na(GENs)==T,"NA", GENs)

COMPsg<-sapply(sort(unique(FAMg)),function(x){
	F1=x
	sapply(sort(unique(FAMs)),function(x){
		F2=x	
	length(intersect(subset(c(1:length(FAMg)),FAMg==F1),
		subset(c(1:length(FAMs)),FAMs==F2)))
	})
})


TAXg<-cbind("Bacteria","Pseudomonadota","Alphaproteobacteria","Hyphomicrobiales", FAMg, GENg)
colnames(TAXg)<-colnames(taxaSILVA)[c(1:6)]

TAXs<-cbind(taxaSILVA[iB,c(1:4)],FAMs, GENs)
colnames(TAXs)<-colnames(taxaSILVA)[c(1:6)]

iO<-setdiff(c(1:length(ASV)),iB)

TAXg <-rbind(TAXg, taxaSILVA[iO,c(1:6)])
TAXs <-rbind(TAXs, taxaSILVA[iO,c(1:6)])

TAXg<-TAXg[order(c(iB,iO)),] #Taxonomy of all ASV from SILVA + Alsobacter, Phreatobacter and Beijerinckiaceae ASVs reclasssified at the Family and Genus level according to 16S sequences extracted from reference genomes
TAXs<-TAXs[order(c(iB,iO)),] #Taxonomy of all ASV from SILVA + Alsobacter, Phreatobacter and Beijerinckiaceae ASVs reasigned to correct famillies (according to recent classification, validated by phylogenomics) but genus not modified (only Methylorubrum assigned to Methylobacterium and "alphaI cluster" and "Neo-b11" replaced by NA)


##################################
#Estimated genus abundance in each BioProject

#Names of bioproject and brief description
ORI<-c("Snow","PRJNA948675_250_250",
"Air","PRJDB13521_250_220",
"Photovoltaic","PRJEB43933_250_250",
"Bark","PRJNA290145_220_200",
"Lichen","PRJNA558995_250_250",
"Moss","PRJEB40635_280_250",
"Lycopodes","PRJNA867490_250_250",
"Forest-Canada","PRJNA729807_280_250",
"Forest-Germany","PRJEB36420_220_220",
"Forest-China","PRJNA1001021-5k_240_240",
"Tillandsia","PRJNA1134710_300_300",
"Kiwi","PRJNA665460_250_250",
"Bamboo","PRJNA718425_280_280",
"Grapevine","PRJEB59055_240_240")

ORI <-matrix(ORI,ncol=2,byrow=T)

#For each Bioproject, only consider Bacteria ASVs and report relative abundance of
# all phylla with relative abundance >= P1% (including Pseudomonadota)
# all Pseudomonadota classes with relative abundance >= P2% (including Alphaproteobacteria)
# all Alphaproteobacteria orders with relative abundance >= P3% (Including Hyphomicrobiales)
# all Hyphomicrobiales famillies with relative abundance >= P4% (Including Beijerinckiaceae, Methylobacteriaceae, Methylocystaceae and Lichenihabitantaceae)
# all Beijerinckiaceae, Methylobacteriaceae, Methylocystaceae and Lichenihabitantaceae genus

P1=0.1 #bacteria phyllq
P2=0.1 #Pseudomonadota classes
P3=0.05 #Alphaproteobacteria orders
P4=0.02 #Hyphomicrobiales families

#Abundance of these taxon per bioproject 
SUMab<-sapply(BP,function(x){
	print(x)
	PBx<-unlist(strsplit(x,split='_'))
	TFx=as.numeric(PBx[2])
	TRx=as.numeric(PBx[3])
	ORIx<-subset(ORI[,1], ORI[,2]==x)
	
	pW<-paste("/Users/jean-baptiste/Desktop/", PBx[1],sep="")

	setwd(pW)	
	
	seqtab.nochimR<-read.table(paste(c(paste(PBx,collapse="_"),
			"-rarefied_NS=",
			NSx,".txt"),collapse=""),header=T)
	
	asvx<-rownames(seqtab.nochimR)
	
	row.names(seqtab.nochimR)<-NULL
	
	#Taxonomy according to SILVA
	taxxs<-t(sapply(c(1:length(asvx)),function(x){
		return(subset(TAXs,ASV== asvx[x]))
	}))
	#Taxonomy according to reference genomes
	taxxg<-t(sapply(c(1:length(asvx)),function(x){
		return(subset(TAXg,ASV== asvx[x]))
	}))
	
	#Only Keep Bacteria ASVs
	seqtab.nochimR<-subset(seqtab.nochimR, unlist(taxxs[,1])=="Bacteria")
	asvx<-subset(asvx, unlist(taxxs[,1])=="Bacteria")
	taxxg <-subset(taxxg, unlist(taxxs[,1])=="Bacteria")
	taxxs <-subset(taxxs, unlist(taxxs[,1])=="Bacteria")
	
	#Number of reads per sample
	NX<-apply(seqtab.nochimR,2,sum)
	
	####Abundance of Bacteria Phylla per sample
	u1<-sort(unique(unlist(taxxg[,2])))
	A1<-sapply(u1,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,2])==x)
		apply(ax,2,sum)/NX
	})
	
	# all phylla with relative abundance >= P1% (including Pseudomonadota)	
	u1<-unique(c("Pseudomonadota",subset(u1,apply(A1,2,mean)>=P1)))

	A1<-sapply(u1,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,2])==x)
		apply(ax,2,sum)/NX
	})
	
	o1<-(1-apply(A1,1,sum))
	
	A1<-cbind(A1,o1)
	
	colnames(A1)<-c(u1,"other_Bacteria")
	
	####Abundance of Pseudomonadota Classes per sample
	u2<-sort(unique(subset(unlist(taxxg[,3]),unlist(taxxg[,2])=="Pseudomonadota")))
	A2<-sapply(u2,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,3])==x)
		apply(ax,2,sum)/NX
	})
	
	# all classes with relative abundance >= P2% (including Alphaproteobacteria)	
	u2<-unique(c("Alphaproteobacteria",subset(u2,apply(A2,2,mean)>=P2)))

	A2<-sapply(u2,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,3])==x)
		apply(ax,2,sum)/NX
	})
	
	o2<-(A1[,grep("Pseudomonadota",colnames(A1))]-apply(A2,1,sum))
	
	A2<-cbind(A2,o2)
	
	colnames(A2)<-c(u2,"other_Pseudomonadota")	

	####Abundance of Alphaproteobacteria orders per sample
	u3<-sort(unique(subset(unlist(taxxg[,4]),unlist(taxxg[,3])=="Alphaproteobacteria")))
	A3<-sapply(u3,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,4])==x)
		apply(ax,2,sum)/NX
	})
	
	# all classes with relative abundance >= P3% (including Hyphomicrobiales)	
	u3<-unique(c("Hyphomicrobiales",subset(u3,apply(A3,2,mean)>=P3)))

	A3<-sapply(u3,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,4])==x)
		apply(ax,2,sum)/NX
	})
	
	o3<-(A2[,grep("Alphaproteobacteria",colnames(A2))]-apply(A3,1,sum))
	
	A3<-cbind(A3,o3)
	
	colnames(A3)<-c(u3,"other_Alphaproteobacteria")	

	####Abundance of Hyphomicrobiales families per sample
	u4g<-sort(unique(subset(unlist(taxxg[,5]),unlist(taxxg[,4])=="Hyphomicrobiales")))
	u4s<-sort(unique(subset(unlist(taxxs[,5]),unlist(taxxs[,4])=="Hyphomicrobiales")))

	A4g<-sapply(u4g,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,5])==x)
		apply(ax,2,sum)/NX
	})
	
	A4s<-sapply(u4s,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxs[,5])==x)
		apply(ax,2,sum)/NX
	})	
	
	
	# all famillies with relative abundance >= P4% (Including Beijerinckiaceae, Methylobacteriaceae, Methylocystaceae and Lichenihabitantaceae)	
	u4g<-unique(c("Beijerinckiaceae","Methylobacteriaceae",
		"Methylocystaceae","Lichenihabitantaceae",subset(u4g,apply(A4g,2,mean)>=P4)))
	u4s<-unique(c("Beijerinckiaceae","Methylobacteriaceae",
		"Methylocystaceae","Lichenihabitantaceae",subset(u4s,apply(A4s,2,mean)>=P4)))

	u4g<-subset(u4g, u4g!="NA")
	u4s<-subset(u4g, u4g!="NA")

	A4g<-sapply(u4g,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,5])==x)
		apply(ax,2,sum)/NX
	})
	
	A4s<-sapply(u4s,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxs[,5])==x)
		apply(ax,2,sum)/NX
	})	
	
	o4s<-(A3[,grep("Hyphomicrobiales",colnames(A3))]-apply(A4s,1,sum))
	o4g<-(A3[,grep("Hyphomicrobiales",colnames(A3))]-apply(A4g,1,sum))

	
	A4g<-cbind(A4g,o4g)
	A4s<-cbind(A4s,o4s)	
	
	colnames(A4g)<-c(u4g,"other_Hyphomicrobiales")	
	colnames(A4s)<-c(u4s,"other_Hyphomicrobiales")	
	
	#Abundance of genera of interest within Beijerinckiaceae, Methylobacteriaceae, Methylocystaceae and Lichenihabitantaceae

	u5<-c("RHAL1","Methylobacterium","Enterovirga","JAJXWB01",
		"Lichenibacterium","Lichenifustis","Lichenihabitans")

	A5g<-sapply(u5,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,6])==x)
		apply(ax,2,sum)/NX
	})
	
	A5s<-sapply(u5,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxs[,6])==x)
		apply(ax,2,sum)/NX
	})	
	
	########
	
	V1<-round(apply(A1,2,mean),3)
	V2<-round(apply(A2,2,mean),3)
	V3<-round(apply(A3,2,mean),3)
	V4g<-round(apply(A4g,2,mean),3)
	V4s<-round(apply(A4s,2,mean),3)
	V5g<-round(apply(A5g,2,mean),3)
	V5s<-round(apply(A5s,2,mean),3)	

	return(list(c(V1,V2,V3,V4g,V5g)))
})

#Taxon names present at least once in a bioproject
SUMab <-sort(summary(as.factor(unlist(sapply(c(1:length(SUMab)),function(x){
	return(names(SUMab[[x]]))
})))),decreasing=T)

T1<-intersect(names(SUMab),TAXg[,2])
T2<-intersect(names(SUMab),TAXg[,3])
T3<-intersect(names(SUMab),TAXg[,4])
T4<-intersect(names(SUMab),TAXg[,5])
T5<-intersect(names(SUMab),TAXg[,6])

#Recalculate abundance for all these taxons
SUMab2<-sapply(BP,function(x){
	print(x)
	PBx<-unlist(strsplit(x,split='_'))
	TFx=as.numeric(PBx[2])
	TRx=as.numeric(PBx[3])
	ORIx<-subset(ORI[,1], ORI[,2]==x)
	
	pW<-paste("/Users/jean-baptiste/Desktop/", PBx[1],sep="")

	setwd(pW)	
	
	seqtab.nochimR<-read.table(paste(c(paste(PBx,collapse="_"),
			"-rarefied_NS=",
			NSx,".txt"),collapse=""),header=T)
	
	asvx<-rownames(seqtab.nochimR)
	
	row.names(seqtab.nochimR)<-NULL
	
	#Taxonomy according to SILVA
	taxxs<-t(sapply(c(1:length(asvx)),function(x){
		return(subset(TAXs,ASV== asvx[x]))
	}))
	#Taxonomy according to reference genomes
	taxxg<-t(sapply(c(1:length(asvx)),function(x){
		return(subset(TAXg,ASV== asvx[x]))
	}))
	
	#Only Keep Bacteria ASVs
	seqtab.nochimR<-subset(seqtab.nochimR, unlist(taxxs[,1])=="Bacteria")
	asvx<-subset(asvx, unlist(taxxs[,1])=="Bacteria")
	taxxg <-subset(taxxg, unlist(taxxs[,1])=="Bacteria")
	taxxs <-subset(taxxs, unlist(taxxs[,1])=="Bacteria")
	
	#Number of reads per sample
	NX<-apply(seqtab.nochimR,2,sum)
	
	####Abundance of Bacteria Phylla per sample
	u1<-T1
	A1<-sapply(u1,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,2])==x)
		apply(ax,2,sum)/NX
	})
	
	o1<-(1-apply(A1,1,sum))
	
	A1<-cbind(A1,o1)
	
	colnames(A1)<-c(u1,"other_Bacteria")
	
	####Abundance of Pseudomonadota Classes per sample
	u2<-T2
	A2<-sapply(u2,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,3])==x)
		apply(ax,2,sum)/NX
	})
	
	o2<-(A1[,grep("Pseudomonadota",colnames(A1))]-apply(A2,1,sum))
	
	A2<-cbind(A2,o2)
	
	colnames(A2)<-c(u2,"other_Pseudomonadota")	

	####Abundance of Alphaproteobacteria orders per sample
	u3<-T3
	A3<-sapply(u3,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,4])==x)
		apply(ax,2,sum)/NX
	})
	
	o3<-(A2[,grep("Alphaproteobacteria",colnames(A2))]-apply(A3,1,sum))
	
	A3<-cbind(A3,o3)
	
	colnames(A3)<-c(u3,"other_Alphaproteobacteria")	

	####Abundance of Hyphomicrobiales families per sample
	u4g<-T4
	u4s<-T4

	A4g<-sapply(u4g,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,5])==x)
		apply(ax,2,sum)/NX
	})
	
	A4s<-sapply(u4s,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxs[,5])==x)
		apply(ax,2,sum)/NX
	})	
		
	o4s<-(A3[,grep("Hyphomicrobiales",colnames(A3))]-apply(A4s,1,sum))
	o4g<-(A3[,grep("Hyphomicrobiales",colnames(A3))]-apply(A4g,1,sum))

	
	A4g<-cbind(A4g,o4g)
	A4s<-cbind(A4s,o4s)	
	
	colnames(A4g)<-c(u4g,"other_Hyphomicrobiales")	
	colnames(A4s)<-c(u4s,"other_Hyphomicrobiales")	
		
	####Abundance of genera of interest
	u5g<-T5
	u5s<-T5

	A5g<-sapply(u5g,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxg[,6])==x)
		apply(ax,2,sum)/NX
	})
	
	A5s<-sapply(u5s,function(x){
		ax<-subset(seqtab.nochimR,unlist(taxxs[,6])==x)
		apply(ax,2,sum)/NX
	})	


	setwd(pS)
	
	colnames(A1)<-paste(colnames(TAXg)[2],colnames(A1),sep="=")
	colnames(A2)<-paste(colnames(TAXg)[3],colnames(A2),sep="=")
	colnames(A3)<-paste(colnames(TAXg)[4],colnames(A3),sep="=")
	colnames(A4g)<-paste(paste(colnames(TAXg)[5],"G",sep=""),colnames(A4g),sep="=")
	colnames(A4s)<-paste(paste(colnames(TAXg)[5],"S",sep=""),colnames(A4s),sep="=")
	colnames(A5g)<-paste(paste(colnames(TAXg)[6],"G",sep=""),colnames(A5g),sep="=")
	colnames(A5s)<-paste(paste(colnames(TAXg)[6],"S",sep=""),colnames(A5s),sep="=")	
	
	write.table(cbind(A1,A2,A3, A4g, A4s, A5g, A5s),paste("SUM_Tax-",x,"-June2025.txt",sep=""))

	V1<-round(apply(A1,2,mean),3)
	V2<-round(apply(A2,2,mean),3)
	V3<-round(apply(A3,2,mean),3)
	V4g<-round(apply(A4g,2,mean),3)
	V4s<-round(apply(A4s,2,mean),3)
	V5g<-round(apply(A5g,2,mean),3)
	V5s<-round(apply(A5s,2,mean),3)
	
	return(c(V1,V2,V3,V4g,V4s,V5g,V5s))
})


setwd(pS)

write.table(SUMab2,"Summary-Per-BioProject-June2025.txt")

########Per category

#Names of bioproject and brief description
ORI<-c("Snow","PRJNA948675_250_250",
"Air","PRJDB13521_250_220",
"Photovoltaic","PRJEB43933_250_250",
"Bark","PRJNA290145_220_200",
"Lichen","PRJNA558995_250_250",
"Moss","PRJEB40635_280_250",
"Lycopodes","PRJNA867490_250_250",
"Forest-Canada","PRJNA729807_280_250",
"Forest-Germany","PRJEB36420_220_220",
"Forest-China","PRJNA1001021-5k_240_240",
"Tillandsia","PRJNA1134710_300_300",
"Kiwi","PRJNA665460_250_250",
"Bamboo","PRJNA718425_280_280",
"Grapevine","PRJEB59055_240_240")

ORI <-matrix(ORI,ncol=2,byrow=T)

setwd(paste(pS,"Bioprojects",sep="/"))


METAdata<-read.csv("Sample-INFO.csv",sep=";")

###Substrate * Environment * Country

LOC<-METAdata$Location

LOC<-gsub('Tsukuba - Japan','Asia',LOC)
LOC<-gsub('Leipzig - Germany','Europe',LOC)
LOC<-gsub('Au_k\002uluhei_i - Iceland','Europe',LOC)
LOC<-gsub('Itatiba - Brazil','South_America',LOC)
LOC<-gsub('Sorocaba - Brazil','South_America',LOC)
LOC<-gsub('Tulln - Austria','Europe',LOC)
LOC<-gsub('Graz - Austria','Europe',LOC)
LOC<-gsub('China: Beijing','Asia',LOC)
LOC<-gsub('China: Guangdong province','Asia',LOC)
LOC<-gsub('China: Hainan province','Asia',LOC)
LOC<-gsub('China: Heilongjiang province','Asia',LOC)
LOC<-gsub('China: Hunan province','Asia',LOC)
LOC<-gsub('China: Jilin province','Asia',LOC)
LOC<-gsub('China: Zhejiang province','Asia',LOC)
LOC<-gsub('Jaboticabal City - Brazil ','South_America',LOC)
LOC<-gsub('Portugal: Northwest near Vila Verde','Europe',LOC)
LOC<-gsub('China:Yaan','Asia',LOC)
LOC<-gsub('China: hunan','Asia',LOC)
LOC<-gsub('Austria:Johnsbach','Europe',LOC)
LOC<-gsub('Colombia: PNN Chingaza','South_America',LOC)
LOC<-gsub('Colombia: PNN Los Nevados','South_America',LOC)
LOC<-gsub('Canada: Northern Quebec - Uapishka','North_America',LOC)
LOC<-ifelse(LOC=='Canada: Station Biologique des Laurentides (QC)','North_America',LOC)
LOC<-ifelse(LOC=='Canada: Mont-Saint-Hilaire (QC)','North_America',LOC)

C12<-paste(METAdata$Substrate,METAdata$Environment)

C12<-gsub(" Forest","",C12)
C12<-gsub(" top layer Boreal forest","",C12)
C12<-gsub(" top layer Skidoo trail","",C12)
C12<-gsub(" top layer Road","",C12)
C12<-gsub(" top layer Tundra","",C12)
C12<-gsub(" top layer Road with boreal forest","",C12)
C12<-gsub(" with boreal forest","",C12)
C12<-gsub("Epiphyte leaves ","Epiphyte ",C12)
C12<-gsub("Epiphyte roots ","Epiphyte ",C12)
C12<-gsub("Leaf endophyte Greenhouse","Endophyte Greenhouse ",C12)
C12<-gsub("Stem endophyte Greenhouse","Endophyte Greenhouse ",C12)
C12<-gsub("Root endophyte Greenhouse","Endophyte Greenhouse ",C12)

C12<-paste(C12,LOC)

C12<-ifelse(is.na(METAdata$Type)==T, C12 ,ifelse(METAdata$Type=="cDNA","cDNA", C12))

#Categories to keep
C12u<-sort(setdiff(C12,c("Positive control Temperate NA","Leaves Temperate NA","cDNA")))
nC12u<-sapply(C12u,function(x){
	return(length(subset(C12, C12==x)))
})

#####Retrieve average relative abundance for remaining categories

CAT<-matrix(unlist(sapply(C12u,function(x){
	
	Rx<-subset(METAdata$Run,C12==x)
	BPx<-subset(METAdata$BioProject,C12==x)
	
	return(list(t(cbind(x,BPx,sapply(BPx,function(x){
		ORI[grep(x,ORI[,2]),1]
	}),Rx))))
})),ncol=4,byrow=T)

colnames(CAT)<-c("Type","BioProject","Source","Run")

CAT<-CAT[order(CAT[,4]),]

#Remove runs for which there is no abundance (filtered out or not available)

iR<-sapply(CAT[,4],function(x){
	print(x)
	BPx<-BP[grep(subset(CAT[,2], CAT[,4]==x),BP)]
	setwd(pS)
	SUMx<-read.table(paste("SUM_Tax-", BPx,"-June2025.txt",sep=""))	
	length(intersect(x,rownames(SUMx)))
})

CAT<-subset(CAT,iR==1)

CAT<-gsub("Forest-Germany","Forest", CAT)
CAT<-gsub("Forest-China","Forest", CAT)
CAT<-gsub("Forest-Canada","Forest", CAT)
CAT<-ifelse(CAT=="Snow","Surface",CAT)
CAT<-ifelse(CAT=="Photovoltaic","Surface",CAT)

#Retrieve abundance for each run 
SUMab3<-t(sapply(CAT[,4],function(x){
	print(x)
	BPx<-BP[grep(subset(CAT[,2], CAT[,4]==x),BP)]
	
	setwd(pS)
	SUMx<-read.table(paste("SUM_Tax-", BPx,"-June2025.txt",sep=""))	
	
	unlist(subset(SUMx,rownames(SUMx)==x))
}))

#Retrieve average abundance for each CAT
CATu<-sort(unique(paste(CAT[,1],CAT[,3],CAT[,2],sep=",")))

SUMab4<-sapply(CATu,function(x){
	SUMx<-subset(SUMab3,paste(CAT[,1],CAT[,3],CAT[,2],sep=",")==x)
	c(length(SUMx[,1]),apply(SUMx,2,mean))
})


setwd(pS)

write.table(SUMab4,"Summary-Per-Biome-June2025.txt")

