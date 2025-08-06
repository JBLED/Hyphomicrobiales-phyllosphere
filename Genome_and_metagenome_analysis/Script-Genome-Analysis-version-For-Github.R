#packages

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msa")


install.packages("seqinr")
install.packages("ape")
install.packages("ips")


library(parallel)
library(seqinr)
library(msa)

#work directory for genomes and annotation files
pG<-"/Users/jean-baptiste/Dropbox/Lichenibacterium/Available-Genome-Ressources/Genomes"

setwd(pG)

data1<-read.table("data4.txt",header=T) #contain information and path for genomes and metagenome files (RAST annotation file and fasta files)


#Genomes statistics 
TYPE<-c("peg","rna","repeat")

#those hidden lines are just to copy the RAST and genome files from different folders to the same one

#sapply(c(1:length(data1[,1])),function(x){
#	print(x)
#	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T)
#	NX<-unlist(strsplit(data1[x,6],split="/"))
#	NX<-NX[length(NX)]
#	write.table(ANNOx, paste("/Users/jean-baptiste/Dropbox/Lichenibacterium/For-Submission-to-SAM/For_Github/Genome_and_metagenome_analysis/RAST/",NX,sep=""),row.names=F)
#	file.copy(paste(gsub("jean-baptisteleducq","jean-baptiste",data1[x,9]),data1[x,5],sep="/"),
#		paste("/Users/jean-baptiste/Dropbox/Lichenibacterium/For-Submission-to-SAM/For_Github/Genome_and_metagenome_analysis/Genomes/",data1[x,5],sep=""))
#})

STAT<-t(sapply(c(1:length(data1[,1])),function(x){
	print(x)
	STx<-data1[x,1]
	setwd(gsub("jean-baptisteleducq","jean-baptiste",data1[x,9])) #you have to modify the path
	seqx<-read.fasta(data1[x,5])
	#number of contigs
	NCx=length(seqx)
	#contig size
	SCx<-sort(as.numeric(summary(seqx)[,1]),decreasing=T)
	#Total size
	TSx=sum(SCx)
	#N50
	N50x=subset(SCx,(cumsum(SCx)/TSx)>0.5)[1]
	#GC
	GCx<-GC(unlist(seqx))
	
	#% dinucleotides
	z1<-unlist(seqx)
	z2<-z1[2:length(z1)]
	z1<-z1[1:length(z2)]
	z1<-summary(as.factor(as.vector(paste(z1,z2,sep=""))))
	rm(z2)
		
	DType<-rbind(c("at",""),
		c("ta",""),
		c("cg",""),
		c("gc",""),
		c("ac","gt"),
		c("ag","ct"),
		c("tc","ga"),
		c("tg","ca"),
		c("aa","tt"),
		c("cc","gg"))
		
	z1<-apply(cbind(sapply(DType[,1],function(x){
		c(subset(z1,names(z1)==x),0)[1]
	}),
	sapply(DType[,2],function(x){
		c(subset(z1,names(z1)==x),0)[1]
	})),1,sum)
	names(z1)<-c(paste(DType[,1],DType[,2],sep="-"))
	
	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T) #you have to modify the path
	#number and type of annotations
	ANx<-sapply(TYPE,function(x){
		return(length(subset(ANNOx[,3],ANNOx[,3]==x)))
	})
	ANx<-c(ANx,dim(ANNOx)[1]-sum(ANx))
	sort(summary(as.factor(ANNOx[,3])),decreasing=T)
	#number of unique annotations
	AUx<-length(unique(ANNOx[,8]))
	#number of hypothetical proteins
	HYPx<-length(subset(ANNOx[,8],ANNOx[,8]=="hypothetical protein"))
	
	return(c(STx, TSx, NCx, N50x, GCx, ANx, AUx, HYPx,round(z1/sum(z1),4)))
	
}))

colnames(STAT)<-c("isolate","size","contigs","N50","GC","peg","rna","repeat","other","unique","hypothetical",colnames(STAT)[12:length(STAT[1,])])

setwd(pG)

write.table(STAT,"Stat-genomes.txt")

#List of unique gene annotations

ANNOf<-sort(unique(unlist(sapply(c(1:length(data1[,1])),function(x){
	print(x)	
	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T) #you have to modify the path
	return(unique(ANNOx[,8]))
}))))

#Frequence of each annotation

ANNOa<-sapply(c(1:length(data1[,1])),function(x){
	print(x)	
	ANNOx<-summary(as.factor(as.vector(read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T)[,8])),
		maxsum=length(ANNOf)) #you have to modify the path
	ANNOo<-setdiff(ANNOf,names(ANNOx))
	ABx<-seq(0,0,length.out=length(ANNOo))
	names(ABx)<-ANNOo
	
	ANNOx<-c(ANNOx, ABx)
	return(ANNOx[order(names(ANNOx))])
	
})

#Retrieve 16s rRNA sequence
X="SSU rRNA ## 16S rRNA, small subunit ribosomal RNA"

SEQ16S<-sapply(c(1:length(data1[,1])),function(x){
	print(x)	
	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T) #you have to modify the path
	#keep the most abundant and the longest sequence
	SEQx<-summary(as.factor(c(subset(ANNOx[,12],ANNOx[,8]==X),"NA")))
	SEQx<-subset(names(SEQx), SEQx==max(SEQx))
	SEQx<-subset(SEQx ,nchar(SEQx)==max(nchar(SEQx)))[1]
	return(list(t(cbind(x, SEQx))))
})
SEQ16S<-matrix(unlist(SEQ16S),ncol=2,byrow=T)
SEQ16S<-subset(SEQ16S, SEQ16S[,2]!="NA")

NX<-sapply(c(1:length(SEQ16S[,1])),function(x){
	Z<-as.numeric(SEQ16S[x,1])
	STx<-paste(data1[Z,1],data1[Z,2],data1[Z,3],sep=";")
	return(STx)	
})

SEQ16S<-sapply(c(1:length(SEQ16S[,1])),function(x){
	return(list(unlist(strsplit(SEQ16S[x,2],split=""))))	
})

setwd(pG)

write.fasta(SEQ16S,names= NX,"SEQ16s_genomics.fasta")
	
#Retrieve rpoB sequence
X="DNA-directed RNA polymerase beta subunit (EC 2.7.7.6)"

SEQrpoB<-sapply(c(1:length(data1[,1])),function(x){
	print(x)	
	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T) #you have to modify the path
	#keep the most abundant and the longest sequence
	SEQx<-summary(as.factor(c(subset(ANNOx[,12],ANNOx[,8]==X),"NA")))
	SEQx<-subset(names(SEQx), SEQx==max(SEQx))
	SEQx<-subset(SEQx ,nchar(SEQx)==max(nchar(SEQx)))[1]
	return(list(t(cbind(x, SEQx))))
})
SEQrpoB <-matrix(unlist(SEQrpoB),ncol=2,byrow=T)
SEQrpoB <-subset(SEQrpoB, SEQrpoB[,2]!="NA")

NX<-sapply(c(1:length(SEQrpoB[,1])),function(x){
	Z<-as.numeric(SEQrpoB[x,1])
	STx<-paste(data1[Z,1],data1[Z,2],data1[Z,3],sep=";")
	return(STx)	
})

SEQrpoB <-sapply(c(1:length(SEQrpoB[,1])),function(x){
	return(list(unlist(strsplit(SEQrpoB[x,2],split=""))))	
})

write.fasta(SEQrpoB,names= NX,"rpoB_genomics.fasta")

####Matrix of pairwise similarity in gene content between each genome

MAT<-sapply(c(1:length(data1[,1])),function(x){
	AN1<-as.numeric(ANNOa[,x])
	AN1<-ifelse(AN1==0,0,1)
	sapply(c(1:length(data1[,1])),function(x){
		AN2<-as.numeric(ANNOa[,x])
		AN2<-ifelse(AN2==0,0,1)
		AN12<-AN1+AN2
		SHARED=length(subset(AN12, AN12==2))
		TOTAL=length(subset(AN12, AN12!=0))
		return(SHARED/sum(AN1))
	})
})

heatmap(MAT,scale="none",labRow=data1[,2],labCol=data1[,2],margin=c(10,10),cexRow=0.6,cexCol=0.6,RowSideColors=ifelse(data1[,7]=="Genomics","black","white"),ColSideColors=ifelse(data1[,7]=="Genomics","black","white"))

#only keep genes present in at least PC% of complete genomes (Genomics)

iC<-subset(c(1:length(data1[,1])),data1[,7]=="Genomics")

hist(apply(ifelse(ANNOa[,iC]==0,0,1),1,sum),breaks=1000)
plot(sort(apply(ifelse(ANNOa[,iC]==0,0,1),1,sum),decreasing=T))

plot(apply(ifelse(ANNOa[,iC]==0,0,1),1,sum), apply(ifelse(ANNOa==0,0,1),1,sum))

PC=0.8

ANNOc<-subset(ANNOa,apply(ifelse(ANNOa[,iC]==0,0,1),1,sum)>= PC*length(iC))

#only keep genes present in average in less than two copies per complete genomes (Genomics)

ANNOc<-subset(ANNOc,as.vector(apply(ANNOc[,iC],1,mean))<2)

# size of each copy of each gene (metagenomes + genomes)
Sizec<-sapply(c(1:length(data1[,1])),function(x){
	print(x)	
	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T) #you have to modify the path
	sapply(c(1:length(ANNOc[,1])),function(x){
		return(paste(nchar(subset(
			as.vector(ANNOx[,12]),
			as.vector(ANNOx[,8])== rownames(ANNOc)[x])),
			collapse="_"))
	})
})


#size of each gene copy normalized by mean size accross genomes 
SizeM<-sapply(c(1:length(Sizec[,1])),function(x){
	sx<-Sizec[x,]
	#mean size
	mx<-mean(as.numeric(
		unlist(strsplit(sx,split="_"))))
		#unlist(strsplit(sx[iC],split="_"))))
	#normalized copy size
	sapply(c(1:length(sx)),function(x){
		return(paste(round(
			as.numeric(unlist(strsplit(
				sx[x],split="_")))/mx,2),collapse="_"
		))
	})	
		
})

hist(as.numeric(unlist(strsplit(SizeM,split="_"))),breaks=1000)

#real copy number per gene
RC<-sapply(c(1:length(SizeM[1,])),function(x){
	sx<-SizeM[,x]
	return(sapply(c(1:length(sx)),function(x){
		return(sum(as.numeric(unlist(
			strsplit(sx[x],split="_")))))
	}))
})

#replace values outside the range 0.7-1.3 as missing data
RC <-ifelse(RC <0.7,0,ifelse(RC>1.3,0,1))

#count observed duplications as missing data
RC <-ifelse(t(ANNOc)==1, RC,0)

plot(as.numeric(STAT[,4])/1000,apply(RC,1,sum),
	log="x",pch=ifelse(data1[,7]=="Genomics",19,1),
	las=1,cex.axis=0.8,xlab="N50 (kB,log scale)",
	ylab="Core genes")

#remove genes with missing data in more than PC2% of genomes

PC2=0.7

ANNOc2<-subset(ANNOc,apply(RC,2,
	sum)>= PC2*length(data1[,1]))
	
RC2<-t(subset(t(RC),apply(RC,2,
	sum)>= PC2*length(data1[,1])))	

rownames(RC2)<-data1[,3]

plot(as.numeric(STAT[,4])/1000,apply(RC2,1,sum),
	log="x",pch=ifelse(data1[,7]=="Genomics",19,1),
	las=1,cex.axis=0.8,xlab="N50 (kB,log scale)",
	ylab="Core genes")


setwd(pG)
write.table(RC2,"Lichenihabitantaceae-core-gene-occ.txt")

#Core gene type
CGtype<-sapply(c(1:length(data1[,1])),function(x){
	print(x)	
	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T) #you have to modify the path
	as.vector(sapply(rownames(ANNOc2),function(x){
		return(c(subset(
			ANNOx[,3],ANNOx[,8]== x),NA)[1])
	}))
})

CGtype<-sapply(c(1:length(CGtype[,1])),function(x){
	return(setdiff(CGtype[x,],NA))
})

CGname<-paste("CG",c(1:length(rownames(ANNOc2))),sep="_")

CGstat<-cbind(CGname,rownames(ANNOc2), CGtype,
	as.numeric(apply(RC2,2,sum)))

colnames(CGstat)<-c("Name","Annotation","Type","Occurence")

write.table(CGstat,
	"Core-gene-code.txt",row.names=F)


pdf("CoreGenome_N50.pdf",width=6,height=6)
plot(as.numeric(STAT[,4])/1000,apply(RC2,1,sum),
	log="x",pch=ifelse(data1[,7]=="Genomics",19,1),
	las=1,cex.axis=0.8,xlab="N50 (kB,log scale)",
	ylab="Core genes")
dev.off()


#write a fasta file with core genes for each genome
dir.create("core-genes")


sapply(c(1:length(data1[,1])),function(x){
	print(x)
	tx<-data1[x,1]
	rcx= RC2[x,]
	ix<-subset(c(1:length(rcx)), rcx!=0)
	print(tx)
	###Annotation file
	ANNOx<-read.table(gsub("jean-baptisteleducq","jean-baptiste",data1[x,6]),header=T) #you have to modify the path

	SEQx<-sapply(ix,function(x){
		return(list(unlist(strsplit(subset(
			as.vector(ANNOx[,12]),
			as.vector(ANNOx[,8])== rownames(ANNOc2)[x]),
			split=""))))
	})
	setwd(paste(pG,"core-genes",sep="/"))
	write.fasta(SEQx,names= CGname[ix],
		file.out=paste(data1[x,1],".fas",sep=""))
	return("")	
})

#Write one fasta per core gene

setwd(paste(pG,"core-genes",sep="/"))

mclapply(c(1:length(CGname)),function(x){
	Gx= CGname[x]
	rcx= RC2[,x]
	ix<-subset(c(1:length(data1[,1])), rcx!=0)
	print(x)
	SEQx<-sapply(ix,function(x){
		#print(x)
		tx<-data1[x,1]
		seqx<-read.fasta(paste(tx,".fas",sep="")	)
		return(list(unlist(
			seqx[[subset(c(1:length(seqx)),
		names(seqx)==Gx)]])))
	})
	write.fasta(SEQx,names= data1[ix,1],
		file.out=paste(Gx,".fas",sep=""))
	return("")
},mc.cores=8)
########################################################################################################################################################################################
#Align each core gene (do separately for protein-encoding genes and rna genes)


####Align non protein coding genes
Irna<-subset(c(1:length(CGname)), CGtype =="rna")

setwd(paste(pG,"core-genes",sep="/"))

mclapply(Irna,function(x){	
	Ix<-CGname[x] #short gene name
	print(Ix)
	SEQw<-msa(paste(Ix,
		".fas",sep=""),method="ClustalW",type="dna",order="input")
	#Retrieve aligned sequence in seqinr format
	SEQn<-msaConvert(SEQw,type="seqinr::alignment")$seq
	#Get sequence names from the unaligned nucleotide sequence
	NX=names(read.fasta(paste(Ix,
		".fas",sep="")))
	#Convert
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(strsplit(SEQn[x],split=""))
	})
	write.fasta(SEQn, names=NX,
		file.out =paste(Ix,
		"_aligned.fas",sep=""))
		
	file.remove(paste(Ix,
		c(".dnd",".aln"),sep=""))
	
		
},mc.cores=8)


####Align protein encoding genes
Ipeg<-subset(c(1:length(CGname)), CGtype =="peg")

setwd(paste(pG,"core-genes",sep="/"))

mclapply(Ipeg,function(x){
	Ix<-CGname[x] #short gene name
	print(Ix)

	#Get the unaligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		".fas",sep=""))

	#Sequence names
	NX=names(SEQ)	
	
	#Remove gaps
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(list(subset(unlist(SEQ[[x]]),
			unlist(SEQ[[x]])!="-")))
	})
	
	#Sequence sizes
	Sx<-as.numeric(summary(SEQ)[,1])
	Sx<-c(Sx,min(Sx)-3)
	
	#Add a vector of "Ns" for potential missing data at the end and begining of each sequence (length: maximum-minimum sequence size)
	nx<-sapply(c(1:(max(Sx)-min(Sx))),function(x){
		return("N")
	})
	
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(list(c(nx,as.vector(SEQ[[x]]),nx)))
	})
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQ)

	#Write a temporary file with amino acid sequences unaligned
	write.fasta(SEQaa,names=NX,paste(Ix,"_AA.fas",sep=""))


	#Align AA sequences together using ClustalW - IMPORTANT: keep same order as input, since sequence names are removed during analysis, and retrieve sequence names in the original unaligned fasta file


	SEQw<-msa(paste(Ix,"_AA.fas",sep=""),method="ClustalW",
		type="protein",order="input")

	#Retrieve aligned sequence in seqinr format
	SEQn<-msaConvert(SEQw,type="seqinr::alignment")$seq
	
	#Convert back in nucleotide sequence

	SEQn<-sapply(c(1:length(NX)),function(x){
		#Unaligned nucleotide sequence
		seqx<-as.vector(SEQ[[x]])
		#Unaligned protein sequence
		sequ<-unlist(strsplit(SEQaa[[x]],split=""))
		#Aligned protein sequence
		seqa<-unlist(strsplit(SEQn[[x]],split=""))
		#Unaligned nucleotide sequence coordinates
		coordx<-c(1:length(seqx))
		#Unaligned proteine sequence coordinates
		coordu<-matrix(coordx,ncol=3,byrow=T)
		#Removed stop codon (because it's automaticaly removed by msa alignement)
		coordu<-subset(coordu, sequ!="*")
		#Aligned protein sequence coordinates
		coorda<-matrix(c(1:(length(seqa)*3)),
			ncol=3,byrow=T)
		#Ungapped positions in Aligned protein sequence coordinates
		coordau<-sort(as.vector(subset(coorda, 
			seqa!="-")))
		#Gapped position in Aligned protein sequence coordinates
	coordag<-sort(as.vector(subset(coorda, 
		seqa=="-")))		
		#Nucleotide sequence (ungapped position)	
		seqx<-seqx[sort(as.vector(coordu))]
		#Nucleotide sequence (Gapped position)
		seqg<-ifelse(seq(0,0,length.out=
			length(coordag))==0,"-","")
		return(list(c(seqg, seqx)
			[order(c(coordag, coordau))]))
	})
	
	#Remove position that only contain missing data
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(as.vector(unlist(	SEQn[[x]])))
	})

	SEQn<-subset(SEQn,sapply(c(1:length(SEQn[,1])),function(x){
		seqx<-SEQn[x,]
		return(length(subset(seqx, seqx=="-"))+
		length(subset(seqx, seqx=="N")))
	})!=length(NX))
	
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(list(SEQn[,x]))
	})
		
	write.fasta(SEQn, names=NX,
		file.out =paste(Ix,
		"_aligned.fas",sep=""))
		
	file.remove(paste(Ix,
		c("_AA.fas","_AA.dnd","_AA.aln"),sep=""))		
		
},mc.cores=8)


#Review each protein-coding alignment and replace "NNN(---)NNN" motifs by "NNN(NNN)NNN"

setwd(paste(pG,"core-genes",sep="/"))

mclapply(c(Ipeg,Irna),function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)

	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		"_aligned.fas",sep=""))
	
	SEQn<-sapply(c(1:length(SEQ)),function(x){
		#print(x)
		seqx<-c(as.vector(SEQ[[x]]),"x","n","-","n","x","n","-","n")
		lx=length(seqx)
		z<-c(0,sapply(c(2: lx),function(x){
			ifelse(seqx[x-1]==seqx[x],1,0)
		}))
		seqxx<-subset(seqx,z==0)
		z<-subset(c(1:length(seqx)),z==0)
		
		y<-t(sapply(c(3:length(seqxx)),function(x){
			return(c(
			paste(z[(x-2):x],collapse="_"),
			paste(seqxx[(x-2):x],collapse="_")
			))
		}))
		
		y<-subset(y,y[,2]=="n_-_n")
		
		#positions to replace by 'Ns'
		pn<-unique(unlist(sapply(c(1:length(y[,1])),function(x){
			yx<-as.numeric(unlist(strsplit(y[x,1],split="_")))
			return(c(min(yx):max(yx)))
		})))
		sn<-sapply(pn,function(x){return("n")})
		#positions to keep
		pk<-setdiff(c(1:lx),pn)
		sk<-seqx[pk]
		
		seqn<-c(sn, sk)[order(c(pn,pk))]
		return(list(toupper(seqn[1:(lx-8)])))
	})
	
	write.fasta(SEQn,paste(Ix,"_alignedN.fas",sep=""),names=names(SEQ))
	
},mc.cores=4)

#Check that gap replacement by Ns did not affected alignments

setwd(paste(pG,"core-genes",sep="/"))

#only "n" should be returned

unique(unlist(
	sapply(c(1:length(CGname)),function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence (before and after N correcting)
	SEQ<-read.fasta(paste(Ix,
		"_aligned.fas",sep=""))	
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	sumx<-unlist(sapply(c(1:length(SEQ)),function(x){
		seq<-as.vector(SEQ[[x]])
		seqn<-as.vector(SEQn[[x]])
		return(subset(seqn,seq!=seqn))
	}))
	print(c(x,unique(sumx))	)
	return(sumx)
})))


#Proportion of missing data (only for protein coding)
setwd(paste(pG,"core-genes",sep="/"))

PPX<-unlist(sapply(Ipeg,function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQn)
	
	#Convert in matrix
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(as.vector(SEQaa[[x]]))
	})
	
	#Count the proportion of X per position
	pn<-sapply(c(1:length(SEQn[,1])),function(x){
		length(subset(SEQn[x,],SEQn[x,]=="X"))
	})/length(NX)
	return(pn)
}))

#Remove aa positions with more than ppn% of missing data (only for protein-coding)
ppn=0.9

setwd(paste(pG,"core-genes",sep="/"))

#Ony "n" should be returned


sumNEW<-mclapply(Ipeg,function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQn)
		
	#Convert both in matrix
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(as.vector(SEQn[[x]]))
	})
	SEQaa <-sapply(c(1:length(SEQaa)),
		function(x){
		return(as.vector(SEQaa[[x]]))
	})
	#Count the proportion of X per position
	pn<-sapply(c(1:length(SEQaa[,1])),
		function(x){
		length(subset(SEQaa[x,], SEQaa[x,]=="X"))
	})/length(NX)
	
	pn <-as.vector(matrix(t(cbind(pn,pn,pn)),
		ncol=1,byrow=T))
	
	SEQn<-subset(SEQn, pn<(1-ppn))
	
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(list(SEQn[,x]))
	})
	
	write.fasta(SEQn,paste(Ix,
		"_alignedNF.fas",sep=""),names=NX)
	return(c(x,Ix,length(SEQn),length(pn),
		length(SEQn[[1]])))	
		
},mc.cores=4)

sumNEW <-matrix(unlist(sumNEW),ncol=5,byrow=T)

#Same with rna

sumNEWrna<-mclapply(Irna,function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	write.fasta(SEQn,paste(Ix,
		"_alignedNF.fas",sep=""),names=NX)
	return(c(x,Ix,length(SEQn),length(SEQn[[1]]),
		length(SEQn[[1]])))	
		
},mc.cores=4)

sumNEWrna <-matrix(unlist(sumNEWrna),ncol=5,byrow=T)

sumNEWrna<-cbind(sumNEWrna,"rna")
sumNEW<-cbind(sumNEW,"peg")
sumNEW<-rbind(sumNEWrna, sumNEW)

#remove intermediate files

sapply(c(Ipeg,Irna),function(x){
	
	Ix<-CGname[x] #short gene name
	print(Ix)
	file.remove(paste(Ix,
		c("_aligned.fas","_alignedN.fas"),sep=""))
})	

setwd(pG)
colnames(sumNEW)<-c("Name1","Name2","Sequences","Size1(bp)","Size2(bp)","Type")
write.table(sumNEW,"Summary-per-core-gene-after-alignment.txt",row.names=F)


###Run each core gene phylogeny separately with RAxML


###		#!/bin/bash
###		#SBATCH --mem=128GB
###		#SBATCH -N 1
###		#SBATCH --time=2-00:00:00
###		#SBATCH -p batch_72h
###		#SBATCH --ntasks-per-node=4
###		#SBATCH --array=1-721%8

###		module load StdEnv/2020  gcc/9.3.0 
###		module load raxml/8.2.12

###		raxmlHPC-MPI -T 4 -f a -s CG_${SLURM_ARRAY_TASK_ID}_alignedNF.fas -x 12345 -p 12345 -# 1000 -m GTRGAMMA -n CG_${SLURM_ARRAY_TASK_ID}-GTRGAMMA-1000


#Reopen each tree in R and format them for ASTRAL analysis

#work directory for genomes and annotation files
pG<-"/Users/jean-baptisteleducq/Dropbox/Lichenibacterium/Available-Genome-Ressources/Genomes"

setwd(pG)

data1<-read.table("data4.txt",header=T)
CGstat<-read.table("Core-gene-code.txt",header=T)
sumNEW<-read.table("Summary-per-core-gene-after-alignment.txt",header=T)

#new tree labels

LABnew<-paste(data1[,1],ifelse(data1$Type=="Genomics",paste(data1[,2],data1[,3],sep="_"),paste("Meta",data1$Source,sep="_")),sep="_")


###Convert RAxML output tree file and pool them in multiphylo file to run with ASTRAL (trees must be converted)
library(ape)
library(ips)


setwd(paste(pG,"RAxML-trees",sep="/"))

mBS=10 #minimum node support to keep node

#Raxml best trees
LX=length(CGstat[,1])

#Reformat and combine all core genes tree in a single multiphylo tree (to use as an input for ASTRAL)

TREE<-sapply(c(1: LX),function(x){
	print(x)
	#read tree as character
	treex<-as.vector(read.table(paste(
		"RAxML_bipartitionsBranchLabels.CG_",x,"-GTRGAMMA-1000",sep="")))
	#split nodes
	treex<-unlist(strsplit(
		as.character(treex),split=":"))
	
	#identify nodes with brackets (node supports)
	bi<-grep("[[]",treex)
	bo<-setdiff(c(1:length(treex)), bi)
	#extract node support and convert in readable tree format: #BranchLength[NodeSupport])# converted in #NodeSupport:BranchLength)# 
	ns<-sapply(bi,function(x){
		nx<-unlist(strsplit(treex[x],split="[[]"))
		nx<-unlist(strsplit(nx,split="[]]"))
		return(paste(nx[2],":",nx[1],nx[3],sep=""))
	})
	treex<-c(ns, treex[bo])[order(c(bi,bo))]
	
	#if last sign is not ")" or ";", add ":"
	bi<-grep("[)]", treex)
	bo<-setdiff(c(1:length(treex)), bi)
	
	ns<-paste(treex[bo],":",sep="")
	
	treex<-c(ns, treex[bi])[order(c(bo,bi))]
	
	treex<-paste(treex,collapse="")
	write.table(treex,
		paste("RAxML_tree.",x,"_reformated",sep=""),
		row.names=F,col.names=F,quote=F)
	
	treex<-read.tree(paste("RAxML_tree.",
		x,"_reformated",sep=""))
		
	#replace tree tip labels	
	TIP<-treex$tip.label
	nTIP<-as.vector(sapply(TIP,function(x){
		return(subset(LABnew ,data1[,1]==x))
	}))		
	treex$tip.label <-nTIP
	
	treex <-collapseUnsupportedEdges(treex, 
		value = "node.label", cutoff= mBS)	
	write.tree(treex,paste("RAxMLtree-CG_",
		x,"-reformated",sep=""))
	return(list(unroot(treex)))
})

class(TREE) <- "multiPhylo"

setwd(pG)

write.tree(TREE,paste("Lichenihabitantaceae-Multi-coregene-feb-2025_mBS=",mBS,".tre",sep=""))

#Then in shell, compute lineage tree with Astral (-Xmx3000M to allow the use of more memory by java)

###		cd Dropbox/Lichenibacterium/Available-Genome-Ressources/Genomes/Astral.5.7.8

###		java -Xmx3000M -jar astral.5.7.8.jar -i ../Lichenihabitantaceae-Multi-coregene-feb-2025_mBS=50.tre -o ../Lichenihabitantaceae-AstralLinTree-feb-2025_mBS=50.tre

#The consensus tree is readable in FigTree but not in MEGA or R (just topology displayed) Convert the consensus tree to be readable in MEGA

	 #NodeSupport:BranchLength)# to convert in #BranchLength[NodeSupport])#

mBS=50


setwd(pG)

########retrieve tree in text format
TREEcons<-as.vector(read.table(paste("Lichenihabitantaceae-AstralLinTree-feb-2025_mBS=",mBS,".tre",sep="")))

#cut the tree
TREEcons <-unlist(strsplit(
	as.character(TREEcons),split=":"))

TREEcons <-strsplit(TREEcons,split=")")

#permut branch support and branch length values, and put branch support in brackets
TREEconsN <-paste(sapply(c(2:(length(TREEcons)-1)),function(x){
	t1<-TREEcons[[x-1]]
	t2<-TREEcons[[x]]
	bx<-as.numeric(t1[2])*100
	lx<-strsplit(t2[1],split=",")[[1]][1]
	ox<-strsplit(t2[1],split=",")[[1]]
	ox<- ox[2:length(ox)]
	ox<-subset(ox,is.na(ox)==F)
	ox<-subset(ox,ox!= lx)
	ox<-paste(ox,collapse=",")
	return(paste(c("):",lx,"[",bx,"],",ox),collapse=""))
	
}),collapse="")

TREEconsN<-paste(TREEcons[[1]][1], TREEconsN,"):0.0);",sep="")

TREEconsN<-gsub(",)",")", TREEconsN)

write.table(TREEconsN,paste("Lichenihabitantaceae-AstralLinTree-feb-2025_mBS=",mBS,".nwk",sep=""),row.names=F,col.names=F,quote=F)

#Retrieve tip names
library(ape)

TC<-read.tree(paste("Lichenihabitantaceae-AstralLinTree-feb-2025_mBS=",mBS,".tre",sep=""))

LAB<-TC$tip.label

#Retrieve branch lengths

BL<-as.vector(read.table(paste("Lichenihabitantaceae-AstralLinTree-feb-2025_mBS=",mBS,".tre",sep="")))

#cut the tree
BL <-unlist(strsplit(
	as.character(BL),split=":"))

BL <-strsplit(BL,split=")")

BL <-sapply(c(2:(length(TREEcons)-1)),function(x){
	t1<-TREEcons[[x-1]]
	t2<-TREEcons[[x]]
	bx<-as.numeric(t1[2])*100
	lx<-strsplit(t2[1],split=",")[[1]][1]
	return(as.numeric(lx))	
})

#minimum branch length that is not null
BL <-min(setdiff(BL,0))

#generate a sequence of small branch length

BL<-seq(BL/10,BL/2,length.out=1000)

#ad random branch length to tips (smaller than the smallest length in the tree)

sapply(c(1:length(LAB)),function(x){
	TX<-as.vector(read.table(paste("Lichenihabitantaceae-AstralLinTree-feb-2025_mBS=",mBS,".nwk",sep="")))
	Lx=LAB[x]
	lx=unlist(strsplit(Lx,split="_"))[1]
	TX<-gsub(Lx,paste(Lx,":",sample(BL,size=1),sep=""),TX)
	write.table(TX,paste("Lichenihabitantaceae-AstralLinTree-feb-2025_mBS=",mBS,".nwk",sep=""),row.names=F,col.names=F,quote=F)
	return("")
})
