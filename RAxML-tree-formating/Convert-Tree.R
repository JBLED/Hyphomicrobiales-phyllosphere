#Script to convert RAXML output tree to be readable in FigTree

library(ape)
library(ips)

MyTree="RAxML_bipartitionsBranchLabels.SEQ16s_genomics-aligned-GTRGAMMA-1000"
MyTree="RAxML_bipartitionsBranchLabels.BLAST-best_aligned-simplified-GTRGAMMA-1000"
MyTree="RAxML_bipartitionsBranchLabels.rpoB_genomics_aligned-1000"
MyTree="RAxML_bipartitionsBranchLabels.ASV+REF-rpoB.fas-1000"
MyFolder="/Users/jean-baptisteleducq/Dropbox/Lichenibacterium/Data-16s"
MyFolder="/Users/jean-baptisteleducq/Dropbox/Lichenibacterium/Available-Genome-Ressources/Genomes"


#Optional: collapse nodes that are poorly supported (below give the minimum edge value to keep)
MinE=60
	
#read tree as character
setwd(MyFolder)
treex<-as.vector(read.table(MyTree))
#split nodes
treex<-unlist(strsplit(as.character(treex),split=":"))
	
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
	paste(MyTree,"_reformated",sep=""),
	row.names=F,col.names=F,quote=F)

#Optional: collapse unsupported edges
treex<-read.tree(paste(MyTree,"_reformated",sep=""))

treex <-collapseUnsupportedEdges(treex, 
	value = "node.label", cutoff= MinE)	
write.tree(treex,paste(MyTree,"_reformated_MinE=",MinE,sep=""))