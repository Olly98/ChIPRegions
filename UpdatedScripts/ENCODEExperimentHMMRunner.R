SetsMeta=read.csv("ExpSets.csv")
AllMeta=read.csv("ENCODEExperimentsMeta.csv")
library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
OptimalModel=list()
OptimalState=c()


for(i in 1:nrow(SetsMeta)){

SubMeta=AllMeta[AllMeta$Experiment.target==SetsMeta$UnfoldedTFs[i],]

SubMeta=SubMeta[SubMeta$Biosample.term.name==SetsMeta$UnfoldedCells[i],]
print(length(SubMeta$File.download.URL))
FileList=list()
FileLengths=c()
for(j in 1:length(SubMeta$File.download.URL)){
BamDownloadCall=paste("curl -L -O",SubMeta$File.download.URL[j],sep=" ")
print(BamDownloadCall)
system(BamDownloadCall)
FileList[[j]]=import(paste(SubMeta$File.accession[j],".bed.gz",sep=""),format="narrowPeak")
FileLengths[j]=length(FileList[[j]])
}

PUS=Unify(FileList)
print(length(PUS))
StoreModel=list()
BICs=c()
for(k in 1:6){
Model=hmm(PUS$OverlapPattern,K=k,itmax=1000)
print(k)
BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
OptimalState[i]=which(BICs==min(BICs))
OptimalModel[[i]]=StoreModel[[which(BICs==min(BICs))]]
system("rm *bed.gz")
}
saveRDS(OptimalModel,file="OptimalHMMsExpSets.rds")
saveRDS(OptimalState,file="OptimalStateExpSets.rds")


Subs=read.csv("/Users/oliverhughes/Downloads/DataFolder/SubSampSets.csv")
Metsa=read_tsv("/Users/oliverhughes/Downloads/DataFolder/ENCODEMeta.tsv")
Metsa$Combs=paste(Metsa$`Experiment target`,Metsa$`Biosample term name`,sep=" ")
OptHMMs=readRDS("/Users/oliverhughes/Downloads/DataFolder/OptimalHMMsExpSets.rds")
ExpSets=read.csv("/Users/oliverhughes/Downloads/DataFolder/ExpSets.csv")

Names=c("Michael","Richard","John","Kevin","Xiang","Peggy","Bradley","Vishwanath","Tim")
Not=c("ENCODE Consortium Analysis Working Group","ENCODE Processing Pipeline","ENCODE Consortium")
AllExp=c()
Count=c()
for(i in 1:nrow(ExpSets)){
  Frem=Metsa[Metsa$Combs==ExpSets$Combs[i],]
  Labbos=Frem$Lab[!Frem$Lab%in%Not]
  AllExp[i]=paste(unique(Labbos),collapse="/")
  Count[i]=length(unique(Labbos))
}

Bias=lengths(OptimalHMMsExpSets)>16


LabTable=matrix(c(sum(Bias[Count>1]),sum(!Bias[Count>1]),sum(Bias[Count==1]),sum(!Bias[Count==1])),nrow=2)
chisq.test(LabTable)
