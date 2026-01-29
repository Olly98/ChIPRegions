SetsMeta=read.csv("ExpSets.csv")
AllMeta=read.csv("ENCODEExperimentsMeta.csv")
library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
OptimalModel=list()
OptimalState=c()
ErrorCatcher=c()

Uns=unique(AllMeta$Experiment.accession)

for(i in 1:length(Uns)){

SubMeta=AllMeta[AllMeta$Experiment.accession==Uns[i],]


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
if(length(PUS)<1000){
next
}


#Chrom=ChromGet(PUS,MinPeaks=3)
#print(names(Chrom))
StoreModel=list()
BICs=c()


StoreModel <- vector(mode = "list", length = 6)
BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE



Model=tryCatch({hmm(sample(PUS$OverlapPattern),K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
ErrorCatcher[i]=sum(is.na(BICs))
BICs[is.na(BICs)]=Inf
OptimalState[i]=which(BICs==min(BICs))
OptimalModel[[i]]=StoreModel[[which(BICs==min(BICs))]]
system("rm *gz")

}
saveRDS(OptimalModel,file="OptimalHMMsRepSetsRand.rds")
saveRDS(OptimalState,file="OptimalStateRepSetsRand.rds")
