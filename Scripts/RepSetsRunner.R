SetsMeta=read.csv("ExpSets.csv")
AllMeta=read.csv("ENCODEExperimentsMeta.csv")
library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
OptimalStateReps=list()
OptimalModelReps=list()
ErrorCatcher=c()
ChromGet=function (UnifiedSequence, MaxGap = NULL, MinPeaks = 3)
{
  ChromsToRemove = as.vector((UnifiedSequence@seqnames@values)[which(UnifiedSequence@seqnames@lengths <
    MinPeaks)])
  if (length(ChromsToRemove) > 0) {
    UnifiedSequence = UnifiedSequence[!UnifiedSequence@seqnames %in%
      ChromsToRemove]
    print(paste(as.character(ChromsToRemove), "contained less than",
      as.character(MinPeaks), "peaks, it has been removed"))
  }
 if (is.null(MaxGap)) {
    UnifiedSequence = GenomeInfoDb::sortSeqlevels((UnifiedSequence))
    UnifiedSequence = sort(UnifiedSequence)
    ChromLengths = as.vector(UnifiedSequence@seqnames@lengths)
    ChromNames = as.vector(UnifiedSequence@seqnames@values)
    ChromSeqs = split(UnifiedSequence$OverlapPattern, rep(seq_along(ChromLengths),
      ChromLengths))
    names(ChromSeqs) = ChromNames
    Output = list(UnifiedSequence, ChromSeqs)
    names(Output) = c("Unified Sequence", "ChromSplit")
    return(Output)
  }}

Unify=function (Bedlist)
{
  PUS = Signac::UnifyPeaks(Bedlist, mode = "reduce")
  PUS$OverlapPattern = do.call(paste, lapply(Bedlist, IRanges::overlapsAny,
    query = PUS))
  return(PUS)
}

Uns=unique(AllMeta$Experiment.accession)

for(i in 1:length(Uns)){

print(i)
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
OptimalStateReps[i]=which(BICs==min(BICs))
OptimalModelReps[[i]]=StoreModel[[which(BICs==min(BICs))]]
system("rm *bed.gz")

}

saveRDS(OptimalModelReps,file="OptimalHMMsRepSetsRand.rds")
saveRDS(OptimalStateReps,file="OptimalStateRepSetsRand.rds")


