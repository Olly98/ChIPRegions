SetsMeta=read.csv("SubSampSets.csv")
AllMeta=read.csv("ENCODEExperimentsMeta.csv")
TreatsAndURLs=read.csv("TreatsAndURLs.csv")
library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
GEMOptimalState=list()
GEMOptimalModel=list()
MACOptimalState=list()
MACOptimalModel=list()
SISSOptimalState=list()
SISSOptimalModel=list()
ErrorCatcherMACS=c()
ErrorCatcherSISS=c()
ErrorCatcherGEM=c()
ErrorCatcherGEMRep=c()
GEMOptimalModelRep=c()
MACOptimalModelRep=c()
ErrorCatcherMACRep=c()
ErrorCatcherSISSRep=c()
SISSOptimalModelRep=c()
Issue=c()

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



for(i in 1:nrow(SetsMeta)){


SubMeta=AllMeta[AllMeta$Experiment.target==SetsMeta$UnfoldedTFs[i],]

SubMeta=SubMeta[SubMeta$Biosample.term.name==SetsMeta$UnfoldedCells[i],]
Sperries=unique(SubMeta$Experiment.accession)
Accs=c()
for(j in 1:length(Sperries)){
Accs=c(Accs,substr(TreatsAndURLs$TreatURLDownload[TreatsAndURLs$TreatExpAccessions==Sperries[j]],start=31,stop=41))

}
MACSFiles=substr(list.files("./MACSData/"),1,11)
AllInMACS=Accs%in%MACSFiles
GEMFiles=list.files("./GEMR/")
AllInGEM=Accs%in%GEMFiles
SISSFiles=substr(list.files("./SISS/"),1,11)
AllInSISS=Accs%in%SISSFiles
AllIn=c(AllInMACS,AllInGEM,AllInSISS)
if(sum(AllIn==FALSE)>0){
Issue[i]="Here"
next
}
GEMNames=c()
BedList=list()
for(k in 1:length(Accs)){
GEMNames[k]=paste("./GEMR/",Accs[k],"/",paste(Accs[k],".GPS_events.txt",sep=""),sep="")
file=read.table(GEMNames[k],header=TRUE)
chr=unlist(strsplit(file$Position,split=":"))[c(T,F)]
Coords=unlist(strsplit(file$Position,split=":"))[c(F,T)]
Start=as.numeric(Coords)-250
Stop=as.numeric(Coords)+250
NewFrame=data.frame(chr,Start,Stop)
BedList[[k]]=makeGRangesFromDataFrame(NewFrame)
}
print(lengths(BedList))



PUS=Unify(BedList)
print(length(PUS))
Chrom=ChromGet(PUS,MinPeaks=3)
print(names(Chrom))
StoreModel=list()
BICs=c()

BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE



Model=tryCatch({hmm(Chrom[["ChromSplit"]],K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
ErrorCatcherGEM[i]=sum(is.na(BICs))
BICs[is.na(BICs)]=Inf

GEMOptimalState[i]=which(BICs==min(BICs))
GEMOptimalModel=c(GEMOptimalModel,setNames(list(StoreModel[[which(BICs==min(BICs))]]),paste(SetsMeta$UnfoldedTFs[i],SetsMeta$UnfoldedCells[i])))

for(k in 1:length(Accs)){
Name=paste("./MACSData/",Accs[k],"_peaks.narrowPeak",sep="")
BedList[[k]]=import(Name,format="narrowPeak")
}

PUS=Unify(BedList)
print(length(PUS))
Chrom=ChromGet(PUS,MinPeaks=3)
print(names(Chrom))
StoreModel=list()
BICs=c()

BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE



Model=tryCatch({hmm(Chrom[["ChromSplit"]],K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
ErrorCatcherMACS[i]=sum(is.na(BICs))
BICs[is.na(BICs)]=Inf




MACOptimalState[i]=which(BICs==min(BICs))
MACOptimalModel=c(MACOptimalModel,setNames(list(StoreModel[[which(BICs==min(BICs))]]),paste(SetsMeta$UnfoldedTFs[i],SetsMeta$UnfoldedCells[i])))

for(k in 1:length(Accs)){
Name=paste("./SISS/",Accs[k],".bed",sep="")
file=read.table(Name,skip=57,header=T,nrow=length(readLines(Name))-59)
file=file[,1:3]
names(file)=c("chr","start","stop")
BedList[[k]]=makeGRangesFromDataFrame(file)
}



PUS=Unify(BedList)
print(length(PUS))
Chrom=ChromGet(PUS,MinPeaks=3)
print(names(Chrom))
StoreModel=list()


BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE


Model=tryCatch({hmm(Chrom[["ChromSplit"]],K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
ErrorCatcherSISS[i]=sum(is.na(BICs))
BICs[is.na(BICs)]=Inf


SISSOptimalState[i]=which(BICs==min(BICs))
SISSOptimalModel=c(SISSOptimalModel,setNames(list(StoreModel[[which(BICs==min(BICs))]]),paste(SetsMeta$UnfoldedTFs[i],SetsMeta$UnfoldedCells[i])))

}
saveRDS(GEMOptimalModel,file="GEMOptimalModelExpSets.rds")
saveRDS(GEMOptimalState,file="GEMOptimalStateExpSets.rds")
saveRDS(MACOptimalModel,file="MACOptimalModelExpSets.rds")
saveRDS(MACOptimalState,file="MACOptimalStateExpSets.rds")
saveRDS(SISSOptimalModel,file="SISSOptimalModelExpSets.rds")
saveRDS(SISSOptimalState,file="SISSOptimalStateExpSets.rds")
saveRDS(ErrorCatcherMACS,file="ErrorMACS.rds")
saveRDS(ErrorCatcherGEM,file="ErrorGEM.rds")
saveRDS(ErrorCatcherSISS,file="ErrorSISS.rds")
saveRDS(Issue,file="IssueFileMissing.rds")


