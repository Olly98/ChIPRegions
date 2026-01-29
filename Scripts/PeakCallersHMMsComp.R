SetsMeta=read.csv("SubSampSets.csv")
AllMeta=read.csv("ENCODEExperimentsMeta.csv")
TreatsAndURLs=read.csv("TreatsAndURLs.csv")
library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
ErrorCatcherReps=c()
RepSetsModel=list()
MACShareOptimalState=c()
MACShareOptimalModel=list()
GEMShareOptimalState=c()
GEMShareOptimalModel=list()
SISSShareOptimalState=c()
SISSShareOptimalModel=list()
PUSKeepMAC=list()
PUSKeepGEM=list()
PUSKeepSISS=list()
RepPUSKeep=list()

ErrorCatcherMAC=c()
ErrorCatcherGEM=c()
ErrorCatcherSISS=c()

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
GEMList=list()
for(k in 1:length(Accs)){
GEMNames[k]=paste("./GEMR/",Accs[k],"/",paste(Accs[k],".GPS_events.txt",sep=""),sep="")
file=read.table(GEMNames[k],header=TRUE)
chr=unlist(strsplit(file$Position,split=":"))[c(T,F)]
chr=paste("chr",chr,sep="")
Coords=unlist(strsplit(file$Position,split=":"))[c(F,T)]
Start=as.numeric(Coords)-250
Stop=as.numeric(Coords)+250
NewFrame=data.frame(chr,Start,Stop)
GEMList[[k]]=makeGRangesFromDataFrame(NewFrame)
}



GEMPUS=Unify(GEMList)

MACList=list()
for(k in 1:length(Accs)){
Name=paste("./MACSData/",Accs[k],"_peaks.narrowPeak",sep="")
MACList[[k]]=import(Name,format="narrowPeak")
}

SISSList=list()
for(k in 1:length(Accs)){
Name=paste("./SISS/",Accs[k],".bed",sep="")
file=read.table(Name,skip=57,header=T,nrow=length(readLines(Name))-59)
file=file[,1:3]
names(file)=c("chr","start","stop")
SISSList[[k]]=makeGRangesFromDataFrame(file)
}


SISSPUS=Unify(SISSList)





MACPUS=Unify(MACList)
GEMPUS2=GEMPUS[overlapsAny(GEMPUS,MACPUS)]
GEMPUS2=GEMPUS2[overlapsAny(GEMPUS2,SISSPUS)]
MACPUSG=MACPUS[overlapsAny(MACPUS,GEMPUS)]
MACPUS2=MACPUSG[overlapsAny(MACPUSG,SISSPUS)]
SISSPUS2=SISSPUS[overlapsAny(SISSPUS,MACPUS2)]


MACChrom=ChromGet(MACPUS2,MinPeaks=3)
PUSKeepMAC=c(PUSKeepMAC,setNames(list(MACPUS2),SetsMeta$Combs[i] ))
PUSKeepSISS=c(PUSKeepSISS,setNames(list(SISSPUS2),SetsMeta$Combs[i] ))
PUSKeepGEM=c(PUSKeepGEM,setNames(list(GEMPUS2),SetsMeta$Combs[i] ))
BICs=c()
StoreModel=list()
BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE


Model=tryCatch({hmm(MACChrom[["ChromSplit"]],K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}

ErrorCatcherMAC[i]=sum(is.na(BICs))
BICs[is.na(BICs)]=Inf

MACShareOptimalState[i]=which(BICs==min(BICs))
MACShareOptimalModel=c(MACShareOptimalModel,setNames(list(StoreModel[[which(BICs==min(BICs))]]),SetsMeta$Combs[i] ))


GEMChrom=ChromGet(GEMPUS2,MinPeaks=3)

BICs=c()
StoreModel=list()
BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

Model=tryCatch({hmm(GEMChrom[["ChromSplit"]],K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
ErrorCatcherGEM[i]=sum(is.na(BICs))
BICs[is.na(BICs)]=Inf

GEMShareOptimalState[i]=which(BICs==min(BICs))
GEMShareOptimalModel=c(GEMShareOptimalModel,setNames(list(StoreModel[[which(BICs==min(BICs))]]),SetsMeta$Combs[i] ))


SISSList=list()
for(k in 1:length(Accs)){
Name=paste("./SISS/",Accs[k],".bed",sep="")
file=read.table(Name,skip=57,header=T,nrow=length(readLines(Name))-59)
file=file[,1:3]
names(file)=c("chr","start","stop")
SISSList[[k]]=makeGRangesFromDataFrame(file)
}

SISSPUS=SISSPUS[overlapsAny(SISSPUS,MACPUS2)]

SISSChrom=ChromGet(SISSPUS,MinPeaks=3)
print(names(SISSChrom))
StoreModel=list()

BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

Model=tryCatch({hmm(SISSChrom[["ChromSplit"]],K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
ErrorCatcherSISS[i]=sum(is.na(BICs))
SISSShareOptimalState[i]=which(BICs==min(BICs))
SISSShareOptimalModel=c(SISSShareOptimalModel,setNames(list(StoreModel[[which(BICs==min(BICs))]]),SetsMeta$Combs[i] ))


for(m in 1:length(Sperries)){
GEMListExp=GEMList[which(SubMeta$Experiment.accession==Sperries[m])]
MACListExp=MACList[which(SubMeta$Experiment.accession==Sperries[m])]
SISSListExp=SISSList[which(SubMeta$Experiment.accession==Sperries[m])]

GEMExp=Unify(GEMListExp)
MACExp=Unify(MACListExp)
SISSExp=Unify(SISSListExp)

OneAll=MACExp[overlapsAny(MACExp,GEMExp)]
OneAll=OneAll[overlapsAny(OneAll,SISSExp)]
RepPUSKeep=c(RepPUSKeep,setNames(list(OneAll),paste(SetsMeta$Combs[i],Sperries[m],sep=" ")))
ChromExp=ChromGet(OneAll)

BICs=c()
StoreModel=list()
BICs=rep(Inf,length=6)
BICs[1]=1
for(k in 1:6){
print(k)
 skip_to_next <- FALSE
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

Model=tryCatch({hmm(ChromExp[["ChromSplit"]],K=k,itmax=1000)} , error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }

BICs[k]=Model$BIC
StoreModel[[k]]=Model
}
ErrorCatcherReps=c(ErrorCatcherReps,sum(is.na(BICs)))
BICs[is.na(BICs)]=Inf

RepSetsModel=c(RepSetsModel,setNames(list(StoreModel[[which(BICs==min(BICs))]]),paste(SetsMeta$Combs[i],Sperries[m],sep=" ")))

}


}
saveRDS(Issue,file="IssueMissing.rds")

saveRDS(RepSetsModel,file="RepSetHMMsShare.rds")
saveRDS(GEMShareOptimalModel,file="GEMShareOptimalModel.rds")
saveRDS(MACShareOptimalModel,file="MACShareOptimalModel.rds")
saveRDS(SISSShareOptimalModel,file="SISSShareOptimalModel.rds")
###SAVESAVE

