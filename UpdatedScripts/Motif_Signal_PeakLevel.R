AllMeta=read.csv("ENCODEExperimentsMeta.csv")
SetsMeta=read.csv("SubSampSets.csv")
TreatsAndURLs=read.csv("TreatsAndURLs.csv")
library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
Issue=c()
library("universalmotif")
library(memes)
library(BSgenome.Hsapiens.UCSC.hg38)
library("universalmotif")
library(memes)
library(BSgenome.Hsapiens.UCSC.hg38)
Chroms=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
hs.genome=BSgenome.Hsapiens.UCSC.hg38
AllMotifs=paste("./JASPAR2024_CORE_redundant_pfms_meme/",list.files("./JASPAR2024_CORE_redundant_pfms_meme"),sep="")
AllPWMs=lapply(AllMotifs,read_meme)
GetTFNames=c()
for(i in 1:length(AllPWMs)){
  GetTFNames[i]=AllPWMs[[i]]@altname
}
names(AllPWMs)=GetTFNames
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

MedianDiff=c()
for(i in 1:nrow(SetsMeta)){
  print(i)
  SubMeta=AllMeta[AllMeta$Experiment.target==SetsMeta$UnfoldedTFs[i],]
  SubMeta=SubMeta[SubMeta$Biosample.term.name==SetsMeta$UnfoldedCells[i],]
  TF=unlist(strsplit(SetsMeta$UnfoldedTFs[i],split="-"))[c(T,F)]
  TF=TF[!(TF%in%c("7","S3"))]
  if(!(TF[1]%in%GetTFNames)){
    next
  }
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
  BedList=list()
  for(k in 1:length(Accs)){
    Name=paste("./MACSData/",Accs[k],"_peaks.narrowPeak",sep="")
    BedList[[k]]=import(Name,format="narrowPeak")
  }
  SigDiff=rep(NA,length(BedList))
  for(j in 1:length(BedList)){
    Chrum=Chroms[Chroms%in%seqnames(BedList[[j]])]
    TrimBed=keepSeqlevels(BedList[[j]],Chrum,pruning.mode = "coarse")
    Seq=get_sequence(TrimBed,hs.genome)
    StateMo=runFimo(Seq,motifs = AllPWMs[names(AllPWMs)==TF[1]],thresh=1e-05)
    if(is.null(StateMo)){
      next
    }
    SigMo=mean((TrimBed$signalValue)[overlapsAny(TrimBed,StateMo)])
    SigNoMo=mean((TrimBed$signalValue)[!overlapsAny(TrimBed,StateMo)])
    SigDiff[j]=SigMo/SigNoMo
  }
  MedianDiff[i]=median(SigDiff)
}

