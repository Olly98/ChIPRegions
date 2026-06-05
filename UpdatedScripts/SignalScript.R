SetsMeta=read.csv("SubSampSets.csv")
AllMeta=read.csv("ENCODEExperimentsMeta.csv")
TreatsAndURLs=read.csv("TreatsAndURLs.csv")


library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)

ErrorCatcherReps=c()
Signal=list()
RepSignal=list()
RepSignals=list()


Model=readRDS("MACOptimalModelExpSets.rds")
Unions=readRDS("MACUnions.rds")
RepModels=readRDS("MACRepSetHMMs.rds")
RepPUS=readRDS("MACUnionsRepSet.rds")
RepViterbis=readRDS("/Users/oliverhughes/Downloads/DataFolder/RepSetViterbis.rds")
MACViterbis=readRDS("/Users/oliverhughes/Downloads/DataFolder/MACViterbis.rds")
NEWPUS=list()

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


ScoreRegions=function (UnionSequence, Viterbi, Bedlist, Scorename = "signalValue",
                       ReferenceState = 1, MinSample = 20)
{
  #Viterbi = unlist(Viterbi)
  NStates = length(unique(unlist(Viterbi)))
  States = sort((unique(Viterbi)))
  if (NStates < 2) {
    print("Only One State")
  }
  NObservations = length(unique(UnionSequence$OverlapPattern))
  OverlapSymbols = unique(UnionSequence$OverlapPattern)
  ScoreMatrix = matrix(1, nrow = NObservations, ncol = NStates)
  colnames(ScoreMatrix) = sort((unique(Viterbi)))
  PreCount = as.matrix(table(UnionSequence$OverlapPattern,
                             Viterbi))
  rownames(ScoreMatrix) = rownames(PreCount)
  ScoreMatrix[PreCount < MinSample] = NA
  ScoreMatrix = na.omit(ScoreMatrix)
  Remove = rownames(PreCount)[!(rownames(PreCount) %in% rownames(ScoreMatrix))]
  if (length(Remove) == length(rownames(PreCount))) {
    print("Broken")
    return("Broken")
  }
  for (i in 1:length(Bedlist)) {
    Mapp = IRanges::findOverlaps(Bedlist[[i]], UnionSequence)
    RepSignal = unname(unlist(Bedlist[[i]]@elementMetadata@listData[Scorename]))
    DataFra = data.frame(RepSignal[Mapp@from], Viterbi[Mapp@to],
                         UnionSequence$OverlapPattern[Mapp@to])
    if(length(unique(DataFra$Viterbi.Mapp.to))!=NStates){
      print("Not all replicates are represented in each state")
      return("Not all replicates are represented in each state")

    }
    for (j in States) {
      StateFra = DataFra[DataFra$Viterbi.Mapp.to. == j,
      ]
      OverlapScore = aggregate(StateFra$RepSignal.Mapp.from.,
                               FUN = "mean", by = list(StateFra$UnionSequence.OverlapPattern.Mapp.to.))
      OverlapScore = OverlapScore[!OverlapScore$Group.1 %in%
                                    Remove, ]
      ScoreMatrix[OverlapScore$Group.1, j] = ScoreMatrix[OverlapScore$Group.1,
                                                         j] * OverlapScore$x
    }
    ScoreMatrix = ScoreMatrix/ScoreMatrix[, ReferenceState]
  }
  ScoreMatrix = ScoreMatrix^(1/rowSums(read.table(text = rownames(ScoreMatrix))))
  return(ScoreMatrix)
}










Unify=function (Bedlist)
{
  PUS = Signac::UnifyPeaks(Bedlist, mode = "reduce")
  PUS$OverlapPattern = do.call(paste, lapply(Bedlist, IRanges::overlapsAny,
                                             query = PUS))
  return(PUS)
}


TFs=unlist(strsplit(names(Model),split=" "))[c(T,F)]
Cells=unlist(strsplit(names(Model),split=" "))[c(F,T)]
for(i in 1:length(Model)){
  SubMeta=AllMeta[AllMeta$Experiment.target==TFs[i],]
  print(i)
  SubMeta=SubMeta[SubMeta$Biosample.term.name==Cells[i],]
  Sperries=unique(SubMeta$Experiment.accession)
  Accs=c()
  for(j in 1:length(Sperries)){
    Accs=c(Accs,substr(TreatsAndURLs$TreatURLDownload[TreatsAndURLs$TreatExpAccessions==Sperries[j]],start=31,stop=41))




  }
  MACList=list()
  for(k in 1:length(Accs)){
    Name=paste("./MACSData/",Accs[k],"_peaks.narrowPeak",sep="")
    MACList[[k]]=import(Name,format="narrowPeak")
  }

  HMM=(Model[[i]])
  if(length(HMM)==19){


    PUD=Unify(MACList)
    Chrom=ChromGet(PUD)


    Viter=MACViterbis[[i]]
    NEWPUS=c(NEWPUS,setNames(list(PUD),names(Model)[i]))
    Sig=ScoreRegions(Chrom[["Unified Sequence"]],unlist(Viter),MACList)
    Signal=c(Signal,setNames(list(Sig),paste(names(Model)[i])))
  }
  for(m in 1:length(Sperries)){
    print(m)
    MACListExp=MACList[which(SubMeta$Experiment.accession==Sperries[m])]
    RepModel=RepModels[[grep(Sperries[m],names(RepModels))]]
    if(length(RepModel)<19){
      next}
    PUS=RepPUS[[grep(Sperries[m],names(RepPUS))]]
    Chrom=ChromGet(PUS)
    RepVit=RepViterbis[[grep(Sperries[m],names(RepViterbis))]]

    RepSignal=ScoreRegions(Chrom[["Unified Sequence"]],unlist(RepVit),MACListExp)
    RepSignals=c(RepSignals,setNames(list(RepSignal),paste(names(Model)[i],Sperries[m],sep=" ")))
  }

}




#saveRDS(NEWPUS,file="NewPUS1.rds")
#saveRDS(RepViterbis,file="RepViterbis1.rds")
#saveRDS(RepSignals,file="RepSignal1.rds")
#saveRDS(Viterbis,file="Viterbis1.rds")
#saveRDS(Signal,file="SignalDataCH1.rds")
###SAVESAVE
