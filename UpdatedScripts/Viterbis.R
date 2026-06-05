SetsMeta=read.csv("SubSampSets.csv")
AllMeta=read.csv("ENCODEExperimentsMeta.csv")
TreatsAndURLs=read.csv("TreatsAndURLs.csv")


library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
ErrorCatcherReps=c()
Signal=list()
NEWPUS=list()
Model=readRDS("MACOptimalModelExpSets.rds")
Unions=readRDS("MACUnions.rds")
RepModels=readRDS("MACRepSetHMMs.rds")
RepPUS=readRDS("MACUnionsRepSet.rds")


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

RepViterbis=list()
for(i in 1:length(RepModels)){
  if(length(RepModels[[i]])<19){
    next
  }
  Chrom=ChromGet(RepPUS[[i]])
  RepViterbis[[i]]=viterbi(Chrom[["ChromSplit"]],RepModels[[i]])
  names(RepViterbis)[i]=names(RepModels)[i]

}

MACViterbis=list()
for(i in 1:length(Model)){
  if(length(Model[[i]])<19){
    next
  }
  print(i)
  Chrom=Unions[[i]]
  MACViterbis[[i]]=viterbi(Chrom[["ChromSplit"]],Model[[i]])
  names(MACViterbis)[i]=names(Model)[i]

}
