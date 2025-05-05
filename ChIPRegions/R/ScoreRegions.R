#'ScoreRegions
#' @description
#'Computes a regional signal score. Computing this score requires stratifying across each overlap pattern
#'and can only be computed for overlap patterns present within each state. The minimum sample size for each state can be set with MinSample.
#'The score is reported as a ratio to a chosen reference state, changing the reference state does not influence the ratio between the scores of any two states.
#' @export
#' @param UnionSequence Output of Unify function
#' @param Viterbi A list of Viterbi sequences. Each Viterbi sequence corresponds to a sequence in ChromSplit.
#' @param Bedlist A list of GRanges objects used to create UnionSequence
#' @param Scorename Character value, the name of the metadata column to be compared across regions
#' @param ReferenceState Numeric value, the HMM state to be used as the baseline
#' @param MinSample Numeric value, overlap patterns will removed from the comparison if any state contains less than MinSample observations

ScoreRegions=function(UnionSequence,Viterbi,Bedlist,Scorename="SignalValue",ReferenceState=1,MinSample=20){

  Viterbi=unlist(Viterbi)
  NStates=length(unique(unlist(Viterbi)))
  States=sort((unique(Viterbi)))
  if(NStates<2){
    print("Only One State")
  }
  NObservations=length(unique(UnionSequence$OverlapPattern))
  OverlapSymbols=unique(UnionSequence$OverlapPattern)
  ScoreMatrix=matrix(1,nrow=NObservations,ncol=NStates)
  colnames(ScoreMatrix)=sort((unique(Viterbi)))
  PreCount=as.matrix(table(UnionSequence$OverlapPattern,Viterbi))
  rownames(ScoreMatrix)=rownames(PreCount)
  ScoreMatrix[PreCount<MinSample]=NA

  ScoreMatrix=na.omit(ScoreMatrix)

  Remove=rownames(PreCount)[!(rownames(PreCount)%in%rownames(ScoreMatrix))]
  if(length(Remove)==length(rownames(PreCount))){
    print("Broken")
    return("Broken")
  }
  for(i in 1:length(Bedlist)){
    Mapp=IRanges::findOverlaps(Bedlist[[i]],UnionSequence)
    RepSignal=unname(unlist(Bedlist[[i]]@elementMetadata@listData[Scorename]))
    DataFra=data.frame(RepSignal[Mapp@from],Viterbi[Mapp@to],UnionSequence$OverlapPattern[Mapp@to])
    for(j in States){
      StateFra=DataFra[DataFra$Viterbi.Mapp.to.==j,]
      OverlapScore=aggregate(StateFra$RepSignal.Mapp.from.,FUN="mean",by=list(StateFra$UnionSequence.OverlapPattern.Mapp.to.))
      OverlapScore=OverlapScore[!OverlapScore$Group.1%in%Remove,]
      ScoreMatrix[OverlapScore$Group.1,j]=ScoreMatrix[OverlapScore$Group.1,j]*OverlapScore$x
    }
    ScoreMatrix=ScoreMatrix/ScoreMatrix[,ReferenceState]
  }
  ScoreMatrix=ScoreMatrix^(1/rowSums(read.table(text=rownames(ScoreMatrix))))
  return(ScoreMatrix)

}
