SetsMeta=read.csv("./SubSampSets.csv")
AllMeta=read.csv("./ENCODEExperimentsMeta.csv")
TreatsAndURLs=read.csv("./TreatsAndURLs.csv")


library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
ErrorCatcherReps=c()
Signal=list()
NEWPUS=list()
Model=readRDS("./MACRepSetHMMs.rds")
Unions=readRDS("./MACUnionsRepSet.rds")



Regions=list()
Issue=c()

K562=read.table("K562LiftOver.bed")
colnames(K562)[1:3]=c("chr","start","stop")
K562=makeGRangesFromDataFrame(K562,keep.extra.columns = T)

HepG2=read.table("HepG2LiftOver.bed")
colnames(HepG2)[1:3]=c("chr","start","stop")
HepG2=makeGRangesFromDataFrame(HepG2,keep.extra.columns = T)

GM12878=read.table("GM12878LiftOver.bed")
colnames(GM12878)[1:3]=c("chr","start","stop")
GM12878=makeGRangesFromDataFrame(GM12878,keep.extra.columns = T)

Refs=list(K562,HepG2,GM12878)
names(Refs)=c("K562","HepG2","GM12878")






ChromGet=function (UnifiedSequence, MaxGap = NULL, MinPeaks = 2)
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
  Viterbi = unlist(Viterbi)
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


TFs=unlist(strsplit(names(Model),split=" "))[c(T,F,F)]
Cells=unlist(strsplit(names(Model),split=" "))[c(F,T,F)]
Accession=unlist(strsplit(names(Model),split=" "))[c(F,F,T)]

RepSVTest=list()
RepSigTest=list()
for(i in 1:length(Model)){
  if(!Cells[i]%in%c("GM12878","K562","HepG2")){
    next
  }
  SubMeta=AllMeta[AllMeta$Experiment.target==TFs[i],]
  print(i)
  SubMeta=SubMeta[SubMeta$Biosample.term.name==Cells[i],]
  Sperries=unique(SubMeta$Experiment.accession)
  Accs=substr(TreatsAndURLs$TreatURLDownload[TreatsAndURLs$TreatExpAccessions==Accession[i]],start=31,stop=41)





  MACList=list()
  for(k in 1:length(Accs)){
    Name=paste("./MACSData/",Accs[k],"_peaks.narrowPeak",sep="")
    MACList[[k]]=import(Name,format="narrowPeak")
  }

  HMM=(Model[[i]])
  if(length(HMM)==19){


    PUD=Unify(MACList)
    Chrom=ChromGet(PUD,MinPeaks = 3)
    Viter=viterbi(Chrom[["ChromSplit"]],HMM)
    Slam=Chrom[["Unified Sequence"]]
    Slam$State=unlist(Viter)
    Furf=findOverlaps(Slam,Refs[[Cells[i]]])
    Frame=data.frame(Slam$State[Furf@from],Refs[[Cells[i]]]$V4[Furf@to])
    colnames(Frame)=c("State","SV")
    if(length(unique(Frame$SV))<2|length(unique(Frame$State))<2){
      RepSVTest=c(RepSVTest,setNames(list(unique(Frame$SV)),paste(names(Model)[i] ) ))
      RepSigTest=c(RepSigTest,setNames(list(unique(Frame$SV)),paste(names(Model)[i] ) ))
    }else{




      colnames(Frame)=c("State","SV")
      print(nrow(Frame))

      print(Cells[i])
      test=chisq.test(Frame$State,Frame$SV)

      tempt=c()
      for(k in 1:10000){
        gem=Rle(unlist(Viter))
        mere=data.frame(gem@values,gem@lengths)
        inds=sample(nrow(mere))
        mere=mere[inds,]
        gem@values=mere$gem.values
        gem@lengths=mere$gem.lengths
        NewVit=rep(gem@values,gem@lengths)
        Mlap=Refs[[Cells[i]]]
        Flap=findOverlaps(Slam,Mlap)
        frame=data.frame(NewVit[Flap@from],Mlap$V4[Flap@to])
        if(length(unique(frame$NewVit.Flap.from.))<2|length(unique(frame$Mlap.V4.Flap.to.))<2){
         next
        }
        tempt[k]=chisq.test(frame$NewVit.Flap.from.,frame$Mlap.V4.Flap.to.)[["statistic"]][["X-squared"]]

      }
      RepSVTest=c(RepSVTest,setNames(list(test),paste(names(Model)[i] ) ))
      RepSigTest=c(RepSigTest,setNames(list(tempt),paste(names(Model)[i] ) ))}

  }


}


saveRDS(RepSVTest,file="RepSVTestNTenK.rds")
saveRDS(RepSigTest,file="RepSigTestNTenK.rds")
#saveRDS(SVTest,file="SVTestNTenK.rds")
#saveRDS(SigTest,file="SigTestNTenK.rds")
###SAVESAVE
