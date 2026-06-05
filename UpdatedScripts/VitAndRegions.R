
SetsMeta=read.csv("./ExpSets.csv")
AllMeta=read.csv("./ENCODEExperimentsMeta.csv")
HMM=readRDS("./OptimalHMMsExpSets.rds")
library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
OptimalModel=list()
OptimalState=c()


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


Viterbis=list()
Regions=list()
for(i in 1:length(HMM)){
if(length(HMM[[i]])<19){
next
}
SubMeta=AllMeta[AllMeta$Experiment.target==SetsMeta$UnfoldedTFs[i],]

SubMeta=SubMeta[SubMeta$Biosample.term.name==SetsMeta$UnfoldedCells[i],]
print(length(SubMeta$File.download.URL))
FileList=list()
FileLengths=c()
for(j in 1:length(SubMeta$File.download.URL)){
BamDownloadCall=paste("wget",SubMeta$File.download.URL[j],sep=" ")
print(BamDownloadCall)
system(BamDownloadCall)
FileList[[j]]=import(paste(SubMeta$File.accession[j],".bed.gz",sep=""),format="narrowPeak")
FileLengths[j]=length(FileList[[j]])
}

PUS=Unify(FileList)
print(length(PUS))
Chrom=ChromGet(PUS,MinPeaks=3)

Viterbis[[i]]=viterbi(Chrom[["ChromSplit"]],HMM[[i]])

Regions[[i]]=GetRegions(Chrom[["ChromSplit"]],Viterbis[[i]],Chrom[["Unified Sequence"]])
system("rm *bed.gz")
}
saveRDS(Viterbis,file="ExpSetViterbis.rds")
saveRDS(Regions,file="ExpSetRegions.rds")
