
SetsMeta=read.csv("./CCSets.csv")
AllMeta=read.csv("./ENCODEExperimentsMeta.csv")
ControlDownloads=read.csv("./ControlsAndDownloadLinks.csv")
TreatsAndURLs=read.csv("./TreatsAndURLs.csv")
AllList=list()
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(tidyverse)

library(rtracklayer)
library(GenomicRanges)
library(ChIPRegions)
library(hmm.discnp)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
ErrorCatcherReps=c()
Sets=read.csv("./CCSets.csv")

Model=readRDS("./CCViterbis.rds")
PUS=readRDS("./NewPUSCCV3.rds")
PUS=PUS[-12]
Model <- Model[!sapply(Model,is.null)]

Viterbis=list()
Regions=list()
Issue=c()
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

Unify=function (Bedlist)
{
  PUS = Signac::UnifyPeaks(Bedlist, mode = "reduce")
  PUS$OverlapPattern = do.call(paste, lapply(Bedlist, IRanges::overlapsAny,
    query = PUS))
  return(PUS)
}



for(i in 1:length(Model)){
if(is.null(Model[[i]])){
  next
}

print(rep(i,15))
Viter=unlist(Model[[i]])
Chrom=ChromGet(PUS[[i]])
Cell=Sets$Cell[i]

RNAFiles=list.files(paste("./",Cell,"RNA","/",sep=""))
RNAFiles=RNAFiles[grep("tsv",RNAFiles)]
States=sort(unique(Viter))

StateList=list()
for(j in (States)){
FileList=list()
Anno=annotatePeak(Chrom[["Unified Sequence"]][Viter==j], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

for(k in 1:length(RNAFiles)){
gg=read_tsv(paste("./",Cell,"RNA","/",RNAFiles[k],sep=""))
gg=gg[grep("EN",gg$gene_id),]
gg$gene_id=substr(gg$gene_id,1,15)
ExpL=gg$pme_TPM[gg$gene_id%in%na.omit(Anno@anno@elementMetadata@listData[["ENSEMBL"]])]
names(ExpL)=gg$gene_id[gg$gene_id%in%na.omit(Anno@anno@elementMetadata@listData[["ENSEMBL"]])]
FileList[[k]]=ExpL

}
StateList[[j]]=FileList


}

AllList[[i]]=StateList

}

saveRDS(AllList,file="RNACC.rds")

###SAVESAVE

