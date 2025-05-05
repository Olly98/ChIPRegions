#' Splits genome into contiguous regions defined by state annotations
#' @export
#' @param ChromSplit A list of sequences of TRUE/FALSE tuples. Obtained after splitting the Peak Union Sequence by chromosome/genomic features.
#' @param Viterbi A list of Viterbi sequences. Each Viterbi sequence corresponds to a sequence in ChromSplit.
#' @param UnifiedSequence a GRanges object containing the unified intervals from which ChromSplit and Viterbi are derived.
GetRegions<-function(ChromSplit,Viterbi,UnifiedSequence){
  UnifiedSequence=GenomeInfoDb::sortSeqlevels(UnifiedSequence)
  UnifiedSequence=sort(UnifiedSequence)
  States=c()
  StateLengths=c()
  ChromLengths=c()
  for(i in 1:length(Viterbi)){
    RLE=S4Vectors::Rle(Viterbi[[i]])
    States=c(States,RLE@values)
    StateLengths=c(StateLengths,RLE@lengths)
    ChromLengths[i]=length(RLE@lengths)
  }

  SegEnds=cumsum(StateLengths)
  SegStarts=SegEnds-(StateLengths)+1
  PeakCoordStart=UnifiedSequence@ranges@start[SegStarts]
  PeakCoordEnd=UnifiedSequence@ranges@start[SegEnds]+UnifiedSequence@ranges@width[SegEnds]
  CleanNames=names(ChromSplit)[is.na(suppressWarnings(as.numeric(unlist(strsplit(names(ChromSplit),split="_")))))]
  Chroms=rep(names(ChromSplit),ChromLengths)

  DataFrame4Granges=data.frame(Chroms,PeakCoordStart,PeakCoordEnd+1)
  names(DataFrame4Granges)=c("chrom","start","stop")
  Regions=GenomicRanges::makeGRangesFromDataFrame(DataFrame4Granges)
  Regions$State=States
  return(Regions)
}


