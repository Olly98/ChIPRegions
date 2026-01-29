#' ChromGet
#' @description
#' Splits the peak union sequence into sub sequences for the HMM, each sub sequence will be considered an
#' I.I.D sample from the same hidden markov process. This is necessary to remove state transitions
#' without biological meaning from the model, such as between the last peak union on chromosome N and the first of chromosome N+1.
#' By default, sub sequences are only generated for each chromosome, but a maximum gap "MaxGap" between peaks can be
#' specified such that the chain is assumed to "start afresh" if the distance to the next peak is > MaxGap.
#' @export
#' @param UnifiedSequence Output of Unify function, Grange object
#' @param MaxGap Numeric, maximum gap allowed between peak unions within the same subsequence
#' @param MinPeaks Numeric, the minimum number of peaks required to form a valid subsequence
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
