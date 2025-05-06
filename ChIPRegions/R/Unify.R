#' Creates the Peak Union Sequence
#' @export
#' @param Bedlist A list of Grange objects
Unify<-function(Bedlist){
  PUS=Signac::UnifyPeaks(Bedlist,mode="reduce")
  PUS$OverlapPattern=do.call(paste,lapply(Bedlist,IRanges::overlapsAny,query=PUS))
  return(PUS)
}
