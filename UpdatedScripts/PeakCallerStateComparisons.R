library(GenomicRanges)
library(ChIPRegions)


AllGEMPUS=readRDS("./GEMSharePUS.rds")
AllMACPUS=readRDS("./MACSharePUS.rds")
AllSISSPUS=readRDS("./SISSSharePUS.rds")
AllGEMModels=readRDS("./GEMShareOptimalModel.rds")
AllMACModels=readRDS("./MACShareOptimalModel.rds")
AllSISSModels=readRDS("./SISSShareOptimalModel.rds")

Keeps=c()

for(i in 1:90){
  GEM=length(AllGEMModels[[i]])
  MAC=length(AllMACModels[[i]])
  SISS=length(AllSISSModels[[i]])
  if(sum(c(MAC,GEM,SISS)<19)<1){
    Keeps[i]=TRUE
  }else{
    Keeps[i]=FALSE
  }
}

GEMVitPUS=list()
MACVitPUS=list()
SISSVitPUS=list()


for(i in 1:90){
  if(Keeps[i]==FALSE){
    next
  }

  GEMPUS=AllGEMPUS[[i]]
  MACPUS=AllMACPUS[[i]]
  SISSPUS=AllSISSPUS[[i]]


  GEMMAC=GEMPUS[overlapsAny(GEMPUS,MACPUS)]

  All=GEMMAC[overlapsAny(GEMMAC,SISSPUS)]


  GEMChrom=ChromGet(All)
  MACChrom=ChromGet(MACPUS[overlapsAny(MACPUS,All)])
  SISSChrom=ChromGet(SISSPUS[overlapsAny(SISSPUS,All)])

  SISSVit=unlist(viterbi(SISSChrom[["ChromSplit"]],AllSISSModels[[i]]))
  SISSOne=SISSChrom[["Unified Sequence"]]
  SISSOne$State=SISSVit
  SISSVitPUS[[i]]=SISSOne

  GEMVit=unlist(viterbi(GEMChrom[["ChromSplit"]],AllGEMModels[[i]]))
  GEMOne=GEMChrom[["Unified Sequence"]]
  GEMOne$State=GEMVit
  GEMVitPUS[[i]]=GEMOne

  MACVit=unlist(viterbi(MACChrom[["ChromSplit"]],AllMACModels[[i]]))
  MACOne=MACChrom[["Unified Sequence"]]
  MACOne$State=MACVit
  MACVitPUS[[i]]=MACOne

}

SISSMACComp=list()
for(i in 1:90){
  if(Keeps[i]==FALSE){
    next
  }

  GERM=SISSVitPUS[[i]]

  MARC=MACVitPUS[[i]]
  MACState=unique(MARC$State)
  if(length(MACState)<2){
    next
  }
  GEMState=unique(GERM$State)
  if(length(GEMState)<2){
    next
  }
  if(length(MACState)<length(GEMState)){
    StateMatProp=matrix(nrow=length(unique(MARC$State)),ncol=length(unique(GERM$State)))
    StateMatCount=matrix(nrow=length(unique(MARC$State)),ncol=length(unique(GERM$State)))
    for(j in MACState){
      MACSub=MARC[MARC$State==j]
      for(k in GEMState){
        GEMSub=GERM[GERM$State==k]
        StateMatProp[as.numeric(j),as.numeric(k)]=sum(overlapsAny(MACSub,GEMSub))/length(MACSub)
        StateMatCount[as.numeric(j),as.numeric(k)]=sum(overlapsAny(MACSub,GEMSub))
      }

    }
  }else{
    StateMatProp=matrix(nrow=length(unique(GERM$State)),ncol=length(unique(MARC$State)))
    StateMatCount=matrix(nrow=length(unique(GERM$State)),ncol=length(unique(MARC$State)))
    for(j in 1:length(GEMState)){
      GEMSub=GERM[GERM$State==GEMState[j]]
      for(k in 1:length(MACState)){
        MACSub=MARC[MARC$State==MACState[k]]
        StateMatProp[j,k]=sum(overlapsAny(GEMSub,MACSub))/length(GEMSub)
        StateMatCount[j,k]=sum(overlapsAny(GEMSub,MACSub))
      }
    }}




  Combined=list(StateMatCount,StateMatProp)
  SISSMACComp[[i]]=Combined
}

AllStateListSISSMAC=list()
OverallListSISSMAC=list()

for(i in 1:length(SISSMACComp)){
  if(Keeps[i]==FALSE){
    next
  }
  GERM=SISSVitPUS[[i]]
  MARC=MACVitPUS[[i]]
  MACState=unique(MARC$State)
  if(length(MACState)<2){
    next
  }
  GEMState=unique(GERM$State)
  if(length(GEMState)<2){
    next
  }
  StateCorrelation=c()
  Combie=SISSMACComp[[i]][[2]]
  GEMEmits=AllSISSModels[[i]]$Rho.matrix
  GEMSymbs=rownames(GEMEmits)
  MACEmits=AllMACModels[[i]]$Rho.matrix
  MACSymbs=rownames(MACEmits)
  AllSymbs=MACSymbs[MACSymbs%in%GEMSymbs]
  GEMEmits=GEMEmits[(AllSymbs),]
  MACEmits=MACEmits[(AllSymbs),]
  if(length(MACState)<=length(GEMState)){
    for(j in 1:nrow(Combie)){
      Comp=which(Combie[j,]==max(Combie[j,]))

      if(ncol(GEMEmits)>2){
        GEMMeanOthers=rowMeans(GEMEmits[,-Comp])
      }else{
        GEMMeanOthers=GEMEmits[,-Comp]
      }
      if(ncol(MACEmits)>2){
        MACMeanOthers=rowMeans(MACEmits[,-j])
      }else{
        MACMeanOthers=MACEmits[,-j]
      }
      MACContrast=MACEmits[,j]
      GEMContrast=GEMEmits[,Comp]
      StateCorrelation[j]=cor(log(MACContrast),log(GEMContrast),method="spearman")
    }}else{
      for(j in 1:nrow(Combie)){
        Comp=which(Combie[j,]==max(Combie[j,]))

        if(ncol(MACEmits)>2){
          MACMeanOthers=rowMeans(MACEmits[,-Comp])
        }else{
          MACMeanOthers=MACEmits[,-Comp]
        }
        if(ncol(GEMEmits)>2){
          GEMMeanOthers=rowMeans(GEMEmits[,-j])
        }else{
          GEMMeanOthers=GEMEmits[,-j]
        }
        MACContrast=MACEmits[,Comp]
        GEMContrast=GEMEmits[,j]
        StateCorrelation[j]=cor(log(MACContrast),log(GEMContrast),method="spearman")
      }

    }
  AllStateListSISSMAC[[i]]=StateCorrelation
  OverallListSISSMAC[[i]]=median(StateCorrelation)
}

GEMMACComp=list()
for(i in 1:90){
  if(Keeps[i]==FALSE){
    next
  }

  GERM=GEMVitPUS[[i]]

  MARC=MACVitPUS[[i]]
  MACState=unique(MARC$State)
  if(length(MACState)<2){
    next
  }
  GEMState=unique(GERM$State)
  if(length(GEMState)<2){
    next
  }
  if(length(MACState)<length(GEMState)){
    StateMatProp=matrix(nrow=length(unique(MARC$State)),ncol=length(unique(GERM$State)))
    StateMatCount=matrix(nrow=length(unique(MARC$State)),ncol=length(unique(GERM$State)))
    for(j in MACState){
      MACSub=MARC[MARC$State==j]
      for(k in GEMState){
        GEMSub=GERM[GERM$State==k]
        StateMatProp[as.numeric(j),as.numeric(k)]=sum(overlapsAny(MACSub,GEMSub))/length(MACSub)
        StateMatCount[as.numeric(j),as.numeric(k)]=sum(overlapsAny(MACSub,GEMSub))
      }

    }
  }else{
    StateMatProp=matrix(nrow=length(unique(GERM$State)),ncol=length(unique(MARC$State)))
    StateMatCount=matrix(nrow=length(unique(GERM$State)),ncol=length(unique(MARC$State)))
    for(j in 1:length(GEMState)){
      GEMSub=GERM[GERM$State==GEMState[j]]
      for(k in 1:length(MACState)){
        MACSub=MARC[MARC$State==MACState[k]]
        StateMatProp[j,k]=sum(overlapsAny(GEMSub,MACSub))/length(GEMSub)
        StateMatCount[j,k]=sum(overlapsAny(GEMSub,MACSub))
      }
    }}




  Combined=list(StateMatCount,StateMatProp)
  GEMMACComp[[i]]=Combined
}

AllStateListGEMMAC=list()
OverallListGEMMAC=list()

for(i in 1:length(GEMMACComp)){
  if(Keeps[i]==FALSE){
    next
  }
  GERM=GEMVitPUS[[i]]
  MARC=MACVitPUS[[i]]
  MACState=unique(MARC$State)
  if(length(MACState)<2){
    next
  }
  GEMState=unique(GERM$State)
  if(length(GEMState)<2){
    next
  }
  StateCorrelation=c()
  Combie=GEMMACComp[[i]][[2]]
  GEMEmits=AllGEMModels[[i]]$Rho.matrix
  GEMSymbs=rownames(GEMEmits)
  MACEmits=AllMACModels[[i]]$Rho.matrix
  MACSymbs=rownames(MACEmits)
  AllSymbs=MACSymbs[MACSymbs%in%GEMSymbs]
  GEMEmits=GEMEmits[(AllSymbs),]
  MACEmits=MACEmits[(AllSymbs),]
  if(length(MACState)<=length(GEMState)){
    for(j in 1:nrow(Combie)){
      Comp=which(Combie[j,]==max(Combie[j,]))

      if(ncol(GEMEmits)>2){
        GEMMeanOthers=rowMeans(GEMEmits[,-Comp])
      }else{
        GEMMeanOthers=GEMEmits[,-Comp]
      }
      if(ncol(MACEmits)>2){
        MACMeanOthers=rowMeans(MACEmits[,-j])
      }else{
        MACMeanOthers=MACEmits[,-j]
      }
      MACContrast=MACEmits[,j]
      GEMContrast=GEMEmits[,Comp]
      StateCorrelation[j]=cor(log(MACContrast),log(GEMContrast),method="spearman")
    }}else{
      for(j in 1:nrow(Combie)){
        Comp=which(Combie[j,]==max(Combie[j,]))

        if(ncol(MACEmits)>2){
          MACMeanOthers=rowMeans(MACEmits[,-Comp])
        }else{
          MACMeanOthers=MACEmits[,-Comp]
        }
        if(ncol(GEMEmits)>2){
          GEMMeanOthers=rowMeans(GEMEmits[,-j])
        }else{
          GEMMeanOthers=GEMEmits[,-j]
        }
        MACContrast=MACEmits[,Comp]
        GEMContrast=GEMEmits[,j]
        StateCorrelation[j]=cor(log(MACContrast),log(GEMContrast),method="spearman")
      }

    }
  AllStateListGEMMAC[[i]]=StateCorrelation
  OverallListGEMMAC[[i]]=median(StateCorrelation)
}


GEMSISSComp=list()
for(i in 1:90){
  if(Keeps[i]==FALSE){
    next
  }

  GERM=GEMVitPUS[[i]]

  MARC=SISSVitPUS[[i]]
  SISSState=unique(MARC$State)
  if(length(SISSState)<2){
    next
  }
  GEMState=unique(GERM$State)
  if(length(GEMState)<2){
    next
  }
  if(length(SISSState)<=length(GEMState)){
    StateMatProp=matrix(nrow=length(unique(MARC$State)),ncol=length(unique(GERM$State)))
    StateMatCount=matrix(nrow=length(unique(MARC$State)),ncol=length(unique(GERM$State)))
    for(j in SISSState){
      SISSSub=MARC[MARC$State==j]
      for(k in GEMState){
        GEMSub=GERM[GERM$State==k]
        StateMatProp[as.numeric(j),as.numeric(k)]=sum(overlapsAny(SISSSub,GEMSub))/length(SISSSub)
        StateMatCount[as.numeric(j),as.numeric(k)]=sum(overlapsAny(SISSSub,GEMSub))
      }

    }
  }else{
    StateMatProp=matrix(nrow=length(unique(GERM$State)),ncol=length(unique(MARC$State)))
    StateMatCount=matrix(nrow=length(unique(GERM$State)),ncol=length(unique(MARC$State)))
    for(j in 1:length(GEMState)){
      GEMSub=GERM[GERM$State==GEMState[j]]
      for(k in 1:length(SISSState)){
        SISSSub=MARC[MARC$State==SISSState[k]]
        StateMatProp[j,k]=sum(overlapsAny(GEMSub,SISSSub))/length(GEMSub)
        StateMatCount[j,k]=sum(overlapsAny(GEMSub,SISSSub))
      }
    }}




  Combined=list(StateMatCount,StateMatProp)
  GEMSISSComp[[i]]=Combined
}



