library("universalmotif")
library(memes)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(ChIPRegions)
Models=readRDS("/Users/oliverhughes/Downloads/DataFolder/MACOptimalModelExpSets.rds")
RepModels=readRDS("/Users/oliverhughes/Downloads/DataFolder/MACRepSetHMMs.rds")
MACVitRep=readRDS("/Users/oliverhughes/Downloads/DataFolder/RepSetViterbis.rds")
MACRepPUS=readRDS("/Users/oliverhughes/Downloads/DataFolder/MACUnionsRepSet.rds")
MACViterbis=readRDS("/Users/oliverhughes/Downloads/DataFolder/MACViterbis.rds")
Unions=readRDS("/Users/oliverhughes/Downloads/DataFolder/MACUnions.rds")
names(Unions)=names(MACViterbis)
Chroms=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
hs.genome=BSgenome.Hsapiens.UCSC.hg38
AllMotifs=paste("./JASPAR2024_CORE_redundant_pfms_meme/",list.files("./JASPAR2024_CORE_redundant_pfms_meme"),sep="")
AllPWMs=lapply(AllMotifs,read_meme)

GetTFNames=c()
for(i in 1:length(AllPWMs)){
  GetTFNames[i]=AllPWMs[[i]]@altname
}
names(AllPWMs)=GetTFNames


RepMotThresh=list()
for(i in 1:length(MACVitRep)){
  print(i)
  RepMotThresh[[i]]="None"
  names(RepMotThresh)[i]=names(MACVitRep)[i]
  if(is.null(MACVitRep[[i]])){
    next
  }
  Rho=RepModels[[names(MACVitRep)[i]]]$Rho.matrix
  StateNumber=ncol(Rho)
  PUS=MACRepPUS[[names(MACVitRep)[i]]]
  Chrom=ChromGet(PUS,MinPeaks = 3)
  PUS=Chrom[["Unified Sequence"]]
  PUS=PUS[seqnames(PUS)%in%Chroms]
  Chrum=Chroms[Chroms%in%seqnames(PUS)]


  PUS=keepSeqlevels(PUS,Chrum,pruning.mode = "coarse")
  Keeps=which(names(Chrom[["ChromSplit"]])%in%Chroms)
  Vit=MACVitRep[[i]][Keeps]
  TF=unlist(strsplit(names(MACRepPUS)[i],split="-"))[c(T,F)]
  TF=TF[!(TF%in%c("7","S3"))]
  if(!(TF[1]%in%GetTFNames)){
    next
  }

  Vit=unlist(Vit)
  StateMot=rep(0.001,StateNumber)
  Seq=get_sequence(PUS,hs.genome)

  for(j in 1:StateNumber){
     if(sum(Vit==j)==0){
       StateMot[j]=NA
       next
     }
    StateMo=runFimo(Seq[Vit==j],motifs = AllPWMs[names(AllPWMs)==TF[1]],thresh=1e-05)
    if(is.null(StateMo)){
      next
    }
    StateMot[j]=sum(overlapsAny(PUS[Vit==j],StateMo))/sum(Vit==j)

  }
  names(StateMot)=colnames(Rho)
  RepMotThresh[[i]]=StateMot
  names(RepMotThresh)[i]=names(MACVitRep)[i]
}



library("universalmotif")
library(memes)
library(BSgenome.Hsapiens.UCSC.hg38)

Chroms=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
hs.genome=BSgenome.Hsapiens.UCSC.hg38
AllMotifs=paste("/Users/oliverhughes/Downloads/JASPAR2024_CORE_non-redundant_pfms_meme/",list.files("/Users/oliverhughes/Downloads/JASPAR2024_CORE_non-redundant_pfms_meme"),sep="")
AllPWMs=lapply(AllMotifs,read_meme)

GetTFNames=c()
for(i in 1:length(AllPWMs)){
  GetTFNames[i]=AllPWMs[[i]]@altname
}
names(AllPWMs)=GetTFNames


ExpMotThresh=list()
for(i in 1:length(MACViterbis)){
  print(i)
  ExpMotThresh[[i]]="None"
  names(ExpMotThresh)[i]=names(MACViterbis)[i]
  if(is.null(MACViterbis[[i]])){
    next
  }

  Rho=Models[[names(MACViterbis)[i]]]$Rho.matrix
  StateNumber=ncol(Rho)
  Chrom=Unions[[names(MACViterbis)[i]]]
  PUS=Chrom[["Unified Sequence"]]
  PUS=PUS[seqnames(PUS)%in%Chroms]
  Chrum=Chroms[Chroms%in%seqnames(PUS)]


  PUS=keepSeqlevels(PUS,Chrum,pruning.mode = "coarse")
  Keeps=which(names(Chrom[["ChromSplit"]])%in%Chroms)
  Vit=MACViterbis[[i]][Keeps]
  TF=unlist(strsplit(names(MACViterbis)[i],split="-"))[c(T,F)]
  TF=TF[!(TF%in%c("7","S3"))]
  if(!(TF[1]%in%GetTFNames)){
    next
  }

  Vit=unlist(Vit)
  StateMot=rep(0.001,StateNumber)
  Seq=get_sequence(PUS,hs.genome)

  for(j in 1:StateNumber){
    if(sum(Vit==j)==0){
      StateMot[j]=NA
      next
    }
    StateMo=runFimo(Seq[Vit==j],motifs = AllPWMs[names(AllPWMs)==TF[1]],thresh=1e-05)
    if(is.null(StateMo)){
      next
    }
    StateMot[j]=sum(overlapsAny(PUS[Vit==j],StateMo ))/sum(Vit==j)

  }
  names(StateMot)=colnames(Rho)
  ExpMotThresh[[i]]=StateMot
  names(ExpMotThresh)[i]=names(MACViterbis)[i]
}


GetTFNames=c()
for(i in 1:length(AllPWMs)){
  GetTFNames[i]=AllPWMs[[i]]@altname
}
names(AllPWMs)=GetTFNames


RepMotGC=list()
for(i in 1:length(MACVitRep)){
  print(i)
  RepMotGC[[i]]="None"
  names(RepMotGC)[i]=names(MACVitRep)[i]
  if(is.null(MACVitRep[[i]])){
    next
  }
  Rho=RepModels[[names(MACVitRep)[i]]]$Rho.matrix
  StateNumber=ncol(Rho)
  PUS=MACRepPUS[[names(MACVitRep)[i]]]
  Chrom=ChromGet(PUS,MinPeaks = 3)
  PUS=Chrom[["Unified Sequence"]]
  PUS=PUS[seqnames(PUS)%in%Chroms]
  Chrum=Chroms[Chroms%in%seqnames(PUS)]


  PUS=keepSeqlevels(PUS,Chrum,pruning.mode = "coarse")
  Keeps=which(names(Chrom[["ChromSplit"]])%in%Chroms)
  Vit=MACVitRep[[i]][Keeps]


  Vit=unlist(Vit)
  StateMot=rep(NA,StateNumber)
  Seq=get_sequence(PUS,hs.genome)

  for(j in 1:StateNumber){



    StateMot[j]=sum(letterFrequency(Seq[Vit==j],letters=c("G","C")))/sum(lengths(Seq[Vit==j]))


  }
  names(StateMot)=colnames(Rho)
  RepMotGC[[i]]=StateMot
  names(RepMotGC)[i]=names(MACVitRep)[i]
}



GetTFNames=c()
for(i in 1:length(AllPWMs)){
  GetTFNames[i]=AllPWMs[[i]]@altname
}
names(AllPWMs)=GetTFNames


ExpMotGC=list()
for(i in 1:length(MACViterbis)){
  print(i)
  ExpMotGC[[i]]="None"
  names(ExpMotGC)[i]=names(MACViterbis)[i]
  if(is.null(MACViterbis[[i]])){
    next
  }
  Rho=Models[[names(MACViterbis)[i]]]$Rho.matrix
  StateNumber=ncol(Rho)

  Chrom=Unions[[names(MACViterbis)[i]]]
  PUS=Chrom[["Unified Sequence"]]
  PUS=PUS[seqnames(PUS)%in%Chroms]
  Chrum=Chroms[Chroms%in%seqnames(PUS)]


  PUS=keepSeqlevels(PUS,Chrum,pruning.mode = "coarse")
  Keeps=which(names(Chrom[["ChromSplit"]])%in%Chroms)
  Vit=MACViterbis[[i]][Keeps]


  Vit=unlist(Vit)
  StateMot=rep(NA,StateNumber)

  Seq=get_sequence(PUS,hs.genome)

  for(j in 1:StateNumber){




    StateMot[j]=sum(letterFrequency(Seq[Vit==j],letters=c("G","C")))/sum(lengths(Seq[Vit==j]))


  }
  names(StateMot)=colnames(Rho)
  ExpMotGC[[i]]=StateMot
  names(ExpMotGC)[i]=names(MACViterbis)[i]
}

saveRDS(RepMotThresh,file="RepMotThresh.rds")
saveRDS(ExpMotThresh,file="ExpMotThresh.rds")
saveRDS(RepMotGC,file="RepMotGC.rds")
saveRDS(ExpMotGC,file="ExpMotGC.rds")


