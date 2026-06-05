
library(memes)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(ChIPRegions)
Models=readRDS("./MACOptimalModelExpSets.rds")
RepModels=readRDS("./MACRepSetHMMs.rds")

Starts=c(1,3,5,7,9)
Stops=c(2,4,6,8,10)
ExpContrasts=list()
for(i in 1:length(Models)){
  if(length(Models[[i]])<19){
    next
  }
  Matr=Models[[i]]$Rho.matrix
  Trable=(read.table(text=rownames(Matr)))
  ExpNum=length(unlist(strsplit(rownames(Matr)[1],split=" ") ))/2
  if(!(ExpNum%in%c(2,3,4,5))){

    next
  }
  SperUn=c()
  for(j in 1:ncol(Matr)){
    Cands=c()
    for(k in 1:ExpNum){
      Sub=Matr[rowSums(Trable[,c(1:ncol(Trable))[-c(Starts[k],Stops[k])]])==0,j]
      Cands[k]=sum(Sub)
    }
    SperUn[j]=max(Cands)
  }
  ExpContrasts[[i]]=SperUn

}

Starts=c("TRUE FALSE","FALSE TRUE","TRUE FALSE FALSE","FALSE TRUE FALSE","FALSE FALSE TRUE")

RepContrasts=list()
for(i in 1:length(RepModels)){
  if(length(RepModels[[i]])<19){
    RepContrasts[[i]]="None"
    next
  }
  Matr=RepModels[[i]]$Rho.matrix
  Trable=(read.table(text=rownames(Matr)))
  #ExpNum=length(unlist(strsplit(rownames(Matr)[1],split=" ") ))
  if((length(rownames(Matr) )>7)){
    print("k")
    next

  }

  SperUn=c()
  Barts=Starts[Starts%in%rownames(Matr)]
  for(j in 1:ncol(Matr)){
    Cands=c()
    for(k in 1:length(Barts)){
      Sub=Matr[Barts[k],j]
      Cands[k]=sum(Sub)
    }
    SperUn[j]=max(Cands)
  }
  RepContrasts[[i]]=SperUn

}

names(RepContrasts)=names(RepModels)
names(ExpContrasts)=names(Models)
saveRDS(ExpContrasts,file="ExpBiasSeverity.rds")
saveRDS(RepContrasts,file="RepBiasSeverity.rds")

