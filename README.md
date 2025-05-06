# ChIPRegions
R package implementing methodology described in the paper 'Systematic Regional Bias: Ubiquitous in ChIP-seq and Widespread Throughout ENCODE'

## Downloading and Using Package
To install the package, Download the folder entitled 'ChIPRegions' to your device 
and run the following line of R code
```
install.packages("Path_To_ChIPRegions",repos=NULL,type="source")
```
you will also need the following packages from Bioconductor  ```GenomicRanges```, ```IRanges```, ```GenomeInfoDb```, ```S4Vectors```, ```rtracklayer``` and from CRAN you will need ```hmm.discnp``` and ```Signac```
## Using ChIPRegions
Use the following code toimport the bedfiles in 'ExampleData':
```
CEBPGK562_Experiment1_Rep1=import("PathToExampleData/CEBPGK562_Experiment1_Rep1.bed.gz",format="narrowPeak")
CEBPGK562_Experiment1_Rep2=import("PathToExampleData/CEBPGK562_Experiment1_Rep2.bed.gz",format="narrowPeak")
CEBPGK562_Experiment2_Rep1=import("PathToExampleData/CEBPGK562_Experiment2_Rep1.bed.gz",format="narrowPeak")
CEBPGK562_Experiment2_Rep2=import("PathToExampleData/CEBPGK562_Experiment2_Rep2.bed.gz",format="narrowPeak")
```
Then create the Peak Union Sequence & split into Chromosomes:
```
library(ChIPRegions)
Bedlist=list(CEBPGK562_Experiment1_Rep1,CEBPGK562_Experiment1_Rep2,CEBPGK562_Experiment2_Rep1,CEBPGK562_Experiment2_Rep2)
UnifiedSequence=Unify(Bedlist)
Chrom=ChromGet(UnifiedSequence)
```
Then Fit the HMM and Viterbi:
```
BiasHMM=hmm.discnp::hmm(Chrom[["ChromSplit"]],K=3),itmax=1000)
BiasViterbi=hmm.discnp::viterbi(Chrom[["ChromSplit"]],BiasHMM)
```
Scoring the Regions:
```
ComputedScores=ScoreRegions(UnifiedSequence,Viterbi = BiasViterbi,Bedlist,Scorename = "signalValue")
# The signalValue here is not fold change from control, results will differ from paper
```
Splitting the Genome into Contigous regions corresponding to each HMM state
```
ContiguousRegions=GetRegions(Chrom[["ChromSplit"]],BiasViterbi,UnifiedSequence)
```
