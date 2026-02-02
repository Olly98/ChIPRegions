# ChIPRegions
R package implementing methodology described in the paper 'Systematic Regional Bias is Widespread In ChIP-seq'

## Downloading and Using Package
To install the package, Download the folder entitled 'ChIPRegions' to your device 
and run the following line of R code
```
install.packages("Path_To_ChIPRegions",repos=NULL,type="source")
```
you will also need the following packages from Bioconductor  ```GenomicRanges```, ```IRanges```, ```GenomeInfoDb```, ```S4Vectors```, ```rtracklayer``` and from CRAN you will need ```hmm.discnp``` and ```Signac```
## Using ChIPRegions
Use the following code to import the bedfiles in 'ExampleData':
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
BiasHMM=hmm.discnp::hmm(Chrom[["ChromSplit"]],K=3,itmax=1000)
BiasViterbi=hmm.discnp::viterbi(Chrom[["ChromSplit"]],BiasHMM)
```
Scoring the Regions:
```
ComputedScores=ScoreRegions(UnifiedSequence,Viterbi = BiasViterbi,Bedlist,Scorename = "signalValue")
# The signalValue here is not fold change from control, results will differ from paper
```
Splitting the Genome into Contiguous regions corresponding to each HMM state
```
ContiguousRegions=GetRegions(Chrom[["ChromSplit"]],BiasViterbi,UnifiedSequence)
```
## Running the R Scripts
This repository contains scripts to run the ENCODE wide bias discovery & Bias Quantification for peaks called with MACS2, GEM and SISSRs seen in the paper 'Systematic Regional Bias is Widespread In ChIP-seq'

All raw data and metadata required to run these scripts have been uploaded to Zenodo (DOI: 10.5281/zenodo.18409232). To run the scripts, move them into a folder containing all of the data downloaded from Zenodo, and run on the command line with:
```
Rscript ScriptName.R
```
Please note that some scripts will iteratively download and delete bed files from ENCODE, so avoid running multiple scripts simultaneously in the same folder.

"ENCODEExperimentHMMRunner.R" and "RepSetRunner.R" perform bias discovery between matched ENCODE human TF ChIP-seq experiments, and the discovery of bias between replicates from the same experiment respectively. "ENCODEExperimentRandomHMM.R" and "RepSetsRandom.R" perform "bias discovery" with randomised versions of ENCODE datasets, serving as a negative control (no "bias" should be discovered in random data).

"PeakCallersHMMs.R" performs bias discovery for a subset of 90 ENCODE experiment sets with peaks called by MACS2, GEMR and SISSRs. "PeakCallersHMMsComp.R" performs bias discovery on the subset of peak loci that are unanimously identified across all three peak callers.

