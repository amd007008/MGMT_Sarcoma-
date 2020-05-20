# # # # # # # # # # # #
# Pipeline MGMT_STP27 #
# # # # # # # # # # # # 
install.packages(c("ade4","MASS"))
BiocManager::install("lumi")
BiocManager::install("methylumi")
BiocManager::install("minfi")

library(ade4)
library(MASS)
library(methylumi)
library(minfi)
library(lumi)
library(devtools)
library(remotes)
library(BiocManager)

install.packages("~/Downloads/mgmtstp27_0.6-3.tar.gz",repos=NULL)
install.packages("~/Downloads/mgmtstp27_0.6-3.zip", repos=NULL)

library(devtools)
#Data Import: 

##### Example: ##### 
# loading R packages
require(mgmtstp27)
require(minfiData)
# preprocessing of the data
dat <- preprocessRaw(RGsetEx)
# computation of M-value
mvalue <- log2((getMeth(dat)+1)/(getUnmeth(dat)+1))
mvalue <- as.data.frame(t(mvalue))
# predictions
pred1 <- MGMTpredict(mvalue)
head(pred1)
# quality control graphics
par(mfrow=c(2,3))
MGMTqc.pop(pred1,which.plot=1:3,mfrow=NULL)
MGMTqc.single(pred1,nsample=1,which.plot=1:3,mfrow=NULL)

##### Procedure for Preprocessing and MGMT Promoter Methylation Prediction: 
#conversion of red/green channel information into signals for methylated and unmethylated: 

## --> MGMTpredict-funktion 
#M-Values: (log2-ratio of methylated and unmethylated intensities corrected by an offset equal to 1) 
#located in the MGMT promoter, cg12434587 and cg12981137: 

#MGMT score: logit transformation of probability that theMGMT promoter is methylated to obtain a quasi-normal score 

#####Effect of normalization #####

#####Gene CNAs (Feber et al.)##### 
#calculation: 

#adapted to 850k? 

#normalization for each sample individually: 

#combined intensities for methylated and unmethylated (total intensity, T) calculated from normalized intensities 

#value log2(R) in Bady et al.: defined as difference of intensity between samples and a synthetic reference 
 #log2(R) = log2(T(observed) + 1) - log2(T(reference) + 1)
#value log2(R) here: ? 

#soothing procedure? 

#corrected with scaling factor method to reduce chemistry type-bias of total intensity: 

#probes with nonsignificant P values(typically >0.01) excluded from analysis when raw data served as input: 

##### Determination of CNA State #####
#R package: CGHcall?  --> uses circular binary segmentation (starting with normalized log 2(R) values for each sample)
#for genomic region: CNA events detected using copy number probe means (CpGs) contained in selected region 

#####Statistics #####
#CIMP+tumours were identified with unsupervised clustering methods: (Ward's algorithm with Euclidean distance) 
#Relations between categorical variables were assessed by X^2 tests with P values computed by Monte Carlo simulation - (cell counts less than five)
#two way ANOVA-like approach instead of two-way analysis of variance for testing effects of CNA and DNA methylation on expression of MGMT from F statistics (robuster for unbalanced data)

#R-Packages: Optimal Cutpoints and epiR for evalutaion of cutoff robustness, including determination of optimal values and performances 





