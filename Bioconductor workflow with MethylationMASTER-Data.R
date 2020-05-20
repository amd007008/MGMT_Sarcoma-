#Import Data: 
library(tidyverse)
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(ade4)
library(MASS)
library(methylumi)
library(minfi)
library(lumi)
library(devtools)
library(remotes)
library(BiocManager)
library(limma)
library(wateRmelon)
library(readxl)
library(RPMM)
library(FlowSorted.Blood.EPIC)
#library(FlowSorted.Blood.450k)
#library(IlluminaHumanMethylation450kmanifest)
#library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

out.dir <- "Users/Anni/MD_thesis/DA/R/MGMT-Sarcoma-/"

mytargets <- read.metharray.sheet(file.path(getwd(), "MASTER_Methylation_raw"))

myrgSet <- read.metharray.exp(targets = mytargets, force = TRUE)

#class: RGChannelSet 
#dim: 1051539 19 
#metadata(0):
#  assays(2): Green Red
#rownames(1051539): 1600101 1600111 ... 99810990 99810992
#rowData names(0):
#  colnames(19): 201465950024_R01C01 201465950024_R07C01 ...
#203257020125_R01C01 203982200140_R02C01
#colData names(9): Sample_Name Sample_Well ... Basename filenames
#Annotation
#array: IlluminaHumanMethylationEPIC
#annotation: ilm10b4.hg19

mytargets$ID <- paste(mytargets$Sample_Name,mytargets$Array,sep=".")
sampleNames(myrgSet) <- mytargets$ID
myrgSet

#class: RGChannelSet 
#dim: 1051539 19 
#metadata(0):
 # assays(2): Green Red
#rownames(1051539): 1600101 1600111 ... 99810990 99810992
#rowData names(0):
#  colnames(19): H021_06H1_M1_D1.R01C01 H021_06H1_M1_D1.R07C01 ...
#H021_XQGYFN_T2_D1.R01C01 H021_XQGYFN_T5_D1.R02C01
#colData names(9): Sample_Name Sample_Well ... Basename filenames
#Annotation
#array: IlluminaHumanMethylationEPIC
#annotation: ilm10b4.hg19

# calculate the detection p-values
mydetP <- detectionP(myrgSet)
head(mydetP)


# examine mean detection p-values across all samples to identify any failed samples
# Scaled Plot [0, 0.03]
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)


plot2 <- ggplot(as.data.frame(colMeans(mydetP)),
                aes(
                  x = names(colMeans(mydetP)),
                  y = as.numeric(colMeans(mydetP)),
                  fill = as.factor(mytargets$Sample_Name)
                )
)+
  xlab(NULL)+
  ylab("Mean detection p-values")+
  scale_y_continuous(limits=c(0, 0.003))+
  geom_col()+
  scale_fill_manual(values=mycolors)+
  theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1))


pdf(
  file=paste0(out.dir, "Barplot_MethylationMASTER2", Sys.Date(), ".pdf"),
  paper = "a4"
) 
plot2
dev.off()

barplot(colMeans(mydetP), col=pal[factor(mytargets$Sample_Name)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")

#Further quality control plots by minfi::qcReport
qcReport(myrgSet, sampNames=mytargets$ID, sampGroups=mytargets$Array, 
         pdf="qcReport.pdf")

# Plotting all detection p-values for each sample to gauge measured effects
#Unscaled Plot
Unscaled_pValue_Sample <- ggplot(as.data.frame(mydetP))+
  geom_histogram(aes(x = H021_06H1_M1_D1.R01C01),binwidth = 0.01, fill = mypal[1])+
  geom_histogram(aes(x = H021_06H1_M1_D1.R07C01),binwidth = 0.01, fill = mypal[2])+
  geom_histogram(aes(x = H021_1KJALE_M1_D1.R01C01),binwidth = 0.01, fill = mypal[3])+
  geom_histogram(aes(x = H021_50B6_T1_D1.R02C01),binwidth = 0.01, fill = mypal[4])+
  geom_histogram(aes(x = H021_6HPT4M_M1_D1.R02C01),binwidth = 0.01, fill = mypal[5])+
  geom_histogram(aes(x = H021_7UDL_T4_D1.R07C01),binwidth = 0.01, fill = mypal[6])+
  geom_histogram(aes(x = H021_7UDL_T4_D1.R01C01),binwidth = 0.01, fill = mypal[7])+
  geom_histogram(aes(x = H021_H7T59N_T1_D1.R04C01),binwidth = 0.01, fill = mypal[8])+
  geom_histogram(aes(x = H021_K4WF_T2_D1.R02C01),binwidth = 0.01, fill = mypal[9])+
  geom_histogram(aes(x = H021_KNSV_T1_D1.R03C01),binwidth = 0.01, fill = mypal[10])+
  geom_histogram(aes(x = H021_KNSV_T1_D2.R03C01),binwidth = 0.01, fill = mypal[11])+
  geom_histogram(aes(x = H021_RCJFWM_T1_D1.R08C01),binwidth = 0.01, fill = mypal[12])+
  geom_histogram(aes(x = H021_RCJFWM_T1_D1.R01C01),binwidth = 0.01, fill = mypal[13])+
  geom_histogram(aes(x = H021_RCJFWM_T2_D1.R01C01),binwidth = 0.01, fill = mypal[14])+
  geom_histogram(aes(x = H021_RF6GGV_M1_D1.R05C01),binwidth = 0.01, fill = mypal[15])+
  geom_histogram(aes(x = H021_RVV7AR_M1_D1.R03C01),binwidth = 0.01, fill = mypal[16])+
  geom_histogram(aes(x = H021_XQGYFN_T2_D1.R07C01),binwidth = 0.01, fill = mypal[17])+
  geom_histogram(aes(x = H021_XQGYFN_T2_D1.R01C01),binwidth = 0.01, fill = mypal[18])+
  geom_histogram(aes(x = H021_XQGYFN_T5_D1.R02C01),binwidth = 0.01, fill = mypal[19])+
  xlab("Detection p-values")


# remove poor quality samples
keep <- colMeans(mydetP) < 0.01
myrgSet <- myrgSet[,keep]

# remove poor quality samples from targets data
mytargets <- mytargets[keep,]
rm(keep)

# Normalization -----------------------------------------------------------

#good rule of thumb within the minfi framework:
#preprocessFunnorm most appropriate for datasets with global methylation differences
#         such as cancer/normal or vastly different tissue types (Fortin et al. 2014)
#preprocessQuantile more suited for datasets where you do not expect global differences
#           between your samples, for example a single tissue (Touleimat and Tost 2012)
#Use Quantro to test whether to use global normalization methods, such as quantile normalization
#           (Quantro essentially runs an ANOVA and compares sample groups)

# normalize the data; this results in a GenomicRatioSet object
mymSetSq <- preprocessQuantile(myrgSet)

# create a MethylSet object from the raw data for plotting
mymSetRaw <- preprocessRaw(myrgSet)

# visualise what the data looks like before and after normalization
setwd("Output")
pdf("data_normalization.pdf", paper = "a4r")
par(mfrow=c(1,2))
densityPlot(myrgSet,
            sampGroups = mytargets$Array,
            main = "Raw",
            legend = FALSE)
legend("topleft",
       legend = levels(factor(mytargets$Array)),
       text.col = mycolors)
densityPlot(getBeta(mymSetSq),
            sampGroups = mytargets$Array,
            main = "Normalized", legend = FALSE)
legend("topleft",
       legend = levels(factor(mytargets$Array)),
       text.col = mycolors)
dev.off()
setwd(old.dir)

#####this is how far I've come the first day######

# Data Exploration --------------------------------------------------------

# MDS plots to look at largest sources of variation
setwd("Output")
pdf("MDS.pdf", paper = "a4")
par(mfrow=c(2,1))
plotMDS(getM(mymSetSq), pch = rep(0:1,c(4,4)), top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)])
legend("bottomright", legend=levels(factor(mytargets$Sample_Group)), pch = 0:1,
       text.col=pal, bg="transparent", cex=0.7)

plotMDS(getM(mymSetSq), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Name)])
legend("bottomright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18, text.col=pal,
       bg="transparent", cex=0.7)
dev.off()
setwd(old.dir)

# Examine higher dimensions to look at other sources of variation
setwd("Output")
pdf("MDS_higher-dim.pdf", paper = "a4r")
par(mfrow=c(1,3))
plotMDS(getM(mymSetSq), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)], dim=c(1,3))
legend("topleft", legend=levels(factor(mytargets$Sample_Group)), text.col=pal,
       cex=0.7, bg="transparent")
legend("topright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18,
       cex=0.7, bg="transparent")

plotMDS(getM(mymSetSq), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)], dim=c(2,3))
legend("top", legend=levels(factor(mytargets$Sample_Group)), text.col=pal,
       cex=0.7, bg="transparent")
legend("topright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18,
       cex=0.7, bg="transparent")

plotMDS(getM(mymSetSq), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)], dim=c(3,4))
legend("topleft", legend=levels(factor(mytargets$Sample_Group)), text.col=pal,
       cex=0.7, bg="transparent")
legend("topright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18,
       cex=0.7, bg="transparent")
dev.off()
setwd(old.dir)


# Filtering ---------------------------------------------------------------

# ensure probes are in the same order in the mymSetSq and mydetP objects
mydetP <- mydetP[match(featureNames(mymSetSq),rownames(mydetP)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(mydetP < 0.01) == ncol(mymSetSq)
mymSetSqFlt <- mymSetSq[keep,]
rm(keep)

# # if your data includes males and females, remove probes on the sex chromosomes
# keep <- !(featureNames(mymSetSqFlt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
# mymSetSqFlt <- mymSetSqFlt[keep,]
# rm(keep)

# remove probes with SNPs at CpG site (what exactly is the default source???)
mymSetSqFlt <- dropLociWithSnps(mymSetSqFlt)

# exclude cross reactive probes (Cartney et al. 2016, Pidsley et al. 2016)
# Better Sources?
xReactiveProbes <- bind_rows(
  read.delim(file = "Cartney_xReactiveProbes.txt", header = FALSE,
             col.names = "featureNames"),
  read.csv(file = "Pidsley_xReactiveProbes.csv", header = FALSE,
           col.names = c("featureNames", rep("V", 5)),colClasses=c(NA, rep("NULL", 5))))
keep <- !(featureNames(mymSetSqFlt) %in% xReactiveProbes$featureNames)
mymSetSqFlt <- mymSetSqFlt[keep,]
rm(keep, xReactiveProbes)

# re-examine the MDS plots to see if the relationship between the samples has changed

# MDS plots to look at largest sources of variation
setwd("Output")
pdf("filtered_MDS.pdf", paper = "a4")
par(mfrow=c(2,1))
plotMDS(getM(mymSetSqFlt), pch = rep(0:1,c(4,4)), top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)])
legend("bottomright", legend=levels(factor(mytargets$Sample_Group)), pch = 0:1,
       text.col=pal, bg="transparent", cex=0.7)

plotMDS(getM(mymSetSqFlt), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Name)])
legend("bottomright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18, text.col=pal,
       bg="transparent", cex=0.7)
dev.off()
setwd(old.dir)

# Examine higher dimensions to look at other sources of variation
setwd("Output")
pdf("filtered_MDS_higher-dim.pdf", paper = "a4r")
par(mfrow=c(1,3))
plotMDS(getM(mymSetSqFlt), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)], dim=c(1,3))
legend("topleft", legend=levels(factor(mytargets$Sample_Group)), text.col=pal,
       cex=0.7, bg="transparent")
legend("topright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18,
       cex=0.7, bg="transparent")

plotMDS(getM(mymSetSqFlt), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)], dim=c(2,3))
legend("top", legend=levels(factor(mytargets$Sample_Group)), text.col=pal,
       cex=0.7, bg="transparent")
legend("topright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18,
       cex=0.7, bg="transparent")

plotMDS(getM(mymSetSqFlt), pch = 0:18, top=1000, gene.selection="common",
        col=pal[factor(mytargets$Sample_Group)], dim=c(3,4))
legend("topleft", legend=levels(factor(mytargets$Sample_Group)), text.col=pal,
       cex=0.7, bg="transparent")
legend("topright", legend=levels(factor(mytargets$Sample_Name)), pch = 0:18,
       cex=0.7, bg="transparent")
dev.off()
setwd(old.dir)


# Statistical Analysis ----------------------------------------------------
# M-values have nicer statistical properties, thus better for use in statistical
# analysis of methylation data

# calculate M-values for statistical analysis
mymVals <- getM(mymSetSqFlt)
mybVals <- getBeta(mymSetSqFlt)

# Plot density plots beta-Values & M-Values
setwd("Output")
pdf("density_plots_m+beta.pdf", paper = "a4r")
par(mfrow=c(1,2))
densityPlot(mybVals, sampGroups=mytargets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(mytargets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mymVals, sampGroups=mytargets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(mytargets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
dev.off()
setwd(old.dir)

# Probe-wise differential methylation analysis as an example
# limma User’s Guide extensively covers other types of designs commonly used


# Graveyard ---------------------------------------------------------------
## The R plot variant of the quality control ggplots
#pdf("mean_detection_p-values.pdf", paper = "a4")
##Colorblind-friendly Palette:
#pal <- brewer.pal(8,"Dark2")
# #fcol, mfrow: A vector of the form c(nr, nc). Subsequent figures will be drawn in
# #an nr-by-nc array on the device by columns (mfcol), or rows (mfrow),
# #respectively.
# par(mfrow=c(1,2))
# #Creates a bar plot with vertical (senkrecht zur Achse, clas=2) bars.
# barplot(colMeans(mydetP),
#         col=pal[factor(mytargets$Sample_Name)],
#         las=2,
#         cex.names=0.8,
#         ylab="Mean detection p-values")
# abline(h=0.05,col="red")
# legend("bottomleft", inset = 0.02, legend=levels(factor(mytargets$Sample_Name)), fill=pal,
#        bg="white")
#
##Same Plot with different y-axis lenght
# barplot(colMeans(mydetP),
#         col=pal[factor(mytargets$Sample_Name)],
#         las=2,
#         cex.names=0.8,
#         ylim=c(0,0.002),
#         ylab="Mean detection p-values")
# abline(h=0.05,col="red")
# legend("topleft", inset = 0.02, legend=levels(factor(mytargets$Sample_Name)), fill=pal,
#        bg="white")
#Close the pdf
#dev.off()