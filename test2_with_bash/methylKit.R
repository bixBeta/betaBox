#!/usr/bin/env Rscript

args <-  commandArgs(trailingOnly = T)
outputPrefix <- args[1]
bedFile <- args[2]
cpgBed <- args[3]


#bedFile <- "hg38.bed"
#cpgBed <- "hg38_CPGs.bed"
#outputPrefix <- "seqtkSamples"
library(methylKit)

##################################################################
##################################################################
## Import coverage files
##################################################################
##################################################################

phenoData <- as.data.frame(read.csv(file = "inputFile.csv", header = T, sep = ","))

sampleNames <- as.character(phenoData$sample.name)
ID <- as.character(phenoData$sample.id)
treatment <- phenoData$treatment

file.list= as.list(sampleNames)
ID = as.list(ID)


myobj = methRead(file.list,
                 sample.id = ID,
                 pipeline = "bismarkCoverage",
                 assembly = "hg38",
                 treatment = treatment, mincov = 1) # control 0

##################################################################
##################################################################
## QC, basic statistics
##################################################################
##################################################################
sink(file = "MethylationSummary.log")
for (i in 1:length(myobj)) {
  print(paste0(phenoData$sample.name[i], "_MethylationSummary:"))
  getMethylationStats(myobj[[i]], plot = F, both.strands = F)
}
sink()

##################################################################
##################################################################
## Plot HISTOGRAMS
##################################################################
##################################################################
library("graphics")
system("mkdir MethylationStatsPlots")
for (i in 1:length(myobj)) {
  jpeg(filename = paste0(phenoData$sample.name[i],".methylationStats.jpg"), width = 1080, height = 1080, res = 155)
  getMethylationStats(myobj[[i]], plot = T, both.strands = F)
  dev.off()
}
system("mv *.methylationStats.jpg MethylationStatsPlots")

system("mkdir CoverageStatsPlots")
for (i in 1:length(myobj)) {
  jpeg(filename = paste0(phenoData$sample.name[i],".coverageStats.jpg"), width = 1080, height = 1080, res = 155)
  getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
  dev.off()
}
system("mv *.coverageStats.jpg CoverageStatsPlots")

## filter out low quality reads (Optional)
##################################################################
# filtered.myobj <- filterByCoverage(myobj, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
##################################################################

##################################################################
##################################################################
##  DIAGNOSTIC PLOTS
##################################################################
##################################################################

# sample corelation
meth <- unite(myobj, destrand = FALSE)
head(meth)

system("mkdir DiagnosticPlots")
jpeg(filename = paste0(outputPrefix,".SampleCorelationPlot" ,".jpg"), width = 1080, height = 1080, res = 155)
getCorrelation(meth, plot = T)
dev.off()
# cluster based on methyl profile
jpeg(filename = paste0(outputPrefix,".ClusterPlot" ,".jpg"), width = 1080, height = 1080, res = 155)
clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
dev.off()
hc <- clusterSamples(meth, dist = "correlation", method = "ward", plot = FALSE)
hc
# PCA
jpeg(filename = paste0(outputPrefix, ".PCA.plot" ,".jpg"), width = 1080, height = 1080, res = 155)
PCASamples(meth, screeplot = TRUE)
dev.off()
system("mv *.jpg DiagnosticPlots")

##################################################################
##################################################################
##  Differential Methylation 
##################################################################
##################################################################

myDiff <- calculateDiffMeth(meth, mc.cores = 4)
write.csv(as(myDiff,"GRanges"), paste0(outputPrefix, "_DiffCalls.csv"), quote = F)

## CUSTOM 
# get all differentially methylated bases
# myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

##################################################################
##################################################################
##  Annotating DMB's 
##################################################################
##################################################################
library(genomation)
gene.obj=readTranscriptFeatures(bedFile)


#
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
#

diffAnn=annotateWithGeneParts(as(myDiff,"GRanges"),gene.obj)
diffAnn
diffAnn_TSS <- getAssociationWithTSS(diffAnn)
head(diffAnn_TSS)
write.csv(diffAnn_TSS, paste0(outputPrefix, "_AssociationsWithTSS.csv"), quote = F)

genomation::getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

jpeg(filename = paste0(outputPrefix, ".methylationAnnotation" ,".jpg"), width = 1080, height = 1080, res = 155)
genomation::plotTargetAnnotation(diffAnn,precedence=TRUE, main="differential methylation annotation")
dev.off()


# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(cpgBed, feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")
diffCpGann

jpeg(filename = paste0(outputPrefix, ".CPG_Annotation" ,".jpg"), width = 1080, height = 1080, res = 155)
genomation::plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="differential methylation annotation")
dev.off()
