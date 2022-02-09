#############################################
### XCMS MS2 processing and exporting script
##
## This script processes MS1 data with XCMS
## and also exports MS2 data to the MGF file.
##
## Note: All files must have MS2 data, otherwise it will fail!
##
#

## clean environment and console for easy tracking of the calculations (do not change)
rm(list=ls())
cat("\014")  



## Parameters for processing
## 
## Working directory (results will be saved here)
workingDirectory = "D:/PeakBot_Data/PHM_comparison/XCMS"

## Path to MSMS samples
dda_path <- "D:/PeakBot_Data/PHM_comparison/XCMS/mzXMLs"

## Boolean specifying if fast polarity switching is used (TRUE) or not (FALSE)
##   If FPS is not used, the positive and negative mode samples must be placed in the subfolders pos and neg in the dda_path
isFPS = FALSE

## Path to export results to (typcially same as working directory)
expPath <- workingDirectory

## Associate samples to groups (via subsets of the samples' names)
##    The first sub-match is used, others are ignored. 
##    e.g. SampleName: "AGs_1" --> AGs group
##                     "QC_" --> QC group
##    Note: Advanced user can use regex here for grepl
sampleNameToGroup <- matrix(c("AOH", "AOH",
                              "AME", "AME"), byrow=TRUE, ncol=2)

## Peak picking parameters (with centWave, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2639432/)
centWavePPM = 5
centWavePeakWidth = c(2, 5)  # seconds
centWavePreFilter = c(1, 1E5)  # times x intensity [counts]
centWaveSNRThresh = 5
centWaveNoise = 5E4

## Feature properties (anything outside these parameters will be removed after peak picking)
startRT = 120 # seconds
endRT = 750 # seconds
SNRCutoff = 5
peakWidthCutoff = c(2, 35)

## Annotation ppm for CAMERA (higher than the ppm for Centwave as some isotopes have a strong mass shift)
CameraAnnotationPPM = 25




























## Processing script
setwd(workingDirectory)

## load necessary packages and define some custom functions
library(xcms)
library(CAMERA)
library(ggplot2)
source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")
ticStart__ptm<<-0
tic <- function(){ticStart__ptm<<-proc.time(); return(ticStart__ptm[3])}
toc <- function(unit="min"){return((proc.time()[3]-ticStart__ptm[3])/list(min=60, sec=1, hours=60*60, days=24*60*60)[unit][[1]])}
tocP <- function(unit="min"){print(sprintf("Time elapsed since last tic: %.1f %s", toc(), unit))}
exportAsFeatureXML <- function(peaks, fileTo, cid = "Num", crt = "rt", cmz = "mz", crtmin = "rtmin", crtmax = "rtmax", cmzmin = "mzmin", cmzmax = "mzmax", cquality = "sn"){
  
  lines = c()
  lines = c(lines, sprintf('<?xml version="1.0" encoding="ISO-8859-1"?>'))
  lines = c(lines, sprintf('  <featureMap version="1.4" id="fm_16311276685788915066" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/FeatureXML_1_4.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'))
  lines = c(lines, sprintf('    <dataProcessing completion_time="2011-02-09T11:23:14">'))
  lines = c(lines, sprintf('      <software name="XCMS" version="3.4.4" />'))
  lines = c(lines, sprintf('    </dataProcessing>'))
  lines = c(lines, sprintf('    <featureList count="%d">', nrow(peaks)))
  fileConn<-file(fileTo, open="w")
  writeLines(lines, fileConn)
  close(fileConn)
  
  for(rowi in 1:nrow(peaks)){
    
    id = ifelse(cid %in% colnames(peaks), peaks[rowi, cid], sprintf("Row%s", rownames(peaks)[rowi]))
    rt = peaks[rowi, crt]
    mz = peaks[rowi, cmz]
    rtmin = ifelse(crtmin %in% colnames(peaks), peaks[rowi, crtmin], rt - 5.)
    rtmax = ifelse(crtmax %in% colnames(peaks), peaks[rowi, crtmax], rt + 5.)
    mzmin = ifelse(cmzmin %in% colnames(peaks), peaks[rowi, cmzmin], mz * (1-5/1E6))
    mzmax = ifelse(cmzmax %in% colnames(peaks), peaks[rowi, cmzmax], mz * (1+5/1E6))
    quality = ifelse(cquality %in% colnames(peaks), peaks[rowi, cquality], -1)
    
    lines = c()
    lines = c(lines, sprintf('<feature id="%s">', id))
    lines = c(lines, sprintf('  <position dim="0">%f</position>', rt))
    lines = c(lines, sprintf('  <position dim="1">%f</position>', mz))
    lines = c(lines, sprintf('  <intensity>%f</intensity>', 1.))
    lines = c(lines, sprintf('  <quality dim="0">0</quality>'))
    lines = c(lines, sprintf('  <quality dim="1">0</quality>'))
    lines = c(lines, sprintf('  <overallquality>%s</overallquality>', quality))
    lines = c(lines, sprintf('  <charge>1</charge>'))
    lines = c(lines, sprintf('  <convexhull nr="0">'))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmin, mzmin))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmin, mzmax))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmax, mzmin))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmax, mzmax))
    lines = c(lines, sprintf('  </convexhull>'))
    lines = c(lines, sprintf('</feature>'))
    fileConn<-file(fileTo, open="a")
    writeLines(lines, fileConn)
    close(fileConn)
  }
  
  lines = c()
  lines = c(lines, sprintf('    </featureList>'))
  lines = c(lines, sprintf('  </featureMap>'))
  fileConn<-file(fileTo, open="a")
  writeLines(lines, fileConn)
  close(fileConn)
}

default <- registered()
multicoreParam <- SerialParam()
register(multicoreParam, default = TRUE)
if(FALSE){
  for (param in rev(default))
    register(param)
}







## Process the two polarities separately (if fast-polarity-switching was employed)
for(pol in c("pos")){  
  files = dir(dda_path, pattern=paste0("*", ".mzXML$"), full.names = TRUE, recursive = TRUE)
  
  ## Create a phenodata data.frame
  grp = sapply(files, function(x){
    for(rowi in 1:nrow(sampleNameToGroup)){
      if(grepl(sampleNameToGroup[rowi,1], x)){
        return(sampleNameToGroup[rowi, 2])
      }
    }
    return("default")
  })
  pd <- data.frame(sample_name = basename(files),
                   sample_group = grp,
                   stringsAsFactors = FALSE)
  
  
  ## Load raw data
  rawData <- readMSData(files, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")
  save.image(file=paste0(expPath, "/0_loaded_data.Rimage"))
  #load(paste0(expPath, "/0_loaded_data.Rimage"))
  
  
  ## Filter for polarity
  cat(paste0("\nProcessing polarity ", pol, "\n"))
  polInd = ifelse(pol == "pos", 1, 0)
  dda_data_pol <- filterPolarity(rawData, polInd)
  dda_data_pol <- removePeaks(dda_data_pol, 3E4)
  save.image(file=paste0(expPath, "/1_filteredData_", pol, ".Rimage"))
  #load(paste0(expPath, "/1_filteredData_", pol, ".Rimage"))
  
  
  ## Show overview of MS and MSMS spectra
  cat("\n\nOverview of number of MS and MSMS spectra\n")
  print(table(msLevel(dda_data_pol)))
  
  
  ## Detect chromatographic peaks
  cat("\n\nDetecting chromatographic peaks\n")
  cwp <- CentWaveParam(snthresh = centWaveSNRThresh, noise = centWaveNoise, 
                       ppm = centWavePPM, mzdiff = -0.0025, 
                       peakwidth = centWavePeakWidth, prefilter = centWavePreFilter, 
                       firstBaselineCheck = FALSE,
                       extendLengthMSW = TRUE,
                       mzCenterFun = "wMeanApex3", verboseColumns = TRUE)
  peaks <- findChromPeaks(dda_data_pol, param = cwp)
  save.image(file=paste0(expPath, "/2_DetectedChromPeaks_", pol, ".Rimage"))
  #load(paste0(expPath, "/2_DetectedChromPeaks_", pol, ".Rimage"))
  
  
  cat("\n\nMerging separated peaks\n")
  mpp <- MergeNeighboringPeaksParam(expandRt = 4)
  peaks <- refineChromPeaks(peaks, mpp)
  save.image(file=paste0(expPath, "/3_RefinedChromPeaks_", pol, ".Rimage"))
  #load(paste0(expPath, "/3_RefinedChromPeaks_", pol, ".Rimage"))
  
  
  cat("\n\nRemoving too early/late peaks and applying SNR filter\n")
  temp = chromPeaks(peaks)
  use = temp[,"rt"]>startRT & temp[,"rt"]<endRT
  cat(paste0("Removing ", sum(!use), " features outside of specified Rt boundaries\n"))
  temp = temp[use,]
  
  chromPeaks(peaks) = temp
  save.image(file=paste0(expPath, "/4_FilteredChromPeaks_", pol, ".Rimage"))
  #load(paste0(expPath, "/4_FilteredChromPeaks_", pol, ".Rimage"))
  
  
  cat("\n\nExporting single sample results\n")
  temp = chromPeaks(peaks)
  for(si in unique(temp[,"sample"])){
    ta = temp[temp[,"sample"]==si,]
    if(!is.null(nrow(ta))){
      cat(paste0("There are ", nrow(ta), " chromatographic peaks in the sample ", files[si], " in polarity ", pol, "\n"))
      exportAsFeatureXML(ta, paste0(expPath, "/", peaks@phenoData@data$sample_name[si], "_", pol, ".featureML"), cid = "Num", crt = "rt", cmz = "mz", crtmin = "rtmin", crtmax = "rtmax", cmzmin = "mzmin", cmzmax = "mzmax")
      write.table(ta, paste0(expPath, "/", peaks@phenoData@data$sample_name[si], "_", pol, ".tsv"), sep="\t", na="", quote=FALSE)
    }else{
      cat(paste0("No chromatographic peaks detected in the sample ", files[si], " in polarity ", pol, "\n"))
    }
  }
  
  
  ## Group results
  cat("\n\nGrouping results\n")
  pdp <- PeakDensityParam(sampleGroups = rawData$sample_group,
                          binSize = 0.005, 
                          minFraction = 0.4, bw = 2)
  grouped <- groupChromPeaks(peaks, param = pdp)
  save.image(file=paste0(expPath, "/5_GroupedResults_", pol, ".Rimage"))
  #load(paste0(expPath, "/5_GroupedResults_", pol, ".Rimage"))
  
  ## Fill missing peaks (semi-targeted re-integration)
  cat("\n\nFilling peaks\n")
  medWidth <- median(chromPeaks(grouped)[, "rtmax"] - chromPeaks(grouped)[, "rtmin"])
  filled <- fillChromPeaks(grouped, param = FillChromPeaksParam(fixedRt = medWidth))
  save.image(file=paste0(expPath, "/6_FilledNAs_", pol, ".Rimage"))
  #load(paste0(expPath, "/6_FilledNAs_", pol, ".Rimage"))
  
  
  ## Annotate features with CAMERA
  cat("\n\nAnnotating features with CAMERA\n")
  xset <- as(filterMsLevel(filled, msLevel = 1L), "xcmsSet")
  sampclass(xset) <- pd$sample_group
  xsa <- xsAnnotate(xset, polarity = ifelse(pol=="pos", "positive", "negative"))
  xsa <- groupFWHM(xsa, sigma = 6, perfwhm = 1)
  #xsa <- groupCorr(xsa, cor_eic_th = 0.6, pval = 0.05, graphMethod = "hcs",
  #                 calcCiS = TRUE, calcCaS = TRUE, calcIso = FALSE)
  xsa <- findIsotopes(xsa, maxcharge = 2, maxiso = 3, minfrac = 0.5,
                      ppm = CameraAnnotationPPM, intval = "maxo", filter = FALSE)
  xsa <- findAdducts(xsa, polarity = "positive", 
                     max_peaks = 100, multiplier = 3, ppm = CameraAnnotationPPM)
  camera_feature_ann <- getFeatureAnnotations(xsa)
  save.image(file=paste0(expPath, "/8_AnnoatedCAMERA_", pol, ".Rimage"))
  #load(paste0(expPath, "/8_AnnoatedCAMERA_", pol, ".Rimage"))
  
  
  ## Export peak abundance table
  cat("\n\nExporting feautre table\n")
  featuresDef <- featureDefinitions(filled)
  featuresIntensities <- featureValues(filled, value = "into")
  ## generate data table
  dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
  dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]
  dataTable <- cbind(dataTable, camera_feature_ann)
  write.table(dataTable, paste0(expPath, "/", "__peaks_", pol, "_quant.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  cat("\n\n\n")
}

message(sprintf("Calculating took %.2f minutes\n\n\n\n\n", toc()))  ## print how long the calculations needed
