## Part 1 - detection of chrom. peaks with XCMS

rm(list=ls())
cat("\014")  
setwd("/home/users/cbueschl/PeakBot/peakbot_example/Data/PHM/PositiveCentroidMode")

# xset@peaks[xset@peaks[,"mz"]<175.1230 & xset@peaks[,"mz"]>175.1100 &
#            xset@peaks[,"rt"]<230 & xset@peaks[,"rt"]>210,]


#############################
### load required packages
library(snow)
library(xcms)
library(CAMERA)

source("/home/users/cbueschl/LV_DataScience/impl/ext__exportAsFeatureML.R")

#############################
### implement some common functions
tic<-function(){ticStart__ptm<<-proc.time(); return(ticStart__ptm[3])}
toc<-function(unit="min"){return((proc.time()[3]-ticStart__ptm[3])/list(min=60, sec=1, hours=60*60, days=24*60*60)[unit][[1]])}
tocP<-function(unit="min"){print(sprintf("Time elapsed since last tic: %.1f %s", toc(), unit))}





cpus=32
totTime=0
overallStartTime=tic()

if(FALSE){
  #############################
  ### process LC-HRMS dataset with XCMS
  ### 1. find features (i.e. chromatographic peaks) in each sample file (mzXML)
  {
    tic()                                             ## start the timer
    xset<-xcmsSet(method="centWave",                  ## XCMS processing with the centWave algorithm
                  files=".",                          ##      files to process
                  peakwidth=c(2, 5),                  ##      set the expected chromatographic peak width (baseline)
                  ppm=5,                              ##      set the expected/allowed ppm deviation
                  snthr=5,                            ##      minimum signal-to-noise ratio for chromatographic peaks in EICs
                  
                  integrate=1,                        ##      specify that the peak borders shall be found via the wavelet transformed data 
                  mzCenterFun="wMeanApex3",           ##      this specifies how the features' mz values should be calculated
                  mzdiff=0.0025,                      ##      minimum m/z difference for RT-overlapping chromatographic peaks
                  prefilter=c(1, 1E5),                ##      filter criterion for initial ROI generation
                  noise=5000,                         ##      noise intensity level
                  fitgauss=TRUE,                      ##      specifies that a Gaussian distribution shall be fitted to each feature
                  
                  verbose.columns=TRUE,               ##      return additional peak meta-data in the results
                  nSlaves=cpus)                       ##      Deprecated: specify that several CPU cores of the PC shall be used for the calculations
    
    xset@peaks = xset@peaks[100 <= xset@peaks[,"rt"] & xset@peaks[,"rt"] <= 750,]
    xset@peaks = xset@peaks[xset@peaks[,"maxo"] >= 1E5,]
  
    for(i in 1:length(xset@filepaths)){
      fi = xset@filepaths[i]
      peaks = xset@peaks[xset@peaks[,"sample"]==i,]
      hist(peaks[,"rtmax"]-peaks[,"rtmin"], breaks=c(-Inf, seq(0,250,0.5), Inf), xlim=c(0,20), main="Peak width (max-min)")
      cat("File", xset@filepaths[i], "has", nrow(peaks), "features\n")
      exportAsFeatureXML(peaks, paste0(gsub(".mzML", "", fi), "_xcms.featureML"))
      write.table(peaks, paste0(gsub(".mzML", "", fi), "_xcms.tsv"), sep="\t", na="", quote=FALSE)
    }
   
    totTime=totTime+toc()                             ## add total time
    message(sprintf("Calculating the xcms set object took %.2f minutes\n\n\n\n\n", toc()))  ## print how long the calculations needed
  }
  
  
  #############################
  ### process LC-HRMS dataset with XCMS
  ### 2. group features from all samples together and align their retention times
  {
    tic()
    xset.aligned <-  group(xset          , method="density", minfrac=0.33, bw=2, mzwid=0.03)
    xset.aligned <- retcor(xset.aligned  , method="loess",   span=.5, plottype="none")
    xset.aligned <-  group(xset.aligned  , method="density", minfrac=0.33, bw=2, mzwid=0.03)
    xset.aligned <-fillPeaks(xset.aligned, method="chrom", nSlaves=cpus)
    
    totTime=totTime+toc()                      ## add total time
    message(sprintf("Aligning the chromatograms took %.2f minutes\n\n\n\n\n", toc()))
  }
  
  
  #############################
  ### process LC-HRMS dataset with XCMS
  ### 3. annotate the detected features with isotopes, adducts and group them
  {
    tic()
    xset.anno <-   xsAnnotate(xset.aligned, nSlaves=cpus, polarity="positive")  ## create xsAnnotate object
    xset.anno <-    groupFWHM(xset.anno)                         ## group according to retention time
    xset.anno <- findIsotopes(xset.anno, maxcharge=2, ppm=25)    ## search for isotopes
    xset.anno <-    groupCorr(xset.anno)                         ## check grouping with EIC peaks
    xset.anno <-  findAdducts(xset.anno, polarity="positive")    ## search for possible adduct combinations
    
    totTime=totTime+toc()                  ## add total time
    message(sprintf("Annotating the results took %.2f minutes\n\n\n\n\n", toc()))
  }
  
  
  #############################
  ### process LC-HRMS dataset with XCMS
  ### 4. extract and save the results to a TSV file and also as a featureXML file so they can easily be visualized with ToppVIEW
  peaks=getPeaklist(xset.anno)
  peaks=cbind(Num=paste0("F", 1:nrow(peaks)), peaks)
  write.table(peaks, file="./peaks.tsv", sep="\t", row.names=FALSE, quote=FALSE, na="")
  
  exportAsFeatureXML(peaks, "peaks.featureML")
  
  message(nrow(peaks), " features have been detected")
  message(sprintf("FINISHED XCMS (%.1f minutes in total, %d cpus used)", totTime, cpus))
}








## Part 2 - comparison of detected chrom. peaks with PeakBot and XCMS

rm(list=ls())
setwd("~/PeakBot/peakbot_example/Data/PHM")

library(ggplot2)
library(grid)
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

source("~/LV_DataScience/impl/ext__exportAsFeatureML.R")

files = c("05_EB3388_AOH_p_0",  "06_EB3389_AOH_p_10", "07_EB3390_AOH_p_20", "08_EB3391_AOH_p_60",
          "16_EB3392_AME_p_0",  "17_EB3393_AME_p_10", "18_EB3394_AME_p_20", "19_EB3395_AME_p_60")

res = matrix(0, ncol=6, nrow=0)
colnames(res)=c("only xcms", "only xcms %", "both", "both %", "only PeakBot", "only PeakBot %")

for(fi in files){
  
  peaksXCMS <- read.table(paste0("./PositiveCentroidMode/", fi, "_xcms.tsv"), sep="\t", quote="")
  peaksXCMS <- cbind(peaksXCMS, foundPB = 0, num = 0)
  #plot(peaksXCMS$rt, peaksXCMS$mz, xlim=c(100, 750), ylim=c(100, 1000))
  
  peaksPB   <- read.table(paste0("./", fi, "_positivePeakBot.tsv"), sep="\t", quote="", header=TRUE)
  peaksPB   <- cbind(peaksPB, foundXCMS = 0, num = 0)
  #plot(peaksPB$RT, peaksPB$MZ, xlim=c(100, 750), ylim=c(100, 1000))
  
  dat = data.frame(rtx=c(), rtp=c(), mzx=c(), mzp=c(), rtdiff=c(), mzdiff=c())
  
  cur = 1
  for(rowi in 1:nrow(peaksXCMS)){
    xrt = peaksXCMS$rt[rowi]
    xmz = peaksXCMS$mz[rowi]
    
    for(rowj in 1:nrow(peaksPB)){
      prt = peaksPB$RT[rowj]
      pmz = peaksPB$MZ[rowj]
      
      if(abs(xrt-prt)<2 && abs(xmz-pmz)/xmz*1E6<10){
        peaksXCMS$foundPB[rowi]=paste0("A", cur)
        peaksPB$foundXCMS[rowj]=paste0("A", cur)
        peaksXCMS$num[rowi]=paste0("A", cur)
        peaksPB$num[rowj]=paste0("A", cur)
        
        dat = rbind(dat, data.frame(rtx=xrt, rtp=prt, mzx=xmz, mzp=pmz, rtdiff=xrt-prt, mzdiff=(xmz-pmz)/xmz*1E6))
        
        cur = cur + 1
      }
    }
    if(peaksXCMS$num[rowi]==0){
      peaksXCMS$num[rowi] = paste0("A", cur)
      peaksXCMS$foundPB[rowi] = paste0("N", cur)
      cur = cur + 1
    }
  }
  for(rowj in 1:nrow(peaksPB)){
    if(peaksPB$num[rowj]==0){
      peaksPB$num[rowj] = paste0("A", cur)
      peaksPB$foundXCMS[rowj] = paste0("N", cur)
      cur = cur + 1
    }
  }
  
  l = list(xcms = peaksXCMS$num, PeakBot = peaksPB$num)
  anums = unique(setdiff(l$xcms, l$PeakBot))
  abnums = unique(intersect(l$xcms, l$PeakBot))
  bnums = unique(setdiff(l$PeakBot, l$xcms))
  a = length(anums)
  ab = length(abnums)
  b = length(bnums)
  res = rbind(res, c(a, round(a/(a+ab+b)*100,2), ab, round(ab/(a+ab+b)*100,2), b, round(b/(a+ab+b)*100,2)))
  rownames(res)[nrow(res)] = fi
  
  grid.newpage()
  grid.draw(venn.diagram(x = l, filename=NULL))
  
  exportAsFeatureXML(peaksXCMS[peaksXCMS$foundPB>0,], gsub("_xcms.tsv", "_both.featureML", fi))
  write.table(peaksXCMS[peaksXCMS$foundPB>0,], gsub("_xcms.tsv", "_both.tsv", fi), col.names=NA)
  
  if(fi=="08_EB3391_AOH_p_60"){
    peaksPB = cbind(peaksPB, annotation = "")
    anno = read.table("/home/users/cbueschl/PeakBot/peakbot_example/Data/PHM/PositiveCentroidMode/__Annotated__08_EB3391_AOH_p_60_PBOnly.tsv", sep="\t", header=TRUE)
    for(rowi in 1:nrow(peaksPB)){
      if(grepl("N", peaksPB$foundXCMS[rowi])){
        rt = peaksPB[rowi,"RT"]
        mz = peaksPB[rowi,"MZ"]
        for(rowj in 1:nrow(anno)){
          if(abs(anno[rowj,"RT"]-rt)<=2 && abs(anno[rowj,"MZ"]-mz)*1E6/anno[rowj,"MZ"]<=10){
            peaksPB[rowi,"annotation"] = anno[rowj, "Type"]
          }
        }
      }
    }
    print(paste0("Features only detected by PeakBot in sample ", fi, ": ", length(peaksPB[grepl("N", peaksPB$foundXCMS),"annotation"])))
    print(table(peaksPB[grepl("N", peaksPB$foundXCMS),"annotation"]))
  }
  exportAsFeatureXML(peaksPB[grepl("N", peaksPB$foundXCMS),], gsub("_xcms.tsv", "_PBOnly.featureML", fi), "RT", "MZ", "RtStart", "RtEnd", "MzStart", "MzEnd")
  write.table(peaksPB[grepl("N", peaksPB$foundXCMS),], gsub("_xcms.tsv", "_PBOnly.tsv", fi), col.names=NA)
  
  
  if(fi=="08_EB3391_AOH_p_60"){
    peaksXCMS = cbind(peaksXCMS, annotation = "")
    anno = read.table("/home/users/cbueschl/PeakBot/peakbot_example/Data/PHM/PositiveCentroidMode/__Annotated__08_EB3391_AOH_p_60_xcmsOnly.tsv", sep="\t", header=TRUE)
    for(rowi in 1:nrow(peaksXCMS)){
      rt = peaksXCMS[rowi,"rt"]
      mz = peaksXCMS[rowi,"mz"]
      for(rowj in 1:nrow(anno)){
        if(abs(anno[rowj,"rt"]-rt)<=2 && abs(anno[rowj,"mz"]-mz)*1E6/mz<=5){
          peaksXCMS[rowi,"annotation"] = anno[rowj, "Type"]
        }
      }
    }
    print(paste0("Features only detected by XCMS in sample ", fi, ": ", length(peaksXCMS[grepl("N", peaksXCMS$foundPB),"annotation"])))
    print(table(peaksXCMS[grepl("N", peaksXCMS$foundPB),"annotation"]))
  }
  exportAsFeatureXML(peaksXCMS[grepl("N", peaksXCMS$foundPB),], gsub("_xcms.tsv", "_xcmsOnly.featureML", fi))
  write.table(peaksXCMS[grepl("N", peaksXCMS$foundPB),], gsub("_xcms.tsv", "_xcmsOnly.tsv", fi), col.names=NA)
}
print(res)

