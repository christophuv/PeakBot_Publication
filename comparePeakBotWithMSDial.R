
rm(list=ls())


exportAsFeatureXML<-function(peaks, fileTo, cid = "Num", crt = "rt", cmz = "mz", crtmin = "rtmin", crtmax = "rtmax", cmzmin = "mzmin", cmzmax = "mzmax", rtConv2Sec = TRUE){
  
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
    
    id = ifelse(cid %in% colnames(peaks), peaks[rowi, cid], sprintf("Row%d", rowi))
    rt = ifelse(rtConv2Sec, peaks[rowi, crt]*60., peaks[rowi, crt])
    mz = peaks[rowi, cmz]
    rtmin = ifelse(crtmin %in% colnames(peaks), ifelse(rtConv2Sec, peaks[rowi, crtmin]*60., peaks[rowi, crtmin]), rt - 5.)
    rtmax = ifelse(crtmax %in% colnames(peaks), ifelse(rtConv2Sec, peaks[rowi, crtmax]*60., peaks[rowi, crtmax]), rt + 5.)
    mzmin = ifelse(cmzmin %in% colnames(peaks), peaks[rowi, cmzmin], mz * (1-5/1E6))
    mzmax = ifelse(cmzmax %in% colnames(peaks), peaks[rowi, cmzmax], mz * (1+5/1E6))
    
    
    lines = c()
    lines = c(lines, sprintf('<feature id="%s">', id))
    lines = c(lines, sprintf('  <position dim="0">%f</position>', rt))
    lines = c(lines, sprintf('  <position dim="1">%f</position>', mz))
    lines = c(lines, sprintf('  <intensity>%f</intensity>', 1.))
    lines = c(lines, sprintf('  <quality dim="0">0</quality>'))
    lines = c(lines, sprintf('  <quality dim="1">0</quality>'))
    lines = c(lines, sprintf('  <overallquality>0</overallquality>'))
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




library(ggplot2)
library(grid)
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

res = matrix(0, ncol=6, nrow=0)
colnames(res)=c("only msdial", "only msdial %", "both", "both %", "only PeakBot", "only PeakBot %")

set.seed(2021)
for(fil in c("HT_SOL1_LYS_010_pos","HT_SOL1_SUP_025_pos","HT_SOL2_LYS_014_pos","HT_SOL2_SUP_029_pos","HT_SOL3_LYS_018_pos","HT_SOL3_SUP_033_pos")){
  cat("Processing dataset", fil, "\n")
  
  peaks <- read.delim2(paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_proc.txt"))
  exportAsFeatureXML(peaks, paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, ".msdial.featureML"), "PeakID", "RT..min.", "Precursor.m.z", "RT.left.min.", "RT.right..min.")
  peaks <- cbind(peaks, foundPB = 0, num = 0)
  
  peaksPB   <- read.table(paste0("/home/cbueschl/Documents/_backups/PeakBot/peakbot_example/Data/MTBLS1358/", fil, "_positivePeakBot.tsv"), sep="\t", quote="", header=TRUE)
  peaksPB   <- cbind(peaksPB, foundXCMS = 0, num = 0)
  
  dat = data.frame(rtx=c(), rtp=c(), mzx=c(), mzp=c(), rtdiff=c(), mzdiff=c())
  
  cur = 1
  for(rowi in 1:nrow(peaks)){
    xrt = peaks[rowi, "RT..min."]*60.
    xmz = peaks[rowi, "Precursor.m.z"]
    
    for(rowj in 1:nrow(peaksPB)){
      prt = peaksPB$RT[rowj]
      pmz = peaksPB$MZ[rowj]
      
      if(abs(xrt-prt)<=4 && abs(xmz-pmz)/xmz*1E6<10){
        peaks$foundPB[rowi]=paste0("A", cur)
        peaksPB$foundXCMS[rowj]=paste0("A", cur)
        peaks$num[rowi]=paste0("A", cur)
        peaksPB$num[rowj]=paste0("A", cur)
        
        dat = rbind(dat, data.frame(rtx=xrt, rtp=prt, mzx=xmz, mzp=pmz, rtdiff=xrt-prt, mzdiff=(xmz-pmz)/xmz*1E6))
        
        cur = cur + 1
      }
    }
    if(peaks$num[rowi]==0){
      peaks$num[rowi] = paste0("A", cur)
      peaks$foundPB[rowi] = paste0("N", cur)
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
  
  l = list(msDial = peaks$num, PeakBot = peaksPB$num)
  anums = unique(setdiff(l$msDial, l$PeakBot))
  abnums = unique(intersect(l$msDial, l$PeakBot))
  bnums = unique(setdiff(l$PeakBot, l$msDial))
  a = length(anums)
  ab = length(abnums)
  b = length(bnums)
  res = rbind(res, c(a, round(a/(a+ab+b)*100,2), ab, round(ab/(a+ab+b)*100,2), b, round(b/(a+ab+b)*100,2)))
  rownames(res)[nrow(res)] = fil
  
  grid.newpage()
  grid.draw(venn.diagram(x = l, filename=NULL))
  
  exportAsFeatureXML(peaks[grepl("A", peaks$foundPB),], paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_MSDial_Peakbot.featureML"), "PeakID", "RT..min.", "Precursor.m.z", "RT.left.min.", "RT.right..min.")
  write.table(peaks[grepl("A", peaks$foundPB),], paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_MSDial_Peakbot.tsv"), col.names=NA)
  
  temp = peaksPB[grepl("N", peaksPB$foundXCMS),]
  samp = sample(1:nrow(temp), 100)
  exportAsFeatureXML(temp[samp,], paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_PeakbotONLY_randomlySelected.featureML"), "Num", "RT", "MZ", "RtStart", "RtEnd", "MzStart", "MzEnd", rtConv2Sec=FALSE)
  write.table(temp[samp,], paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_PeakbotONLY_randomlySelected.tsv"), col.names=NA)
  exportAsFeatureXML(temp, paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_PeakbotONLY.featureML"), "Num", "RT", "MZ", "RtStart", "RtEnd", "MzStart", "MzEnd", rtConv2Sec=FALSE)
  write.table(temp, paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_PeakbotONLY.tsv"), col.names=NA)
  
  temp = peaks[grepl("N", peaks$foundPB),]
  samp = sample(1:nrow(temp), 100)
  exportAsFeatureXML(temp[samp,], paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_MSDialONLY_randomlySelected.featureML"), "PeakID", "RT..min.", "Precursor.m.z", "RT.left.min.", "RT.right..min.")
  write.table(temp[samp,], paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_MSDialONLY_randomlySelected.tsv"), col.names=NA)
  exportAsFeatureXML(temp, paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_MSDialONLY.featureML"), "PeakID", "RT..min.", "Precursor.m.z", "RT.left.min.", "RT.right..min.")
  write.table(temp, paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_MSDialONLY.tsv"), col.names=NA)
  
  if(fil=="HT_SOL1_LYS_010_pos"){
    temp = read.delim2(paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_MSDialONLY_randomlySelected_annotated.tsv"))
    print("MS Dial annotated features")
    print(table(temp[,"Comment"]))
    
    temp = read.delim2(paste0("/home/cbueschl/Documents/Projects/PeakBot/Data/MTBLS1358/", fil, "_PeakbotONLY_randomlySelected_annotated.tsv"))
    print("PeakBot annotated features")
    print(table(temp[,"comment"]))
  }
}
print(res)
